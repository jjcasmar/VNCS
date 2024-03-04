#include "NBESolver.h"

#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/OdeSolver.h>

#include <VNCS/Logger.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace Log
{
struct NBESolver {
    constexpr static const char *name = "NBESolver";
    static int createLogger()
    {
        std::shared_ptr<spdlog::logger> logger = spdlog::stdout_color_mt(Log::NBESolver::name);
        logger->set_level(spdlog::level::debug);
        return 0;
    }

    inline static const int logger = NBESolver::createLogger();
};
}  // namespace Log

namespace VNCS
{
using sofa::core::VecId;
using namespace sofa::defaulttype;
using namespace sofa::core::behavior;

NBESolver::NBESolver()
    : sofa::core::behavior::OdeSolver()
    , m_rayleighStiffness(0)
    , m_rayleighMass(0)
    , m_newtonThreshold(1e-5)
    , m_maxNewtonIterations(1)
    , m_lineSearchMaxIterations(1)
{
}

void NBESolver::init()
{
    sofa::core::behavior::OdeSolver::init();
}

void NBESolver::cleanup()
{
}

void NBESolver::solve(const sofa::core::ExecParams *params,
                      SReal dt,
                      sofa::core::MultiVecCoordId xResult,
                      sofa::core::MultiVecDerivId vResult)
{
    Logger<Log::NBESolver>::info("Newton solver");
    sofa::simulation::common::VectorOperations vop(params, this->getContext());
    sofa::simulation::common::MechanicalOperations mop(params, this->getContext());
    MultiVecCoord pos(&vop, sofa::core::VecCoordId::position());
    MultiVecDeriv vel(&vop, sofa::core::VecDerivId::velocity());
    MultiVecDeriv f(&vop, sofa::core::VecDerivId::force());
    MultiVecDeriv b(&vop);
    MultiVecCoord newPos(&vop, xResult);
    MultiVecDeriv newVel(&vop, vResult);

    b.realloc(&vop, false, true);

    // To store the values at the start of the timestep
    MultiVecCoord pos0(&vop);
    MultiVecDeriv vel0(&vop);
    pos0.realloc(&vop, false, true);
    vel0.realloc(&vop, false, true);

    // To store the previous values before testing new values in the line search
    MultiVecDeriv vel_i(&vop);
    vel_i.realloc(&vop, false, true);

    /// inform the constraint parameters about the position and velocity id
    mop.cparams.setX(xResult);
    mop.cparams.setV(vResult);

    // dx is no longer allocated by default (but it will be deleted automatically by the mechanical objects)
    MultiVecDeriv dx(&vop, sofa::core::VecDerivId::dx());
    dx.realloc(&vop, false, true);

    mop->setImplicit(true);  // this solver is implicit

    // Store pos and vel at the start of the solve
    pos0.eq(pos);
    vel0.eq(vel);

    // Take an explicit step to bootstrap Newton solver
    pos.eq(pos0, vel0, dt);
    mop.propagateXAndV(pos, vel);

    const auto rhsEvaluator = [this, &b, &f, &vel0, &vel, &mop, &dt]() {
        b.clear();                                 // b = 0
        f.clear();                                 // f = 0;
        f.eq(vel0);                                // f = v_0
        f.peq(vel, -(1.0 + dt * m_rayleighMass));  // f = (v_0 - (1+h*a) * v)
        mop.addMdx(b, f, -1.0);                    // b = M * (v_0 - (1+h*a) * v)
        f.clear();
        mop.computeForce(f);  // f = FORCES
        b.peq(f, dt);         // b = h*f + M * (v_0 - (1+h*a) * v)

        mop.projectResponse(b);
    };

    const auto energyEvaluator = [&b, &f, &vel0, &vel, &mop, &dt]() {
        b.clear();              // b = 0
        f.clear();              // f = 0;
        f.eq(vel, vel0, -1.0);  // f = (x - x*)
        f.teq(dt);
        mop.addMdx(b, f, 1.0);  // b = M * (x - x*)
        const auto kinetic = 1.0 / (2.0 * dt * dt) * f.dot(b);

        VNCS::Real kinetic_ = 0;
        VNCS::Real potential = 0;
        mop.computeEnergy(kinetic_, potential);

        return kinetic + potential;
    };

    // Compute the rhs after, as the errorEvaluator uses the b vector as placeholder!
    rhsEvaluator();
    VNCS::Real errorAtIteration = b.dot(b);

    int currentNewtonIteration = 0;

    Logger<Log::NBESolver>::info("Initial error = {}", errorAtIteration);
    bool lineSearchConverged = false;

    while ((currentNewtonIteration == 0) ||  // Forces at least one Newton iteration
           (errorAtIteration > m_newtonThreshold && currentNewtonIteration < m_maxNewtonIterations)) {
        sofa::core::behavior::MultiMatrix<sofa::simulation::common::MechanicalOperations> matrix(&mop);
        matrix = MechanicalMatrix(1.0 + dt * m_rayleighMass,  // mass component
                                  -dt,                        // damping component
                                  -dt * dt                    // stiffness component
        );

        dx.clear();
        mop->setDx(dx);
        matrix.solve(dx, b);

        VNCS::Real lineSearchStep = 2.0;
        int lineSearchCurrentIteration = 0;
        VNCS::Real error_i1 = 0.0;

        // Perform a line search on the direction of dx

        // Save the values at the beginning of the line search
        // This is needed in case we need to reduce the line search step!
        vel_i.eq(vel);

        lineSearchConverged = false;
        VNCS::Real initialEnergy = energyEvaluator();
        while (lineSearchCurrentIteration < m_lineSearchMaxIterations) {
            lineSearchStep *= 0.5;
            vel.eq(vel_i, dx, lineSearchStep);  // v_i1 = v_i + step * dx
            pos.eq(pos0, vel, dt);              // pos_i1 = pos_0 + h * v_i1

            // Propagate this temptative position and velocity
            mop.propagateXAndV(pos, vel);

            // error_i1 = b.dot(b);
            const auto energy = energyEvaluator();

            Logger<Log::NBESolver>::info(
                "a = {}. New energy = {}. Initial energy = {}", lineSearchStep, energy, initialEnergy);

            if (energy <= initialEnergy) {
                lineSearchConverged = true;
                break;
            }

            lineSearchCurrentIteration++;
        }

        if (lineSearchConverged) {
            // Prepare next Newton iteration
            rhsEvaluator();
            errorAtIteration = b.dot(b);
            Logger<Log::NBESolver>::info("Line search reduces energy. Current Newton error is = {}", errorAtIteration);
        } else {
            pos.eq(pos0);
            vel.eq(vel0);
            Logger<Log::NBESolver>::info("Line search failed. Cancelling Newton");
            break;
        }

        currentNewtonIteration++;
    }

    if (errorAtIteration > m_newtonThreshold) {
        // The Newton solver failed as it didn't reduce the error of  h * F = M (v - v0) enough
        pos.eq(pos0);
        vel.eq(vel0);
        Logger<Log::NBESolver>::info("Newton failed after {} iterations", currentNewtonIteration);
    } else {
        Logger<Log::NBESolver>::info("Finished Newton after {} iterations", currentNewtonIteration);
        Logger<Log::NBESolver>::info("Current error = {}", errorAtIteration);
    }

    newPos.eq(pos);
    newVel.eq(vel);
}

VNCS::Real NBESolver::rayleighStiffness() const
{
    return m_rayleighStiffness;
}

void NBESolver::setRayleighStiffness(VNCS::Real newRayleighStiffness)
{
    m_rayleighStiffness = newRayleighStiffness;
}

VNCS::Real NBESolver::rayleighMass() const
{
    return m_rayleighMass;
}

void NBESolver::setRayleighMass(VNCS::Real newRayleighMass)
{
    m_rayleighMass = newRayleighMass;
}

VNCS::Real NBESolver::newtonThreshold() const
{
    return m_newtonThreshold;
}

void NBESolver::setNewtonThreshold(VNCS::Real newNewtonThreshold)
{
    m_newtonThreshold = newNewtonThreshold;
}

int NBESolver::maxNewtonIterations() const
{
    return m_maxNewtonIterations;
}

void NBESolver::setMaxNewtonIterations(int newMaxNewtonIterations)
{
    m_maxNewtonIterations = newMaxNewtonIterations;
}

int NBESolver::lineSearchMaxIterations() const
{
    return m_lineSearchMaxIterations;
}

void NBESolver::setLineSearchMaxIterations(int newLineSearchMaxIterations)
{
    m_lineSearchMaxIterations = newLineSearchMaxIterations;
}
}  // namespace VNCS
