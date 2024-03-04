#include "ConjugateGradientSolver.h"

#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <SofaEigen2Solver/EigenVectorWrapper.h>

#include <VNCS/Logger.h>

#include <chrono>
#include <fstream>

namespace Log
{
struct CG {
    constexpr static const char *name = "CG";

    static int createLogger()
    {
        std::shared_ptr<spdlog::logger> logger = spdlog::stdout_color_mt(name);
        return 0;
    }

    inline static const int logger = createLogger();
};  // namespace CG
}  // namespace Log

namespace VNCS
{
ConjugateGradientSolver::ConjugateGradientSolver()
    : m_iterations(25)
    , m_toleranceThreshold(1e-5)
{
    Logger<::Log::CG>::logger()->set_level(spdlog::level::off);
}

void ConjugateGradientSolver::resetSystem()
{
}

void ConjugateGradientSolver::setSystemMBKMatrix(const sofa::core::MechanicalParams *mparams)
{
    m_mparams = *mparams;
}

void ConjugateGradientSolver::setSystemRHVector(sofa::core::MultiVecDerivId b_id)
{
    m_mparams.setDf(b_id);
    m_bId = b_id;
}

void ConjugateGradientSolver::setSystemLHVector(sofa::core::MultiVecDerivId x_id)
{
    m_mparams.setDx(x_id);
    m_xId = x_id;
}

void ConjugateGradientSolver::solve(sofa::core::behavior::MultiVecDeriv &b, sofa::core::behavior::MultiVecDeriv &x)
{
    m_squaredResiduals.clear();
    m_squaredResiduals.reserve(m_iterations);

    sofa::simulation::common::VectorOperations vop(&m_mparams, this->getContext());
    sofa::simulation::common::MechanicalOperations mop(&m_mparams, this->getContext());

    // Create temporary vectors needed for the method
    sofa::core::behavior::MultiVecDeriv r(&vop);
    sofa::core::behavior::MultiVecDeriv p(&vop);
    sofa::core::behavior::MultiVecDeriv s(&vop);
    sofa::core::behavior::MultiVecDeriv h(&vop);

    // Get the method parameters
    const auto &maximum_number_of_iterations = m_iterations;
    const auto &residual_tolerance_threshold = m_toleranceThreshold;

    m_squaredResiduals.clear();
    m_squaredResiduals.reserve(m_iterations);

    // Get the matrices coefficient m, b and k : A = (mM + bB + kK)
    const auto m_coef = m_mparams.mFactor();
    const auto b_coef = m_mparams.bFactor();
    const auto k_coef = m_mparams.kFactor();

    // Declare the method variables
    VNCS::Real b_delta;
    VNCS::Real delta;            // RHS and residual squared norms
    VNCS::Real rho0, rho1 = 0.;  // Stores r*r as it is used two times per iterations
    VNCS::Real alpha, beta;      // Alpha and Beta coefficients
    VNCS::Real threshold;        // Residual threshold
    int currentIteration = 0;    // Current iteration number
    bool converged = false;

    sofa::component::linearsolver::DefaultMultiMatrixAccessor accessor;
    mop.getMatrixDimension(nullptr, nullptr, &accessor);
    const auto n = static_cast<Eigen::Index>(accessor.getGlobalDimension());

    mop.projectResponse(b);

    if (m_preconditioner) {
        Eigen::VectorXd bEigen;
        bEigen.resize(n);
        sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> bEigenWrapped(bEigen);
        mop.multiVector2BaseVector(b.id(), &bEigenWrapped, &accessor);

        b_delta = bEigen.transpose() * m_preconditioner.value() * bEigen;
    } else {
        b_delta = b.dot(b);
    }

    threshold = residual_tolerance_threshold * residual_tolerance_threshold * b_delta;

    // INITIAL RESIDUAL
    // Do the A*x(0) with visitors since we did not construct the matrix A
    mop.propagateDxAndResetDf(x, s);          // Calls applyJ(x) on every mechanical mappings
    mop.addMBKdx(s, m_coef, b_coef, k_coef);  // s = (m M + b B + k K) x
    r.eq(b, s, -1.0);                         // r = b - A*x

    // Do the projection of the result in the constrained space
    mop.projectResponse(r);  // r = S (b - Ax)

    Eigen::VectorXd rEigen;
    sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> rEigenWrapped(rEigen);
    rEigen.resize(n);
    if (m_preconditioner) {
        rEigen.resize(n);

        mop.multiVector2BaseVector(r.id(), &rEigenWrapped, &accessor);

        rEigen = m_preconditioner.value() * rEigen;  // It doesnt alias because m_preconditioner is diagonal

        mop.baseVector2MultiVector(&rEigenWrapped, p.id(), &accessor);  // Save it directly in p, (p = P * r)
        mop.projectResponse(p);                                         // p = S * P * r

        mop.multiVector2BaseVector(p.id(), &rEigenWrapped, &accessor);
    } else {
        p.eq(r);
    }

    delta = r.dot(p);

    // ITERATIONS
    m_squaredResiduals.push_back(delta);
    while (delta > threshold && currentIteration < maximum_number_of_iterations) {
        const auto start = std::chrono::steady_clock::now();
        // 1. Computes q(k+1) = A*p(k)
        mop.propagateDxAndResetDf(p, s);  // calls applyJ(p) on mechanical mappings

        mop.addMBKdx(s, m_coef, b_coef, k_coef);  // s = (m M + b B + k K) p

        // We need to project the residual in the constrained space since the constraints haven't been added to the
        // matrix (the matrix is never constructed) and addMBKdx of the forcefields do not take constraints into
        // account.
        mop.projectResponse(s);  // s = S*A*p

        // 2. Computes x(k+1) and r(k+1)
        alpha = delta / p.dot(s);
        x.peq(p, alpha);   // x = x + alpha*p
        r.peq(s, -alpha);  // r = r - alpha*s

        if (m_preconditioner) {
            mop.multiVector2BaseVector(r.id(), &rEigenWrapped, &accessor);

            rEigen = m_preconditioner.value() * rEigen;

            mop.baseVector2MultiVector(&rEigenWrapped, h.id(), &accessor);  // Save it directly in h, (h = P * r)
        } else {
            h.eq(r);
        }

        // 3. Computes the new residual norm
        VNCS::Real deltaOld = delta;
        delta = r.dot(h);

        // 6. Compute the next search direction
        p.eq(h, p, delta / deltaOld);  // p = r + beta*p
        mop.projectResponse(p);

        ++currentIteration;

        m_squaredResiduals.push_back(delta);
        m_durations.push_back(
            std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count());
    }

    // Check positive definiteness
    if (isVerbose()) {
        mop.propagateDxAndResetDf(x, s);          // calls applyJ(p) on mechanical mappings
        mop.addMBKdx(s, m_coef, b_coef, k_coef);  // s = (m M + b B + k K) p

        const VNCS::Real positiveDefiniteness = x.dot(s);

        s.peq(b, -1.0);
        const VNCS::Real error = s.norm();

        Logger<Log::CG>::info("CG finished:");
        Logger<Log::CG>::info("Iterations: {} / {}", currentIteration, maximum_number_of_iterations);
        Logger<Log::CG>::info("Threshold: {} / {}", delta, threshold);
        Logger<Log::CG>::info("Error ||Ax - b||: {}", error);
        Logger<Log::CG>::info("PSD: {}", positiveDefiniteness);
    }
}

bool ConjugateGradientSolver::isVerbose() const
{
    return Logger<::Log::CG>::logger()->level() == spdlog::level::info;
}

void ConjugateGradientSolver::setIsVerbose(bool newIsVerbose)
{
    if (newIsVerbose)
        Logger<::Log::CG>::logger()->set_level(spdlog::level::info);
    else
        Logger<::Log::CG>::logger()->set_level(spdlog::level::off);
}

// Eigen::SparseMatrix<VNCS::Real> ConjugateGradientSolver::systemMatrix()
//{
//    sofa::simulation::common::VectorOperations vop(&m_mparams, this->getContext());
//    sofa::simulation::common::MechanicalOperations mop(&m_mparams, this->getContext());

//    const auto m_coef = m_mparams.mFactor();
//    const auto b_coef = m_mparams.bFactor();
//    const auto k_coef = m_mparams.kFactor();

//    sofa::core::behavior::MultiVecDeriv dx(&vop);
//    sofa::core::behavior::MultiVecDeriv df(&vop);
//    dx.realloc(&vop, true, true);

//    sofa::component::linearsolver::DefaultMultiMatrixAccessor accessor;

//    mop.getMatrixDimension(nullptr, nullptr, &accessor);
//    const auto n = static_cast<Eigen::Index>(accessor.getGlobalDimension());

//    Eigen::VectorXd dxEigen;
//    dxEigen.resize(n);
//    dxEigen.setZero();

//    std::vector<Eigen::Triplet<VNCS::Real>> triplets;

//    sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> dxEigenWrapped(dxEigen);
//    for (int i = 0; i < n; ++i) {
//        dxEigen.setZero();
//        dxEigen[i] = 1;

//        mop.baseVector2MultiVector(&dxEigenWrapped, dx.id(), &accessor);

//        mop.propagateDxAndResetDf(dx.id(), df.id());
//        mop.addMBKdx(df.id(), m_coef, b_coef, k_coef);

//        mop.multiVector2BaseVector(df.id(), &dxEigenWrapped, &accessor);

//        for (int j = 0; j < n; ++j) {
//            if (dxEigen[j] > 0.0)
//                triplets.emplace_back(i, j, dxEigen[j]);
//        }
//    }

//    Eigen::SparseMatrix<VNCS::Real> M;
//    M.resize(n, n);

//    M.setFromTriplets(std::begin(triplets), std::end(triplets));
//    return M;
//}

const std::optional<Eigen::SparseMatrix<Real>> &ConjugateGradientSolver::preconditioner() const
{
    return m_preconditioner;
}

void ConjugateGradientSolver::setPreconditioner(std::optional<Eigen::SparseMatrix<Real>> newPreconditioner)
{
    m_preconditioner = newPreconditioner;
}

const std::vector<Real> &ConjugateGradientSolver::squaredResiduals() const
{
    return m_squaredResiduals;
}

// Eigen::SparseMatrix<Real> ConjugateGradientSolver::massMatrix()
//{
//    sofa::core::MechanicalParams mparams;
//    sofa::simulation::common::VectorOperations vop(&m_mparams, this->getContext());
//    sofa::simulation::common::MechanicalOperations mop(&m_mparams, this->getContext());

//    sofa::core::behavior::MultiVecDeriv dx(&vop);
//    sofa::core::behavior::MultiVecDeriv df(&vop);
//    mparams.setDx(dx);
//    mparams.setDf(df);
//    mparams.setX(sofa::core::ConstVecCoordId::position());
//    mparams.setV(sofa::core::ConstVecDerivId::velocity());
//    dx.realloc(&vop, true, true);

//    sofa::component::linearsolver::DefaultMultiMatrixAccessor accessor;

//    mop.getMatrixDimension(nullptr, nullptr, &accessor);
//    const auto n = static_cast<Eigen::Index>(accessor.getGlobalDimension());

//    Eigen::VectorXd dxEigen;
//    dxEigen.resize(n);
//    dxEigen.setZero();

//    std::vector<Eigen::Triplet<VNCS::Real>> triplets;

//    sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> dxEigenWrapped(dxEigen);
//    for (int i = 0; i < n; ++i) {
//        dxEigen.setZero();
//        dxEigen[i] = 1;

//        mop.baseVector2MultiVector(&dxEigenWrapped, dx.id(), &accessor);

//        mop.propagateDxAndResetDf(dx.id(), df.id());
//        mop.addMBKdx(df.id(), 1.0, 0.0, 0.0);

//        mop.multiVector2BaseVector(df.id(), &dxEigenWrapped, &accessor);

//        for (int j = 0; j < n; ++j) {
//            if (dxEigen[j] > 0.0)
//                triplets.emplace_back(i, j, dxEigen[j]);
//        }
//    }

//    Eigen::SparseMatrix<VNCS::Real> M;
//    M.resize(n, n);

//    M.setFromTriplets(std::begin(triplets), std::end(triplets));
//    return M;
//}

// Eigen::SparseMatrix<Real> ConjugateGradientSolver::stiffnessMatrix()
//{
//    sofa::core::MechanicalParams mparams;
//    sofa::simulation::common::VectorOperations vop(&m_mparams, this->getContext());
//    sofa::simulation::common::MechanicalOperations mop(&m_mparams, this->getContext());

//    sofa::core::behavior::MultiVecDeriv dx(&vop);
//    sofa::core::behavior::MultiVecDeriv df(&vop);
//    mparams.setDx(dx);
//    mparams.setDf(df);
//    mparams.setX(sofa::core::ConstVecCoordId::position());
//    mparams.setV(sofa::core::ConstVecDerivId::velocity());
//    dx.realloc(&vop, true, true);

//    sofa::component::linearsolver::DefaultMultiMatrixAccessor accessor;

//    mop.getMatrixDimension(nullptr, nullptr, &accessor);
//    const auto n = static_cast<Eigen::Index>(accessor.getGlobalDimension());

//    Eigen::VectorXd dxEigen;
//    dxEigen.resize(n);
//    dxEigen.setZero();

//    std::vector<Eigen::Triplet<VNCS::Real>> triplets;

//    sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> dxEigenWrapped(dxEigen);
//    for (int i = 0; i < n; ++i) {
//        dxEigen.setZero();
//        dxEigen[i] = 1;

//        mop.baseVector2MultiVector(&dxEigenWrapped, dx.id(), &accessor);

//        mop.propagateDxAndResetDf(dx.id(), df.id());
//        mop.addMBKdx(df.id(), 0.0, 0.0, 1.0);

//        mop.multiVector2BaseVector(df.id(), &dxEigenWrapped, &accessor);

//        for (int j = 0; j < n; ++j) {
//            if (dxEigen[j] > 0.0)
//                triplets.emplace_back(i, j, dxEigen[j]);
//        }
//    }

//    Eigen::SparseMatrix<VNCS::Real> M;
//    M.resize(n, n);

//    M.setFromTriplets(std::begin(triplets), std::end(triplets));
//    return M;
//}

VNCS::Real ConjugateGradientSolver::toleranceThreshold() const
{
    return m_toleranceThreshold;
}

void ConjugateGradientSolver::setToleranceThreshold(VNCS::Real newToleranceThreshold)
{
    m_toleranceThreshold = newToleranceThreshold;
}

int ConjugateGradientSolver::iterations() const
{
    return m_iterations;
}

void ConjugateGradientSolver::setIterations(int newIterations)
{
    m_iterations = newIterations;
}

void ConjugateGradientSolver::solveSystem()
{
    // Gather the x and b vector identifiers
    sofa::simulation::common::VectorOperations vop(&m_mparams, this->getContext());
    sofa::core::behavior::MultiVecDeriv x(&vop, m_xId);
    sofa::core::behavior::MultiVecDeriv b(&vop, m_bId);

    solve(b, x);
}

void ConjugateGradientSolver::clearDurations()
{
    m_durations.clear();
}

const std::vector<int> &ConjugateGradientSolver::durations() const
{
    return m_durations;
}

}  // namespace VNCS
