#include "ConjugateGradient.h"

#include <VNCS/ConjugateGradientSolver.h>
#include <VNCS/NBESolver.h>
#include <pybind11/stl.h>

#include <pybind11/pybind11.h>

#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <SofaEigen2Solver/EigenVectorWrapper.h>
#include <SofaBaseLinearSolver/CGLinearSolver.h>

#include <memory>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::py::conjugateGradient(pybind11::module &m)
{
    {
        using SofaCG =
            sofa::component::linearsolver::CGLinearSolver<sofa::component::linearsolver::GraphScatteredMatrix,
                                                          sofa::component::linearsolver::GraphScatteredVector>;
        using VNCSCG = VNCS::ConjugateGradientSolver;
        m.def("residuals", [](const sofa::core::objectmodel::BaseObject *cgBaseObject) -> std::vector<VNCS::Real> {
            {
                auto cg = dynamic_cast<const SofaCG *>(cgBaseObject);
                if (cg) {
                    auto error = cg->f_graph.getValue().at("Error");
                    return std::vector<VNCS::Real>(error.begin(), error.end());
                }
            }

            {
                auto cg = dynamic_cast<const VNCSCG *>(cgBaseObject);
                if (cg)
                    return cg->squaredResiduals();
            }
            return {};
        });
    }

    {
        using Inherit = sofa::core::objectmodel::BaseObject;
        using CG = VNCS::ConjugateGradientSolver;

        ::py::class_<CG, Inherit, sofa::core::sptr<CG>>(m, "ConjugateGradient")
            .def(::py::init([]() { return sofa::core::objectmodel::New<CG>(); }))
            .def_property("iterations", &CG::iterations, &CG::setIterations)
            .def_property("tolerance", &CG::toleranceThreshold, &CG::setToleranceThreshold)
            .def_property("preconditioner", &CG::preconditioner, &CG::setPreconditioner)
            .def_property_readonly("residuals", &CG::squaredResiduals)
            .def_property("verbose", &CG::isVerbose, &CG::setIsVerbose)
            .def_property_readonly("durations", &CG::durations)
            .def("clearDurations", &CG::clearDurations);

        ::py::class_<VNCS::NBESolver, Inherit, sofa::core::sptr<VNCS::NBESolver>>(m, "NBESolver")
            .def(::py::init<>())
            .def_property("newtonIterations", &NBESolver::maxNewtonIterations, &NBESolver::setMaxNewtonIterations)
            .def_property("newtonThreshold", &NBESolver::newtonThreshold, &NBESolver::setNewtonThreshold)
            .def_property("rayleighMass", &NBESolver::rayleighMass, &NBESolver::setRayleighMass)
            .def_property("rayleighStiffness", &NBESolver::rayleighStiffness, &NBESolver::setRayleighStiffness)
            .def_property(
                "lineSearchIterations", &NBESolver::lineSearchMaxIterations, &NBESolver::setLineSearchMaxIterations);

        m.def("residuals", [](const CG *cg) {
            std::cout << "Using CG approach\n";
            return cg->squaredResiduals();
        });
    }
}
