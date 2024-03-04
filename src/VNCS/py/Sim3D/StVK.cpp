#include "StVK.h"
#include <VNCS/Sim3D/StableNeoHookean.h>
#include <VNCS/Sim3D/StVKForceField.h>
#include <VNCS/Sim3D/TriangleFabricBendingEnergyCommon.h>
#include <VNCS/Sim3D/StretchForceField.h>
#include <VNCS/Sim3D/BendingForceField.h>
#include <VNCS/py/SamplingPoints.h>
#include <sofa/core/sptr.h>

#include <VNCS/Spaces.h>

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim3D::py::stVK(pybind11::module &m)
{
    pybind11::class_<VNCS::Sim3D::StableNeoHookean,
                     sofa::core::objectmodel::BaseObject,
                     sofa::core::sptr<VNCS::Sim3D::StableNeoHookean>>(m, "StableNeoHookean")
        .def(pybind11::init<>())
        .def("set_lame_mu", &VNCS::Sim3D::StableNeoHookean::setMu)
        .def("set_lame_lambda", &VNCS::Sim3D::StableNeoHookean::setLambda)
        .def("setSamplingPoints",
             [](VNCS::Sim3D::StableNeoHookean &ff,
                const VNCS::py::SamplingPointsHolder<VNCS::Space3D> &samplingPointsHolder) {
                 ff.setSamplingPoints(samplingPointsHolder.samplingPoints);
             });

    pybind11::class_<VNCS::Sim3D::StVKForceField,
                     sofa::core::objectmodel::BaseObject,
                     sofa::core::sptr<VNCS::Sim3D::StVKForceField>>(m, "StVK")
        .def(pybind11::init<>())
        .def_property("lame_mu", &VNCS::Sim3D::StVKForceField::mu, &VNCS::Sim3D::StVKForceField::setMu)
        .def_property("lame_lambda", &VNCS::Sim3D::StVKForceField::lambda, &VNCS::Sim3D::StVKForceField::setLambda)
        .def_property("beta", &VNCS::Sim3D::StVKForceField::beta, &VNCS::Sim3D::StVKForceField::setBeta)
        .def("setSamplingPoints",
             [](VNCS::Sim3D::StVKForceField &ff,
                const VNCS::py::SamplingPointsHolder<VNCS::Space2D> &samplingPointsHolder) {
                 ff.setSamplingPoints(samplingPointsHolder.samplingPoints);
             });

    pybind11::class_<VNCS::Sim3D::TriangleFabricBendingEnergy,
                     sofa::core::objectmodel::BaseObject,
                     sofa::core::sptr<VNCS::Sim3D::TriangleFabricBendingEnergy>>(m, "DiscreteShells")
        .def(pybind11::init<>())
        .def_property("stiffness",
                      &VNCS::Sim3D::TriangleFabricBendingEnergy::bendingStiffness,
                      &VNCS::Sim3D::TriangleFabricBendingEnergy::setBendingStiffness)
        .def_property("beta", &VNCS::Sim3D::TriangleFabricBendingEnergy::beta, &VNCS::Sim3D::TriangleFabricBendingEnergy::setBeta)
        .def_property(
            "meshPath",
            [](const VNCS::Sim3D::TriangleFabricBendingEnergy &e) -> std::string { return e.meshFilepath(); },
            [](VNCS::Sim3D::TriangleFabricBendingEnergy &e, const std::string &s) { return e.setMeshFilepath(s); });
}
