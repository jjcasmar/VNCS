#include "StVK.h"
#include <VNCS/Sim2D/StVKForceField.h>
#include <sofa/core/sptr.h>
#include <VNCS/py/SamplingPoints.h>

#include <VNCS/Spaces.h>

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim2D::py::stVK(pybind11::module &m)
{
    pybind11::class_<VNCS::Sim2D::StVKForceField,
                     sofa::core::objectmodel::BaseObject,
                     sofa::core::sptr<VNCS::Sim2D::StVKForceField>>(m, "StVK")
        .def(pybind11::init<>())
        .def_property("poisson", &VNCS::Sim2D::StVKForceField::poisson, &VNCS::Sim2D::StVKForceField::setPoisson)
        .def_property("young", &VNCS::Sim2D::StVKForceField::young, &VNCS::Sim2D::StVKForceField::setYoung)
        .def_property("beta", &VNCS::Sim2D::StVKForceField::beta, &VNCS::Sim2D::StVKForceField::setBeta)
        .def("setSamplingPoints",
             [](VNCS::Sim2D::StVKForceField &ff,
                const VNCS::py::SamplingPointsHolder<VNCS::Space2D> &samplingPointsHolder) {
                 ff.setSamplingPoints(samplingPointsHolder.samplingPoints);
             });
}
