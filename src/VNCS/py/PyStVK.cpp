#include "PyStVK.h"
#include <VNCS/LinearStVKForceField.h>
#include <VNCS/StVKForceField.h>
#include <sofa/core/sptr.h>

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::py::module::stVK(pybind11::module &m)
{
    pybind11::class_<VNCS::StVKForceField, sofa::core::sptr<VNCS::StVKForceField>>(m, "StVK")
        .def(pybind11::init([](pybind11::args &args, pybind11::kwargs &kwargs) { return new VNCS::StVKForceField(); }))
        .def_property("poisson", &StVKForceField::poisson, &StVKForceField::setPoisson)
        .def_property("young", &StVKForceField::young, &StVKForceField::setYoung);

    pybind11::class_<VNCS::LinearStVKForceField, sofa::core::sptr<VNCS::LinearStVKForceField>>(m, "LinearStVK")
        .def(pybind11::init(
            [](pybind11::args &args, pybind11::kwargs &kwargs) { return new VNCS::LinearStVKForceField(); }))
        .def_property("poisson", &LinearStVKForceField::poisson, &LinearStVKForceField::setPoisson)
        .def_property("young", &LinearStVKForceField::young, &LinearStVKForceField::setYoung);
}
