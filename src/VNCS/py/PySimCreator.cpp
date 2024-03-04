#include "PySimCreator.h"
#include <VNCS/SimCreator.h>
#include <VNCS/DeformationGradientMap.h>

#include <pybind11/pybind11.h>
#include <sofa/core/sptr.h>
#include <pybind11/eigen.h>

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::py::module::simCreator(pybind11::module &m)
{
    pybind11::class_<VNCS::SimCreator, sofa::core::sptr<VNCS::SimCreator>>(m, "SimCreator")
        .def(pybind11::init([](pybind11::args &args, pybind11::kwargs &kwargs) { return new VNCS::SimCreator(); }))
        .def_property("blendingField", &VNCS::SimCreator::blendingField, &VNCS::SimCreator::setBlendingField)
        .def_property_readonly("C", &SimCreator::C);
}
