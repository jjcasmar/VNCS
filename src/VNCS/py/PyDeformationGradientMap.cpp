#include "PyDeformationGradientMap.h"

#include <VNCS/DeformationGradientMap.h>

#include <pybind11/pybind11.h>
#include <sofa/core/sptr.h>
#include <pybind11/eigen.h>

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::py::module::deformationGradientMap(pybind11::module &m)
{
    pybind11::class_<VNCS::DeformationGradientMap, sofa::core::sptr<VNCS::DeformationGradientMap>>(
        m, "DeformationGradientMap")
        .def_property_readonly("phiF", &VNCS::DeformationGradientMap::phiF)
        .def_property_readonly("phiC", &VNCS::DeformationGradientMap::phiC);
}
