#include "Mass.h"

#include <VNCS/Mass.h>
#include <VNCS/SamplingPoints.h>
#include <VNCS/py/SamplingPoints.h>
#include <sofa/core/sptr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim2D::py::mass(pybind11::module &m)
{
    using Inherit = sofa::core::objectmodel::BaseObject;
    ::py::class_<VNCS::Sim2D::Mass22, Inherit, sofa::core::sptr<VNCS::Sim2D::Mass22>>(m, "Mass22")
        .def(::py::init<>())
        .def("setDensity", [](VNCS::Sim2D::Mass22 &mass, Mass22::Real density) { mass.setDensity(density); })
        .def("setSamplingPoints",
             [](VNCS::Sim2D::Mass22 &mass, const VNCS::py::SamplingPointsHolder<VNCS::Space2D> &samplingPoints) {
                 mass.setSamplingPoints(samplingPoints.samplingPoints);
             });

    ::py::class_<VNCS::Sim2D::Mass21, Inherit, sofa::core::sptr<VNCS::Sim2D::Mass21>>(m, "Mass21")
        .def(::py::init<>())
        .def("setDensity", [](VNCS::Sim2D::Mass21 &mass, Mass21::Real density) { mass.setDensity(density); })
        .def("setSamplingPoints",
             [](VNCS::Sim2D::Mass21 &mass, const VNCS::py::SamplingPointsHolder<VNCS::Space1D> &samplingPoints) {
                 mass.setSamplingPoints(samplingPoints.samplingPoints);
             });
}
