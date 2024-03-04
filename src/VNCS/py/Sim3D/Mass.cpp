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

void VNCS::Sim3D::py::mass(pybind11::module &m)
{
    using Inherit = sofa::core::objectmodel::BaseObject;
    ::py::class_<VNCS::Sim3D::Mass33, Inherit, sofa::core::sptr<VNCS::Sim3D::Mass33>>(m, "Mass33")
        .def(::py::init<>())
        .def("setDensity", [](VNCS::Sim3D::Mass33 &mass, Mass33::Real density) { mass.setDensity(density); })
        .def("setSamplingPoints",
             [](VNCS::Sim3D::Mass33 &mass, const VNCS::py::SamplingPointsHolder<VNCS::Space3D> &samplingPoints) {
                 mass.setSamplingPoints(samplingPoints.samplingPoints);
             });

    ::py::class_<VNCS::Sim3D::Mass32, Inherit, sofa::core::sptr<VNCS::Sim3D::Mass32>>(m, "Mass32")
        .def(::py::init<>())
        .def("setDensity", [](VNCS::Sim3D::Mass32 &mass, Mass32::Real density) { mass.setDensity(density); })
        .def("setSamplingPoints",
             [](VNCS::Sim3D::Mass32 &mass, const VNCS::py::SamplingPointsHolder<VNCS::Space2D> &samplingPoints) {
                 mass.setSamplingPoints(samplingPoints.samplingPoints);
             });

    ::py::class_<VNCS::Sim3D::Mass31, Inherit, sofa::core::sptr<VNCS::Sim3D::Mass31>>(m, "Mass31")
        .def(::py::init<>())
        .def("setDensity", [](VNCS::Sim3D::Mass31 &mass, Mass31::Real density) { mass.setDensity(density); })
        .def("setSamplingPoints",
             [](VNCS::Sim3D::Mass31 &mass, const VNCS::py::SamplingPointsHolder<VNCS::Space1D> &samplingPoints) {
                 mass.setSamplingPoints(samplingPoints.samplingPoints);
             });
}
