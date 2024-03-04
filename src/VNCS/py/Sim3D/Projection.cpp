#include "Projection.h"

#include <VNCS/Projection.h>
#include <sofa/core/sptr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/chrono.h>
#include <memory>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim3D::py::projection(pybind11::module &m)
{
    using Inherit = sofa::core::objectmodel::BaseObject;
    ::py::class_<VNCS::Sim3D::Projection, Inherit, sofa::core::sptr<VNCS::Sim3D::Projection>>(m, "Projection")
        .def(pybind11::init<>())
        .def("setClusterMatrix",
             [](VNCS::Sim3D::Projection &projection, const Eigen::SparseMatrix<VNCS::Sim3D::Projection::Real> &m) {
                 projection.setClusterMatrix(m);
             })
        .def("setBarycentricMatrix",
             [](VNCS::Sim3D::Projection &projection, const Eigen::SparseMatrix<VNCS::Sim3D::Projection::Real> &m) {
                 projection.setBarycentricMatrix(m);
             })
        .def_property_readonly("durations", &VNCS::Sim3D::Projection::durations)
        .def("clearDurations", &VNCS::Sim3D::Projection::clearDurations);
}
