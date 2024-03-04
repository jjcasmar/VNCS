#include "Projection.h"

#include <VNCS/Projection.h>
#include <sofa/core/sptr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <memory>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim2D::py::projection(pybind11::module &m)
{
    using Inherit = sofa::core::objectmodel::BaseObject;
    ::py::class_<VNCS::Sim2D::Projection, Inherit, sofa::core::sptr<VNCS::Sim2D::Projection>>(m, "Projection")
        .def(pybind11::init<>())
        .def("setClusterMatrix",
             [](VNCS::Sim2D::Projection &projection, const Eigen::SparseMatrix<VNCS::Sim2D::Projection::Real> &m) {
                 projection.setClusterMatrix(m);
             })
        .def("setBarycentricMatrix",
             [](VNCS::Sim3D::Projection &projection, const Eigen::SparseMatrix<VNCS::Sim3D::Projection::Real> &m) {
                 projection.setBarycentricMatrix(m);
             });
}
