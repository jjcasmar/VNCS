#include "MechanicalObject.h"
#include <sofa/core/sptr.h>

#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>

#include <VNCS/MechanicalObject.h>
#include <Eigen/Dense>

#include <VNCS/Preconditioner.h>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim3D::py::mechanicalObject(::py::module &m)
{
    ::py::class_<VNCS::Sim3D::MO, sofa::core::objectmodel::BaseObject, sofa::core::sptr<VNCS::Sim3D::MO>>(m, "MO")
        .def(::py::init([]() { return sofa::core::objectmodel::New<VNCS::Sim3D::MO>(); }))
        .def_property_readonly("restPosition",
                               [](const VNCS::Sim3D::MO &mo) {
                                   auto positions = mo.readRestPositions();
                                   return Eigen::Map<const Eigen::VectorXd>(&positions[0][0], 3 * positions.size());
                               })
        .def_property_readonly("displacement", [](const VNCS::Sim3D::MO &mo) {
            auto positions = mo.readPositions();
            return Eigen::Map<const Eigen::VectorXd>(&positions[0][0], 3 * positions.size());
        });
}
