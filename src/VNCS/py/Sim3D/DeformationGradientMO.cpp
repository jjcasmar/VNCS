#include "DeformationGradientMO.h"
#include <sofa/core/sptr.h>

#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <VNCS/MechanicalObject.h>
#include <Eigen/Dense>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim3D::py::deformationGradientMO(::py::module &m)
{
    ::py::class_<VNCS::Sim3D::FMO33, sofa::core::objectmodel::BaseObject, sofa::core::sptr<VNCS::Sim3D::FMO33>>(m,
                                                                                                                "FMO33")
        .def(::py::init([]() { return sofa::core::objectmodel::New<VNCS::Sim3D::FMO33>(); }))
        .def_property_readonly("F",
                               [](const VNCS::Sim3D::FMO33 &fMo) {
                                   std::vector<VNCS::Sim3D::FMO33::Coord> gradients;
                                   std::copy(std::begin(fMo.readPositions()),
                                             std::end(fMo.readPositions()),
                                             std::back_inserter(gradients));
                                   return gradients;
                               })
        .def_property_readonly("S", [](const VNCS::Sim3D::FMO33 &fMo) {
            std::vector<VNCS::Sim3D::FMO33::Coord> gradients;
            std::copy(std::begin(fMo.readForces()), std::end(fMo.readForces()), std::back_inserter(gradients));
            return gradients;
        });

    ::py::class_<VNCS::Sim3D::FMO32, sofa::core::objectmodel::BaseObject, sofa::core::sptr<VNCS::Sim3D::FMO32>>(m,
                                                                                                                "FMO32")
        .def(::py::init([]() { return sofa::core::objectmodel::New<VNCS::Sim3D::FMO32>(); }))
        .def_property_readonly("F",
                               [](const VNCS::Sim3D::FMO32 &fMo) {
                                   std::vector<VNCS::Sim3D::FMO32::Coord> gradients;
                                   std::copy(std::begin(fMo.readPositions()),
                                             std::end(fMo.readPositions()),
                                             std::back_inserter(gradients));
                                   return gradients;
                               })
        .def_property_readonly("S", [](const VNCS::Sim3D::FMO32 &fMo) {
            std::vector<VNCS::Sim3D::FMO32::Coord> gradients;
            std::copy(std::begin(fMo.readForces()), std::end(fMo.readForces()), std::back_inserter(gradients));
            return gradients;
        });

    ::py::class_<VNCS::Sim3D::FMO31, sofa::core::objectmodel::BaseObject, sofa::core::sptr<VNCS::Sim3D::FMO31>>(m,
                                                                                                                "FMO31")
        .def(::py::init([]() { return sofa::core::objectmodel::New<VNCS::Sim3D::FMO31>(); }))
        .def_property_readonly("F",
                               [](const VNCS::Sim3D::FMO31 &fMo) {
                                   std::vector<VNCS::Sim3D::FMO31::Coord> gradients;
                                   std::copy(std::begin(fMo.readPositions()),
                                             std::end(fMo.readPositions()),
                                             std::back_inserter(gradients));
                                   return gradients;
                               })
        .def_property_readonly("S", [](const VNCS::Sim3D::FMO31 &fMo) {
            std::vector<VNCS::Sim3D::FMO31::Coord> gradients;
            std::copy(std::begin(fMo.readForces()), std::end(fMo.readForces()), std::back_inserter(gradients));
            return gradients;
        });
}
