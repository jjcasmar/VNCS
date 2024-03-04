#include "DeformationGradientMO.h"
#include <sofa/core/sptr.h>

#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <VNCS/MechanicalObject.h>
#include <Eigen/Dense>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim2D::py::deformationGradientMO(::py::module &m)
{
    ::py::class_<VNCS::Sim2D::FMO22, sofa::core::objectmodel::BaseObject, sofa::core::sptr<VNCS::Sim2D::FMO22>>(m,
                                                                                                                "FMO22")
        .def(::py::init([]() { return sofa::core::objectmodel::New<VNCS::Sim2D::FMO22>(); }))
        .def_property_readonly("F",
                               [](const VNCS::Sim2D::FMO22 &fMo) {
                                   std::vector<VNCS::Sim2D::FMO22::Coord> gradients;
                                   std::copy(std::begin(fMo.readPositions()),
                                             std::end(fMo.readPositions()),
                                             std::back_inserter(gradients));
                                   return gradients;
                               })
        .def_property_readonly("S", [](const VNCS::Sim2D::FMO22 &fMo) {
            std::vector<VNCS::Sim2D::FMO22::Coord> gradients;
            std::copy(std::begin(fMo.readForces()), std::end(fMo.readForces()), std::back_inserter(gradients));
            return gradients;
        });

    ::py::class_<VNCS::Sim2D::FMO21, sofa::core::objectmodel::BaseObject, sofa::core::sptr<VNCS::Sim2D::FMO21>>(m,
                                                                                                                "FMO21")
        .def(::py::init([]() { return sofa::core::objectmodel::New<VNCS::Sim2D::FMO21>(); }))
        .def_property_readonly("F",
                               [](const VNCS::Sim2D::FMO21 &fMo) {
                                   std::vector<VNCS::Sim2D::FMO21::Coord> gradients;
                                   std::copy(std::begin(fMo.readPositions()),
                                             std::end(fMo.readPositions()),
                                             std::back_inserter(gradients));
                                   return gradients;
                               })
        .def_property_readonly("S", [](const VNCS::Sim2D::FMO21 &fMo) {
            std::vector<VNCS::Sim2D::FMO21::Coord> gradients;
            std::copy(std::begin(fMo.readForces()), std::end(fMo.readForces()), std::back_inserter(gradients));
            return gradients;
        });
}
