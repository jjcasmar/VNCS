#include "AdaptiveMeshCriteriaFunctor.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <memory>

std::pair<VNCS::Generators::py::AdaptiveMeshCriteriaFunctor::Real,
          VNCS::Generators::py::AdaptiveMeshCriteriaFunctor::Real>
VNCS::Generators::py::AdaptiveMeshCriteriaFunctor::operator()(const Space2D::Point &p) const
{
    using Inherit = VNCS::Generators::AdaptiveMeshCriteriaFunctor;
    using RType = std::pair<VNCS::Space2D::Real, VNCS::Space2D::Real>;
    const std::array<double, 2> pyP{p[0], p[1]};
    PYBIND11_OVERLOAD_PURE_NAME(RType,      /* Return type */
                                Inherit,    /* Parent class */
                                "__call__", /* Name of function in C++ (must match Python name) */
                                operator(),
                                pyP /* Argument(s) */
    );
}

void VNCS::Generators::py::adaptiveMeshCriteriaFunctor(pybind11::module &m)
{
    using Inherit = VNCS::Generators::AdaptiveMeshCriteriaFunctor;
    pybind11::class_<Inherit, VNCS::Generators::py::AdaptiveMeshCriteriaFunctor, std::shared_ptr<Inherit>>(
        m, "AdaptiveMeshCriteria")
        .def(pybind11::init<>())
        .def("__call__", &Inherit::operator());
}
