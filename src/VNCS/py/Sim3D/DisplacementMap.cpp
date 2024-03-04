#include "DisplacementMap.h"

#include <VNCS/UMap.h>
#include <VNCS/SamplingPoints.h>
#include <VNCS/py/SamplingPoints.h>
#include <sofa/core/sptr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

template <class T, class Space>
void registerUMap(pybind11::module &m, const std::string &name)
{
    using Inherit = sofa::core::objectmodel::BaseObject;
    ::py::class_<T, Inherit, sofa::core::sptr<T>>(m, name.c_str())
        .def(::py::init<>())
        .def_property("outputMO", &T::output, &T::setOutput)
        .def_property("coarseMO", &T::coarseInput, &T::setCoarseInput)
        .def_property("fineMO", &T::fineInput, &T::setFineInput)
        .def("setSamplingPoints",
             [](T &fMap, const VNCS::py::SamplingPointsHolder<Space> &samplingPointsHolder) {
                 fMap.setSamplingPoints(samplingPointsHolder.samplingPoints);
             })
        .def_property_readonly("phiC", [](const T &fMap) { return fMap.phiC(); })
        .def_property_readonly("phiF", [](const T &fMap) { return fMap.phiF(); });
}

void VNCS::Sim3D::py::displacementMap(pybind11::module &m)
{
    registerUMap<VNCS::Sim3D::UMap33, VNCS::Space3D>(m, "UMap33");
    registerUMap<VNCS::Sim3D::UMap32, VNCS::Space2D>(m, "UMap32");
    registerUMap<VNCS::Sim3D::UMap31, VNCS::Space1D>(m, "UMap31");
}
