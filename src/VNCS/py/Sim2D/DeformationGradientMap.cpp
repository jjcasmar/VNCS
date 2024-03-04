#include "DeformationGradientMap.h"

#include <VNCS/FMap.h>
#include <VNCS/SamplingPoints.h>
#include <VNCS/py/SamplingPoints.h>
#include <sofa/core/sptr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

template <class T, class Space>
void registerFMap(pybind11::module &m, const std::string &name)
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

void VNCS::Sim2D::py::deformationGradientMap(pybind11::module &m)
{
    {
        using FMap = VNCS::FMap<VNCS::Space2D::F22>;
        registerFMap<FMap, VNCS::Space2D::F22::MaterialSpace>(m, "FMap22");
    }

    {
        using FMap = VNCS::FMap<VNCS::Space2D::F21>;
        registerFMap<FMap, VNCS::Space2D::F21::MaterialSpace>(m, "FMap21");
    }
}
