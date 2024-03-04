#include "FinalPositionMap.h"

#include <sofa/core/BaseState.h>
#include <sofa/core/objectmodel/BaseObject.h>

#include <VNCS/Sim3D/FinalPositionMap.h>
#include <sofa/core/sptr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim3D::py::finalPositionMap(pybind11::module &m)
{
    using Inherit = sofa::core::objectmodel::BaseObject;

    ::py::class_<VNCS::Sim3D::FinalPositionMap, Inherit, sofa::core::sptr<VNCS::Sim3D::FinalPositionMap>>(
        m, "FinalPositionMap")
        .def(::py::init<>())
        .def("setInput",
             [](VNCS::Sim3D::FinalPositionMap &dMap, sofa::core::objectmodel::BaseObject *input) {
                 dMap.setFrom(dynamic_cast<sofa::core::BaseState *>(input));
             })
        .def("setOutput", [](VNCS::Sim3D::FinalPositionMap &dMap, sofa::core::objectmodel::BaseObject *output) {
            dMap.setTo(dynamic_cast<sofa::core::BaseState *>(output));
        });
}
