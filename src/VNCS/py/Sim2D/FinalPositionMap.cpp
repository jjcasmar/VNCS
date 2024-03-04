#include "FinalPositionMap.h"

#include <sofa/core/BaseState.h>
#include <sofa/core/objectmodel/BaseObject.h>

#include <VNCS/Sim2D/FinalPositionMap.h>
#include <sofa/core/sptr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim2D::py::finalPositionMap(pybind11::module &m)
{
    using Inherit = sofa::core::objectmodel::BaseObject;

    ::py::class_<VNCS::Sim2D::FinalPositionMap, Inherit, sofa::core::sptr<VNCS::Sim2D::FinalPositionMap>>(
        m, "FinalPositionMap")
        .def(::py::init<>())
        .def("setInput",
             [](VNCS::Sim2D::FinalPositionMap &dMap, sofa::core::objectmodel::BaseObject *input) {
                 dMap.setFrom(dynamic_cast<sofa::core::BaseState *>(input));
             })
        .def("setOutput", [](VNCS::Sim2D::FinalPositionMap &dMap, sofa::core::objectmodel::BaseObject *output) {
            dMap.setTo(dynamic_cast<sofa::core::BaseState *>(output));
        });
}
