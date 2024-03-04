#include "ConcatMap.h"

#include <VNCS/ConcatMap.h>
#include <sofa/core/sptr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim2D::py::concatMap(pybind11::module &m)
{
    using Inherit = sofa::core::objectmodel::BaseObject;
    ::py::class_<VNCS::Sim2D::ConcatMap, Inherit, sofa::core::sptr<VNCS::Sim2D::ConcatMap>>(m, "ConcatMap")
        .def(::py::init<>())
        .def("add_input",
             [](VNCS::Sim2D::ConcatMap &concatMap, sofa::core::objectmodel::BaseObject *input) {
                 concatMap.addInputModel(dynamic_cast<sofa::core::BaseState *>(input));
             })
        .def("setOutput", [](VNCS::Sim2D::ConcatMap &concatMap, sofa::core::objectmodel::BaseObject *output) {
            concatMap.addOutputModel(dynamic_cast<sofa::core::BaseState *>(output));
        });
}
