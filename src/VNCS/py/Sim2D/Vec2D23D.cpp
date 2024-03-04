#include "Vec2D23D.h"

#include <sofa/core/BaseState.h>
#include <sofa/core/objectmodel/BaseObject.h>

#include <VNCS/Sim2D/Vec2D23D.h>
#include <sofa/core/sptr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim2D::py::vec2D23D(pybind11::module &m)
{
    using Inherit = sofa::core::objectmodel::BaseObject;

    ::py::class_<VNCS::Sim2D::Vec2D23D, Inherit, sofa::core::sptr<VNCS::Sim2D::Vec2D23D>>(m, "Vec2D23D")
        .def(::py::init<>())
        .def("setInput",
             [](VNCS::Sim2D::Vec2D23D &dMap, sofa::core::objectmodel::BaseObject *input) {
                 dMap.setFrom(dynamic_cast<sofa::core::BaseState *>(input));
             })
        .def("setOutput", [](VNCS::Sim2D::Vec2D23D &dMap, sofa::core::objectmodel::BaseObject *output) {
            dMap.setTo(dynamic_cast<sofa::core::BaseState *>(output));
        });
}
