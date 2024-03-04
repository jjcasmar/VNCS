#include "UfMap.h"

#include <sofa/core/BaseState.h>
#include <sofa/core/objectmodel/BaseObject.h>

#include <VNCS/UfMap.h>
#include <sofa/core/sptr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim3D::py::ufMap(pybind11::module &m)
{
    using Inherit = sofa::core::objectmodel::BaseObject;
    ::py::class_<VNCS::Sim3D::UfMap, Inherit, sofa::core::sptr<VNCS::Sim3D::UfMap>>(m, "UfMap")
        .def(::py::init<>())
        .def("setClusterMatrix",
             [](VNCS::Sim3D::UfMap &ufMap, const Eigen::SparseMatrix<VNCS::Sim3D::UfMap::Real> &m) {
                 ufMap.setClusterMatrix(m);
             })
        .def("add_input",
             [](VNCS::Sim3D::UfMap &ufMap, ::py::handle input) {
                 sofa::core::BaseState *inputBS =
                     dynamic_cast<sofa::core::BaseState *>(input.cast<sofa::core::objectmodel::BaseObject *>());
                 ufMap.addInputModel(inputBS);
             })
        .def("setOutput", [](VNCS::Sim3D::UfMap &ufMap, ::py::handle output) {
            sofa::core::BaseState *outputBS =
                dynamic_cast<sofa::core::BaseState *>(output.cast<sofa::core::objectmodel::BaseObject *>());
            ufMap.addOutputModel(outputBS);
        });
}
