#include "Gen202.h"

#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <VNCS/Generators/Gen202/Gen202.h>
#include "BlendingField.h"
#include "AdaptiveMeshCriteriaFunctor.h"
#include <Eigen/Dense>

namespace py = pybind11;

void VNCS::Generators::py::gen202(::py::module &m)
{
    ::py::class_<VNCS::Generators::Gen202::Generator>(m, "Gen202")
        .def(pybind11::init<>())
        .def_property(
            "coarseInputMeshPath",
            [](const VNCS::Generators::Gen202::Generator &mo) {},
            [](VNCS::Generators::Gen202::Generator &Gen202, const std::string &meshPath) {
                Gen202.setCoarseMeshPath(meshPath);
            })
        .def_property("remeshCoarseInputMesh",
                      &VNCS::Generators::Gen202::Generator::remeshCoarseMesh,
                      &VNCS::Generators::Gen202::Generator::setRemeshCoarseMesh)
        .def_property(
            "coarseCriteria",
            [](const VNCS::Generators::Gen202::Generator &mo) { return mo.coarseCriteria(); },
            [](VNCS::Generators::Gen202::Generator &Gen202,
               std::shared_ptr<VNCS::Generators::py::AdaptiveMeshCriteriaFunctor> criteria) {
                Gen202.setCoarseCriteria(criteria);
            })
        .def_property(
            "coarseSamplersFilePath",
            [](const VNCS::Generators::Gen202::Generator &mo) {},
            [](VNCS::Generators::Gen202::Generator &Gen202, std::string &coarseSamplersFilePath) {
                Gen202.setCoarseSamplersFilePath(coarseSamplersFilePath);
            })
        .def_property(
            "dofFilePath",
            [](const VNCS::Generators::Gen202::Generator &mo) {},
            [](VNCS::Generators::Gen202::Generator &Gen202, std::string &dofFilePath) {
                Gen202.setDofFilePath(dofFilePath);
            })
        .def("create", &VNCS::Generators::Gen202::Generator::operator());
}
