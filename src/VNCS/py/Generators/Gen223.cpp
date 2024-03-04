#include "Gen223.h"

#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include <VNCS/Generators/Gen223/Gen223.h>
#include <VNCS/Generators/TetraMeshGenerator.h>
#include <Eigen/Dense>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>

#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>

namespace py = pybind11;

void VNCS::Generators::py::gen223(::py::module &m)
{
    ::py::class_<VNCS::Generators::Gen223::Generator>(m, "Gen223")
        .def(pybind11::init<>())
        .def_property(
            "coarseInputMeshPath",  //
            &VNCS::Generators::Gen223::Generator::coarseMeshPath,
            [](VNCS::Generators::Gen223::Generator &Gen323, std::string path) { Gen323.setCoarseMeshPath(path); })
        .def_property("coarseCriteria",  //
                      &VNCS::Generators::Gen223::Generator::coarseCriteria,
                      &VNCS::Generators::Gen223::Generator::setCoarseCriteria)
        .def_property(
            "fineInputMeshPath",  //
            &VNCS::Generators::Gen223::Generator::fineMeshPath,
            [](VNCS::Generators::Gen223::Generator &Gen323, std::string path) { Gen323.setFineMeshPath(path); })
        .def_property("fineCriteria",  //
                      &VNCS::Generators::Gen223::Generator::fineCriteria,
                      &VNCS::Generators::Gen223::Generator::setFineCriteria)
        .def_property(
            "coarseSamplersFilePath",
            [](const VNCS::Generators::Gen223::Generator &mo) {},
            [](VNCS::Generators::Gen223::Generator &Gen323, std::string &coarseSamplersFilePath) {
                Gen323.setCoarseSamplersFilePath(coarseSamplersFilePath);
            })
        .def_property(
            "fineSamplersFilePath",
            [](const VNCS::Generators::Gen223::Generator &mo) {},
            [](VNCS::Generators::Gen223::Generator &Gen323, std::string &fineSamplersFilePath) {
                Gen323.setFineSamplersFilePath(fineSamplersFilePath);
            })
        .def_property(
            "clusterMatrixFilePath",
            [](const VNCS::Generators::Gen223::Generator &mo) {},
            [](VNCS::Generators::Gen223::Generator &Gen323, std::string &clusterMatrixFilePath) {
                Gen323.setClusterMatrixFilePath(clusterMatrixFilePath);
            })
        .def_property(
            "dofFilePath",
            [](const VNCS::Generators::Gen223::Generator &mo) {},
            [](VNCS::Generators::Gen223::Generator &Gen323, std::string &dofFilePath) {
                Gen323.setDofFilePath(dofFilePath);
            })
        .def("__call__", &VNCS::Generators::Gen223::Generator::operator());
}
