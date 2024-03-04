#include "Gen323.h"

#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include <VNCS/Generators/Gen323/Gen323.h>
#include <VNCS/Generators/TetraMeshGenerator.h>
#include <Eigen/Dense>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>

#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>

namespace py = pybind11;

void VNCS::Generators::py::gen323(::py::module &m)
{
    ::py::class_<VNCS::Generators::Gen323::Generator>(m, "Gen323")
        .def(pybind11::init<>())
        .def_property(
            "coarseInputMeshPath",  //
            &VNCS::Generators::Gen323::Generator::coarseMeshPath,
            [](VNCS::Generators::Gen323::Generator &Gen323, std::string path) { Gen323.setCoarseMeshPath(path); })
        .def_property("coarseCriteria",  //
                      &VNCS::Generators::Gen323::Generator::coarseCriteria,
                      &VNCS::Generators::Gen323::Generator::setCoarseCriteria)
        .def_property(
            "fineInputMeshPath",  //
            &VNCS::Generators::Gen323::Generator::fineMeshPath,
            [](VNCS::Generators::Gen323::Generator &Gen323, std::string path) { Gen323.setFineMeshPath(path); })
        .def_property("fineCriteria",  //
                      &VNCS::Generators::Gen323::Generator::fineCriteria,
                      &VNCS::Generators::Gen323::Generator::setFineCriteria)
        .def_property(
            "coarseSamplersFilePath",
            [](const VNCS::Generators::Gen323::Generator &mo) {},
            [](VNCS::Generators::Gen323::Generator &Gen323, std::string &coarseSamplersFilePath) {
                Gen323.setCoarseSamplersFilePath(coarseSamplersFilePath);
            })
        .def_property(
            "fineSamplersFilePath",
            [](const VNCS::Generators::Gen323::Generator &mo) {},
            [](VNCS::Generators::Gen323::Generator &Gen323, std::string &fineSamplersFilePath) {
                Gen323.setFineSamplersFilePath(fineSamplersFilePath);
            })
        .def_property(
            "clusterMatrixFilePath",
            [](const VNCS::Generators::Gen323::Generator &mo) {},
            [](VNCS::Generators::Gen323::Generator &Gen323, std::string &clusterMatrixFilePath) {
                Gen323.setClusterMatrixFilePath(clusterMatrixFilePath);
            })
        .def_property(
            "dofFilePath",
            [](const VNCS::Generators::Gen323::Generator &mo) {},
            [](VNCS::Generators::Gen323::Generator &Gen323, std::string &dofFilePath) {
                Gen323.setDofFilePath(dofFilePath);
            })
        .def("exportFunction",
             [](VNCS::Generators::Gen323::Generator &Gen323, pybind11::function f) {
                 Gen323.setExportFunction(
                     [f](const VNCS::Space3D::TetraMesh &mesh,
                         const std::string &filename,
                         const std::unordered_map<std::string, std::vector<VNCS::Space3D::Real>> &attributes) {
                         std::vector<std::array<VNCS::Space3D::Real, 3>> points =
                             mesh.points | ranges::views::transform([](const auto &p) {
                                 return std::array<VNCS::Space3D::Real, 3>{p[0], p[1], p[2]};
                             }) |
                             ranges::to_vector;

                         f(points, mesh.tetras, attributes, filename);
                     });
             })
        .def("__call__", &VNCS::Generators::Gen323::Generator::operator());
}
