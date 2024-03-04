#include "Gen333.h"

#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include <VNCS/Generators/Gen333/Gen333.h>
#include <VNCS/Generators/TetraMeshGenerator.h>
#include <Eigen/Dense>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>

#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>

#include <VNCS/Generators/Gen333/VisualSamplersExporter.h>

namespace py = pybind11;

void VNCS::Generators::py::gen333(::py::module &m)
{
    ::py::class_<VNCS::Generators::C3t3Criteria>(m, "C3t3Criteria")
        .def(pybind11::init<>())
        .def_readwrite("facetAngle", &VNCS::Generators::C3t3Criteria::facet_angle)
        .def_readwrite("facetSize", &VNCS::Generators::C3t3Criteria::facet_size)
        .def_readwrite("facetDistance", &VNCS::Generators::C3t3Criteria::facet_distance)
        .def_readwrite("cellSize", &VNCS::Generators::C3t3Criteria::cell_size)
        .def_readwrite("cellRatio", &VNCS::Generators::C3t3Criteria::cell_radius_edge_ratio)
        .def_readwrite("sharpAngle", &VNCS::Generators::C3t3Criteria::sharpAngle)
        .def_readwrite("sharpEdgeSize", &VNCS::Generators::C3t3Criteria::sharpEdgeSize);

    ::py::class_<VNCS::Generators::RemeshCriteria>(m, "RemeshCriteria")
        .def(pybind11::init<>())
        .def_readwrite("target_edge_length", &VNCS::Generators::RemeshCriteria::target_edge_length)
        .def_readwrite("sharpAngle", &VNCS::Generators::RemeshCriteria::sharpAngle)
        .def_readwrite("iterations", &VNCS::Generators::RemeshCriteria::iterations);

    ::py::class_<VNCS::Generators::Gen333::Generator>(m, "Gen333")
        .def(pybind11::init<>())
        .def_property(
            "coarseInputMeshPath",  //
            &VNCS::Generators::Gen333::Generator::coarseMeshPath,
            [](VNCS::Generators::Gen333::Generator &gen333, std::string path) { gen333.setCoarseMeshPath(path); })
        .def_property("coarseCriteria",  //
                      &VNCS::Generators::Gen333::Generator::coarseCriteria,
                      &VNCS::Generators::Gen333::Generator::setCoarseCriteria)
        .def_property(
            "fineInputMeshPath",  //
            &VNCS::Generators::Gen333::Generator::fineMeshPath,
            [](VNCS::Generators::Gen333::Generator &gen333, std::string path) { gen333.setFineMeshPath(path); })
        .def_property("fineCriteria",  //
                      &VNCS::Generators::Gen333::Generator::fineCriteria,
                      &VNCS::Generators::Gen333::Generator::setFineCriteria)
        .def_property(
            "visualInputMeshPath",  //
            &VNCS::Generators::Gen333::Generator::visualMeshPath,
            [](VNCS::Generators::Gen333::Generator &gen333, std::string path) { gen333.setVisualMeshPath(path); })
        .def_property("visualCriteria",  //
                      &VNCS::Generators::Gen333::Generator::visualCriteria,
                      &VNCS::Generators::Gen333::Generator::setVisualCriteria)
        .def_property(
            "blendingField",
            [](const VNCS::Generators::Gen333::Generator &mo) {},
            [](VNCS::Generators::Gen333::Generator &gen333,
               std::shared_ptr<VNCS::Generators::BlendingField<VNCS::Space3D::Real, VNCS::Space3D::Point>>
                   blendingField) { gen333.setBlendingField(blendingField); })
        .def_property(
            "coarseSamplersFilePath",
            [](const VNCS::Generators::Gen333::Generator &mo) {},
            [](VNCS::Generators::Gen333::Generator &gen333, std::string &coarseSamplersFilePath) {
                gen333.setCoarseSamplersFilePath(coarseSamplersFilePath);
            })
        .def_property(
            "fineSamplersFilePath",
            [](const VNCS::Generators::Gen333::Generator &mo) {},
            [](VNCS::Generators::Gen333::Generator &gen333, std::string &fineSamplersFilePath) {
                gen333.setFineSamplersFilePath(fineSamplersFilePath);
            })
        .def_property(
            "visualSamplersFilePath",
            [](const VNCS::Generators::Gen333::Generator &mo) {},
            [](VNCS::Generators::Gen333::Generator &gen333, std::string &visualSamplersFilePath) {
                gen333.setVisualSamplersFilePath(visualSamplersFilePath);
            })
        .def_property(
            "clusterMatrixFilePath",
            [](const VNCS::Generators::Gen333::Generator &mo) {},
            [](VNCS::Generators::Gen333::Generator &gen333, std::string &clusterMatrixFilePath) {
                gen333.setClusterMatrixFilePath(clusterMatrixFilePath);
            })
        .def_property(
            "dofFilePath",
            [](const VNCS::Generators::Gen333::Generator &mo) {},
            [](VNCS::Generators::Gen333::Generator &gen333, std::string &dofFilePath) {
                gen333.setDofFilePath(dofFilePath);
            })
        .def("exportFunction",
             [](VNCS::Generators::Gen333::Generator &gen333, pybind11::function f) {
                 gen333.setExportFunction(
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
        .def("__call__", &VNCS::Generators::Gen333::Generator::operator());

    ::py::class_<VNCS::Generators::Gen333::VisualSamplersExporter>(m, "RenderGenerator")
        .def(pybind11::init<>())
        .def_property("visualMeshPath",
                      nullptr,
                      [](VNCS::Generators::Gen333::VisualSamplersExporter &e, const std::string &path) {
                          e.setVisualMeshPath(path);
                      })
        .def_property(
            "coarseMesh",
            nullptr,
            [](VNCS::Generators::Gen333::VisualSamplersExporter &e,
               const std::pair<std::vector<std::array<VNCS::Real, 3>>, std::vector<std::array<int, 4>>> &coarseMesh) {
                std::vector<VNCS::Space3D::Point> points;
                for (const auto &p : coarseMesh.first)
                    points.emplace_back(p[0], p[1], p[2]);

                VNCS::Space3D::TetraMesh mesh;
                mesh.points = std::move(points);
                mesh.tetras = coarseMesh.second;
                e.setCoarseTetraMesh(mesh);
            })
        .def_property(
            "fineMesh",
            nullptr,
            [](VNCS::Generators::Gen333::VisualSamplersExporter &e,
               const std::pair<std::vector<std::array<VNCS::Real, 3>>, std::vector<std::array<int, 4>>> &coarseMesh) {
                std::vector<VNCS::Space3D::Point> points;
                for (const auto &p : coarseMesh.first)
                    points.emplace_back(p[0], p[1], p[2]);

                VNCS::Space3D::TetraMesh mesh;
                mesh.points = std::move(points);
                mesh.tetras = coarseMesh.second;
                e.setFineTetraMesh(mesh);
            })
        .def("matrices", &VNCS::Generators::Gen333::VisualSamplersExporter::matrices)
        .def_property("blendingField", nullptr, &VNCS::Generators::Gen333::VisualSamplersExporter::setBlendingField);
}
