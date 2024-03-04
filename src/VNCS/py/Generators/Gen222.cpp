#include "Gen222.h"

#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include <VNCS/Generators/Gen222/Gen222.h>
#include "BlendingField.h"
#include "AdaptiveMeshCriteriaFunctor.h"
#include <Eigen/Dense>

#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>

namespace py = pybind11;

void VNCS::Generators::py::gen222(::py::module &m)
{
    ::py::class_<VNCS::Generators::Gen222::Generator>(m, "Gen222")
        .def(pybind11::init<>())
        .def_property(
            "coarseInputMeshPath",
            [](const VNCS::Generators::Gen222::Generator &mo) {},
            [](VNCS::Generators::Gen222::Generator &gen222, const std::string &meshPath) {
                gen222.setCoarseMeshPath(meshPath);
            })
        .def_property(
            "coarseCriteria",
            [](const VNCS::Generators::Gen222::Generator &mo) { return mo.coarseCriteria(); },
            [](VNCS::Generators::Gen222::Generator &gen222,
               std::shared_ptr<VNCS::Generators::py::AdaptiveMeshCriteriaFunctor> criteria) {
                gen222.setCoarseCriteria(criteria);
            })
        .def_property(
            "fineInputMeshPath",
            [](const VNCS::Generators::Gen222::Generator &mo) {},
            [](VNCS::Generators::Gen222::Generator &gen222, const std::string &meshPath) {
                gen222.setFineMeshPath(meshPath);
            })
        .def_property(
            "fineCriteria",
            [](const VNCS::Generators::Gen222::Generator &mo) { return mo.fineCriteria(); },
            [](VNCS::Generators::Gen222::Generator &gen222,
               std::shared_ptr<VNCS::Generators::py::AdaptiveMeshCriteriaFunctor> criteria) {
                gen222.setFineCriteria(criteria);
            })
        .def_property(
            "visualInputMeshPath",
            [](const VNCS::Generators::Gen222::Generator &mo) {},
            [](VNCS::Generators::Gen222::Generator &gen222, const std::string &meshPath) {
                gen222.setVisualMeshPath(meshPath);
            })
        .def_property(
            "visualCriteria",
            [](const VNCS::Generators::Gen222::Generator &mo) { return mo.visualCriteria(); },
            [](VNCS::Generators::Gen222::Generator &gen222,
               std::shared_ptr<VNCS::Generators::py::AdaptiveMeshCriteriaFunctor> criteria) {
                gen222.setVisualCriteria(criteria);
            })
        .def_property(
            "blendingField",
            [](const VNCS::Generators::Gen222::Generator &mo) {},
            [](VNCS::Generators::Gen222::Generator &gen222,
               std::shared_ptr<VNCS::Generators::BlendingField<VNCS::Space2D::Real, VNCS::Space2D::Point>>
                   blendingField) { gen222.setBlendingField(blendingField); })
        .def_property(
            "coarseSamplersFilePath",
            [](const VNCS::Generators::Gen222::Generator &mo) {},
            [](VNCS::Generators::Gen222::Generator &gen222, std::string &coarseSamplersFilePath) {
                gen222.setCoarseSamplersFilePath(coarseSamplersFilePath);
            })
        .def_property(
            "fineSamplersFilePath",
            [](const VNCS::Generators::Gen222::Generator &mo) {},
            [](VNCS::Generators::Gen222::Generator &gen222, std::string &fineSamplersFilePath) {
                gen222.setFineSamplersFilePath(fineSamplersFilePath);
            })
        .def_property(
            "fineNodeSamplersFilePath",
            [](const VNCS::Generators::Gen222::Generator &mo) {},
            [](VNCS::Generators::Gen222::Generator &gen222, std::string &fineNodeSamplersFilePath) {
                gen222.setFineNodeSamplersFilePath(fineNodeSamplersFilePath);
            })
        .def_property(
            "visualSamplersFilePath",
            [](const VNCS::Generators::Gen222::Generator &mo) {},
            [](VNCS::Generators::Gen222::Generator &gen222, std::string &visualSamplersFilePath) {
                gen222.setVisualSamplersFilePath(visualSamplersFilePath);
            })
        .def_property(
            "clusterMatrixFilePath",
            [](const VNCS::Generators::Gen222::Generator &mo) {},
            [](VNCS::Generators::Gen222::Generator &gen222, std::string &clusterMatrixFilePath) {
                gen222.setClusterMatrixFilePath(clusterMatrixFilePath);
            })
        .def_property(
            "dofFilePath",
            [](const VNCS::Generators::Gen222::Generator &mo) {},
            [](VNCS::Generators::Gen222::Generator &gen222, std::string &dofFilePath) {
                gen222.setDofFilePath(dofFilePath);
            })
        .def_property(
            "gridSize",
            [](const VNCS::Generators::Gen222::Generator &mo) {},
            [](VNCS::Generators::Gen222::Generator &gen222, std::size_t gridSize) { gen222.setGridSize(gridSize); })
        .def_property(
            "allowBoundary",
            [](const VNCS::Generators::Gen222::Generator &mo) { return mo.allowBoundary(); },
            [](VNCS::Generators::Gen222::Generator &gen222, bool allowBoundary) {
                gen222.setAllowBoundary(allowBoundary);
            })
        .def_property(
            "exportFunction",
            [](const VNCS::Generators::Gen222::Generator &mo) {},
            [](VNCS::Generators::Gen222::Generator &gen222, pybind11::function f) {
                gen222.setExportFunction([f](const VNCS::Space2D::Mesh &mesh,
                                             const std::vector<VNCS::Real> &b,
                                             const std::string &filename) {
                    std::vector<std::array<VNCS::Space2D::Real, 2>> points =
                        mesh.vertices() | ranges::views::transform([&mesh](const auto vId) {
                            const auto &p = mesh.point(vId);
                            return std::array<VNCS::Space2D::Real, 2>{p[0], p[1]};
                        }) |
                        ranges::to_vector;

                    const auto triangles = mesh.faces() | ranges::views::transform([&mesh](const auto fId) {
                                               const auto triangle =
                                                   mesh.vertices_around_face(mesh.halfedge(fId)) |
                                                   ranges::views::transform([](const auto vId) { return vId.idx(); }) |
                                                   ranges::to_vector;
                                               return triangle;
                                           }) |
                                           ranges::to_vector;

                    f(filename, points, triangles, b);
                });
            })
        .def("create", &VNCS::Generators::Gen222::Generator::operator());
}
