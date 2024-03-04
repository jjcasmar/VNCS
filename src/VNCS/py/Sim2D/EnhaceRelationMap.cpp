#include "EnhaceRelationMap.h"

#include <pybind11/stl.h>

#include <VNCS/Generators/PMPTriangleIntersection.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/IO/print_wavefront.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>

namespace VNCS
{
namespace Sim2D
{
namespace py
{
using TriangleCoordinates = CGAL::Barycentric_coordinates::Triangle_coordinates_2<VNCS::Space2D::K>;

std::vector<std::pair<int, std::array<VNCS::Space2D::Real, 3>>> computeEnhaceMap(
    const std::vector<VNCS::Space2D::Point> &meshPoints,
    const std::vector<std::array<int, 3>> &meshTriangles,
    const std::vector<VNCS::Space2D::Point> &newSimulationPoints)
{
    namespace PMP = CGAL::Polygon_mesh_processing;

    VNCS::Space2D::Mesh mesh;
    PMP::polygon_soup_to_polygon_mesh(meshPoints, meshTriangles, mesh);

    VNCS::Generators::PMPTriangleIntersection intersector(mesh);

    std::vector<std::pair<int, std::array<VNCS::Space2D::Real, 3>>> map;
    for (const auto &point : newSimulationPoints) {
        const auto closestTriangle = intersector.query(point, true);
        if (closestTriangle) {
            const auto points = mesh.vertices_around_face(mesh.halfedge(closestTriangle.value())) |
                                ranges::views::transform([&mesh](const auto vId) { return mesh.point(vId); }) |
                                ranges::to_vector;

            TriangleCoordinates triangleCoordinates(points[0], points[1], points[2]);

            std::array<VNCS::Space2D::Real, 3> coordinates;
            triangleCoordinates(point, std::begin(coordinates));

            map.push_back(std::make_pair(closestTriangle->idx(), coordinates));
        }
    }

    return map;
}

void enhaceMap(pybind11::module &m)
{
    m.def("enhaceMap",
          [](const std::vector<std::array<VNCS::Space2D::Real, 2>> &meshPoints,
             const std::vector<std::array<int, 3>> &meshTriangles,
             const std::vector<std::array<VNCS::Space2D::Real, 2>> &newSimulationPoints) {
              const auto meshPointsCGAL =
                  meshPoints |
                  ranges::views::transform([](const auto &p) { return VNCS::Space2D::Point(p[0], p[1]); }) |
                  ranges::to_vector;

              const auto newSimulationPointsCGAL =
                  newSimulationPoints |
                  ranges::views::transform([](const auto &p) { return VNCS::Space2D::Point(p[0], p[1]); }) |
                  ranges::to_vector;

              return computeEnhaceMap(meshPointsCGAL, meshTriangles, newSimulationPointsCGAL);
          });
}

}  // namespace py
}  // namespace Sim2D
}  // namespace VNCS
