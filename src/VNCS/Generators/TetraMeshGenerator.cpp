#include "TetraMeshGenerator.h"

#include <CGAL/make_mesh_3.h>

#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>

using namespace CGAL::parameters;

VNCS::Generators::C3t3 VNCS::Generators::makeC3T3FromMesh(const VNCS::Space3D::Mesh &mesh, C3t3Criteria criteria)
{
    const auto points = mesh.vertices() |
                        ranges::views::transform([&mesh](const auto vId) { return mesh.point(vId); }) |
                        ranges::to_vector;

    const auto triangles = mesh.faces() | ranges::views::transform([&mesh](const auto fId) {
                               std::array<std::size_t, 3> vertices;
                               const auto verticesAroundFace = mesh.vertices_around_face(mesh.halfedge(fId));
                               int i = 0;
                               for (auto vertex : verticesAroundFace) {
                                   vertices[i] = vertex.idx();
                                   i++;
                               }
                               return vertices;
                           }) |
                           ranges::to_vector;

    VNCS::Generators::detail::Polyhedron polyhedron;
    namespace PMP = CGAL::Polygon_mesh_processing;
    PMP::polygon_soup_to_polygon_mesh(points, triangles, polyhedron);

    VNCS::Generators::detail::Mesh_domain domain(polyhedron);

    if (criteria.sharpAngle > 0.0)
        domain.detect_features(criteria.sharpAngle);

    VNCS::Generators::detail::Mesh_criteria meshCriteria(facet_angle = criteria.facet_angle,
                                                         facet_size = criteria.facet_size,
                                                         facet_distance = criteria.facet_distance,
                                                         cell_radius_edge_ratio = criteria.cell_radius_edge_ratio,
                                                         cell_size = criteria.cell_size,
                                                         edge_size = criteria.sharpEdgeSize);

    VNCS::Generators::C3t3 c3t3 =
        CGAL::make_mesh_3<VNCS::Generators::C3t3>(domain, meshCriteria, no_perturb(), no_exude());

    return c3t3;
}

std::pair<std::vector<VNCS::Space3D::Point>, std::vector<std::array<int, 4>>> VNCS::Generators::c3t3ToTetrahedra(
    const VNCS::Generators::C3t3 &c3t3)
{
    const VNCS::Generators::detail::Tr &tr = c3t3.triangulation();

    std::map<VNCS::Generators::detail::Tr::Vertex_handle, int> Vnbe;

    for (auto cit = c3t3.cells_begin(); cit != c3t3.cells_end(); ++cit) {
        for (int i = 0; i < 4; i++)
            ++Vnbe[cit->vertex(i)];
    }

    std::map<VNCS::Generators::detail::Tr::Vertex_handle, int> V;
    std::vector<VNCS::Space3D::Point> points;
    auto inum = 0;
    for (auto vit = tr.finite_vertices_begin(); vit != tr.finite_vertices_end(); ++vit) {
        if (!(Vnbe.find(vit) == Vnbe.end() || Vnbe[vit] <= 0)) {
            const auto &pointCgal = vit->point();
            points.emplace_back(CGAL::to_double(pointCgal.x()),  //
                                CGAL::to_double(pointCgal.y()),
                                CGAL::to_double(pointCgal.z()));
            V[vit] = inum++;
        }
    }
    std::vector<std::array<int, 4>> tetrahedra;
    for (auto cit = c3t3.cells_begin(); cit != c3t3.cells_end(); ++cit) {
        std::array<int, 4> tetra;
        for (int i = 0; i < 4; i++)
            tetra[i] = V[cit->vertex(i)];
        tetrahedra.push_back(tetra);
    }

    return std::make_pair(points, tetrahedra);
}
