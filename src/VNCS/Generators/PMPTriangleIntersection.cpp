#include "PMPTriangleIntersection.h"

#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>
#include <CGAL/polygon_mesh_processing.h>

using BarycentricCoordinates = CGAL::Barycentric_coordinates::Triangle_coordinates_2<VNCS::Space2D::K>;
namespace PMP = CGAL::Polygon_mesh_processing;

VNCS::Generators::PMPTriangleIntersection::PMPTriangleIntersection(const Mesh &mesh)
    : m_mesh(std::cref(mesh))
{
    if (m_mesh.get().num_vertices()) {
        m_point3vpm = Point3VPM(m_mesh.get());

        PMP::build_AABB_tree(m_mesh.get(), m_tree, CGAL::parameters::vertex_point_map(m_point3vpm));
    }
}

std::optional<VNCS::Generators::PMPTriangleIntersection::Mesh::Face_index>
VNCS::Generators::PMPTriangleIntersection::query(const VNCS::Generators::PMPTriangleIntersection::Point &p,
                                                 bool closest) const
{
    if (m_mesh.get().num_vertices()) {
        auto location = PMP::locate_with_AABB_tree(p, m_tree, m_mesh.get());

        if (closest)
            return location.first;

        // Build triangle
        const auto vIndices = m_mesh.get().vertices_around_face(m_mesh.get().halfedge(location.first));
        const auto v = vIndices | ranges::views::transform([this](const auto v) { return m_mesh.get().point(v); }) |
                       ranges::to_vector;

        const VNCS::Space2D::Triangle triangle(v[0], v[1], v[2]);
        if (!triangle.has_on_unbounded_side(p)) {
            return location.first;
        }
    }
    return {};
}

VNCS::Generators::PMPTriangleIntersection::Point3VPM::Point3VPM()
    : meshPtr(nullptr)
{
}

VNCS::Generators::PMPTriangleIntersection::Point3VPM::Point3VPM(const Mesh &m)
    : meshPtr(&m)
{
}

VNCS::Generators::PMPTriangleIntersection3D::PMPTriangleIntersection3D(const Mesh &mesh)
    : m_mesh(std::cref(mesh))
{
    if (m_mesh.get().num_vertices())
        PMP::build_AABB_tree(m_mesh.get(), m_tree);
}

std::optional<VNCS::Generators::PMPTriangleIntersection3D::Mesh::Face_index>
VNCS::Generators::PMPTriangleIntersection3D::query(const Point &p) const
{
    if (m_mesh.get().num_vertices()) {
        auto location = PMP::locate_with_AABB_tree(p, m_tree, m_mesh.get());
        return location.first;
    }
    return {};
}
