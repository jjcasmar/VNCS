#ifndef VNCS_GENERATORS_PMPTRIANGLEINTERSECTION_H
#define VNCS_GENERATORS_PMPTRIANGLEINTERSECTION_H

#include <VNCS/Spaces.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>

#include <range/v3/algorithm/all_of.hpp>

namespace VNCS
{
namespace Generators
{
class PMPTriangleIntersection
{
public:
    using Point = VNCS::Space2D::Point;
    using Point3 = VNCS::Space2D::K::Point_3;
    using Mesh = VNCS::Space2D::Mesh;
    using Real = VNCS::Space2D::Real;

    struct Point3VPM {
        using key_type = Mesh::Vertex_index;
        using value_type = Point3;
        using reference = value_type;
        using category = boost::readable_property_map_tag;

        Point3VPM();

        Point3VPM(const Mesh &m);

        friend Point3VPM::value_type get(const Point3VPM &map, Mesh::Vertex_index idx)
        {
            Point p = map.meshPtr->point(idx);
            return {p[0], p[1], 0};
        }

        const Mesh *meshPtr;
    };

    using AABBTreeTraits = CGAL::AABB_traits<VNCS::Space2D::K,  //
                                             CGAL::AABB_face_graph_triangle_primitive<Mesh, Point3VPM>>;
    using AABBTree = CGAL::AABB_tree<AABBTreeTraits>;

    PMPTriangleIntersection(const Mesh &mesh);

    std::optional<Mesh::Face_index> query(const Point &p, bool closest = false) const;

private:
    std::reference_wrapper<const Mesh> m_mesh;
    Point3VPM m_point3vpm;
    AABBTree m_tree;
};

class PMPTriangleIntersection3D
{
public:
    using Point = VNCS::Space3D::Point;
    using Point3 = VNCS::Space3D::K::Point_3;
    using Mesh = VNCS::Space3D::Mesh;
    using Real = VNCS::Space3D::Real;

    using AABBTreeTraits = CGAL::AABB_traits<VNCS::Space3D::K,  //
                                             CGAL::AABB_face_graph_triangle_primitive<Mesh>>;
    using AABBTree = CGAL::AABB_tree<AABBTreeTraits>;

    PMPTriangleIntersection3D(const Mesh &mesh);

    std::optional<Mesh::Face_index> query(const Point &p) const;

private:
    std::reference_wrapper<const Mesh> m_mesh;
    AABBTree m_tree;
};
}  // namespace Generators
}  // namespace VNCS

#endif  // PMPTRIANGLEINTERSECTION_H
