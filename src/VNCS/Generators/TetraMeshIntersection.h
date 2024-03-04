#ifndef VNCS_GENERATORS_TETRAMESHINTERSECTION_H
#define VNCS_GENERATORS_TETRAMESHINTERSECTION_H

#include <VNCS/Spaces.h>
#include <VNCS/Generators/TetraAABBTree.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

namespace VNCS
{
namespace Generators
{
class TetraMeshIntersection
{
    using TriangleId = std::vector<VNCS::Space3D::Triangle>::const_iterator;
    using TrianglePrimitive = CGAL::AABB_triangle_primitive<VNCS::Space3D::K, TriangleId>;
    using AABBTriangleTraits = CGAL::AABB_traits<VNCS::Space3D::K, TrianglePrimitive>;
    using TriangleAABBTree = CGAL::AABB_tree<AABBTriangleTraits>;

public:
    using Point = VNCS::Space3D::Point;
    using Mesh = VNCS::Space3D::TetraMesh;
    using Real = VNCS::Space3D::Real;

    TetraMeshIntersection(const VNCS::Space3D::TetraMesh &tetraMesh);

    bool doIntersect(const CGAL::Bbox_3 &bbox) const;
    std::optional<int> query(const Point &p) const;
    std::optional<int> query(const Point &p, bool closest) const;

    const std::vector<VNCS::Space3D::Tetra> &tetras() const;
    const TetraAABBTree &tetraAABBTree() const;

private:
    std::reference_wrapper<const VNCS::Space3D::TetraMesh> m_mesh;
    std::vector<VNCS::Space3D::Tetra> m_tetrahedra;
    std::vector<VNCS::Space3D::Triangle> m_boundaryTriangles;
    std::vector<int> m_boundaryTrianglesTetraId;

    TetraAABBTree m_tetraAABBTree;
    TriangleAABBTree m_triangleAABBTree;
};
}  // namespace Generators
}  // namespace VNCS

#endif  //  VNCS_GENERATORS_TETRAMESHINTERSECTION_H
