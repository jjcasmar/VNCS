#ifndef VNCS_TETRAAABBTREE_H
#define VNCS_TETRAAABBTREE_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <VNCS/Spaces.h>

namespace VNCS
{
namespace detail
{
class AABBTetrahedronPrimitive
{
public:
    using Id = std::vector<VNCS::Space3D::Tetra>::const_iterator;
    using Point = const VNCS::Space3D::Point;
    using Datum = const VNCS::Space3D::Tetra;

private:
    Id m_tetrahedron;

public:
    AABBTetrahedronPrimitive() = default;
    AABBTetrahedronPrimitive(Id id);

    const Id &id() const;
    const Datum &datum() const;
    Point reference_point() const;
};

using AABBTraitsTetrahedron = CGAL::AABB_traits<VNCS::Space2D::K, AABBTetrahedronPrimitive>;
}  // namespace detail

using TetraAABBTree = CGAL::AABB_tree<detail::AABBTraitsTetrahedron>;
}  // namespace VNCS

#endif  // TETRAAABBTREE_H
