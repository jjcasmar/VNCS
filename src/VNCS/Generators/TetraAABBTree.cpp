#include "TetraAABBTree.h"

VNCS::detail::AABBTetrahedronPrimitive::AABBTetrahedronPrimitive(VNCS::detail::AABBTetrahedronPrimitive::Id id)
    : m_tetrahedron(id)
{
}

const VNCS::detail::AABBTetrahedronPrimitive::Id &VNCS::detail::AABBTetrahedronPrimitive::id() const
{
    return m_tetrahedron;
}

const VNCS::detail::AABBTetrahedronPrimitive::Datum &VNCS::detail::AABBTetrahedronPrimitive::datum() const
{
    return *m_tetrahedron;
}

VNCS::detail::AABBTetrahedronPrimitive::Point VNCS::detail::AABBTetrahedronPrimitive::reference_point() const
{
    return (*m_tetrahedron)[0];
}
