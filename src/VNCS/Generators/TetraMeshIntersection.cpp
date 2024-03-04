#include "TetraMeshIntersection.h"

#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <algorithm>

#include <spdlog/spdlog.h>

namespace
{
struct STriangleLessComparator {
    bool operator()(const std::array<int, 3> &a, const std::array<int, 3> &b) const
    {
        auto a_ = a;
        auto b_ = b;
        std::sort(std::begin(a_), std::end(a_));
        std::sort(std::begin(b_), std::end(b_));
        return a_ < b_;
    }
};

struct STriangleEqualComparator {
    bool operator()(const std::pair<std::array<int, 3>, int> &a, const std::pair<std::array<int, 3>, int> &b) const
    {
        auto a_ = a.first;
        auto b_ = b.first;
        std::sort(std::begin(a_), std::end(a_));
        std::sort(std::begin(b_), std::end(b_));
        return a_[0] == b_[0] &&  //
               a_[1] == b_[1] &&  //
               a_[2] == b_[2];
    }
};
}  // namespace

VNCS::Generators::TetraMeshIntersection::TetraMeshIntersection(const VNCS::Space3D::TetraMesh &tetraMesh)
    : m_mesh(tetraMesh)
{
    const auto &tetras = tetraMesh.tetras;
    m_tetrahedra = tetras | ranges::views::transform([&tetraMesh](const auto &tetra) {
                       return VNCS::Space3D::Tetra{tetraMesh.points[tetra[0]],
                                                   tetraMesh.points[tetra[1]],
                                                   tetraMesh.points[tetra[2]],
                                                   tetraMesh.points[tetra[3]]};
                   }) |
                   ranges::to_vector;

    m_tetraAABBTree.rebuild(std::begin(m_tetrahedra), std::end(m_tetrahedra));

    {
        std::vector<std::pair<std::array<int, 3>, int>> tetraTrianglesAndTetraId;
        for (const auto [tetraId, tetra] : ranges::views::enumerate(tetraMesh.tetras)) {
            tetraTrianglesAndTetraId.emplace_back(std::make_pair<std::array<int, 3>, int>({tetra[0],  //
                                                                                           tetra[1],
                                                                                           tetra[2]},
                                                                                          int(tetraId)));
            tetraTrianglesAndTetraId.emplace_back(std::make_pair<std::array<int, 3>, int>({tetra[1],  //
                                                                                           tetra[2],
                                                                                           tetra[3]},
                                                                                          int(tetraId)));
            tetraTrianglesAndTetraId.emplace_back(std::make_pair<std::array<int, 3>, int>({tetra[0],  //
                                                                                           tetra[2],
                                                                                           tetra[3]},
                                                                                          int(tetraId)));
            tetraTrianglesAndTetraId.emplace_back(std::make_pair<std::array<int, 3>, int>({tetra[0],  //
                                                                                           tetra[1],
                                                                                           tetra[3]},
                                                                                          int(tetraId)));
        }

        // Order the triangles and remove repeated ones, no matter the order of the indices: (a,b,c) is the same
        // triangle as (b,a,c)
        ranges::sort(std::begin(tetraTrianglesAndTetraId),
                     std::end(tetraTrianglesAndTetraId),
                     STriangleLessComparator{},
                     &std::pair<std::array<int, 3>, int>::first);
        for (auto it = std::adjacent_find(
                 tetraTrianglesAndTetraId.begin(), tetraTrianglesAndTetraId.end(), STriangleEqualComparator{});  //
             it != tetraTrianglesAndTetraId.end();                                                               //
             it = std::adjacent_find(
                 tetraTrianglesAndTetraId.begin(), tetraTrianglesAndTetraId.end(), STriangleEqualComparator{}))  //
            tetraTrianglesAndTetraId.erase(it, it + 2);

        // The remaining triangles are the boundary
        for (const auto &boundaryTriangle : tetraTrianglesAndTetraId) {
            m_boundaryTriangles.emplace_back(tetraMesh.points[boundaryTriangle.first[0]],
                                             tetraMesh.points[boundaryTriangle.first[1]],
                                             tetraMesh.points[boundaryTriangle.first[2]]);
            m_boundaryTrianglesTetraId.push_back(boundaryTriangle.second);
        }

        spdlog::get("Gen333")->info("Accelerate distance queries()");
        m_triangleAABBTree.rebuild(m_boundaryTriangles.begin(), m_boundaryTriangles.end());
        m_triangleAABBTree.accelerate_distance_queries();
    }
}

bool VNCS::Generators::TetraMeshIntersection::doIntersect(const CGAL::Bbox_3 &bbox) const
{
    return m_tetraAABBTree.do_intersect(bbox);
}

std::optional<int> VNCS::Generators::TetraMeshIntersection::query(
    const VNCS::Generators::TetraMeshIntersection::Point &p) const
{
    return query(p, false);
}

std::optional<int> VNCS::Generators::TetraMeshIntersection::query(
    const VNCS::Generators::TetraMeshIntersection::Point &p,
    bool closest) const
{
    const auto intersection = m_tetraAABBTree.any_intersected_primitive(p);
    if (intersection)
        return {std::distance(std::begin(m_tetrahedra), intersection.value())};

    if (closest) {
        auto closestPrimitive = m_triangleAABBTree.closest_point_and_primitive(p);
        int index = std::distance(m_boundaryTriangles.begin(), closestPrimitive.second);
        return m_boundaryTrianglesTetraId[index];
    }
    return {};
}

const std::vector<VNCS::Space3D::Tetra> &VNCS::Generators::TetraMeshIntersection::tetras() const
{
    return m_tetrahedra;
}

const VNCS::TetraAABBTree &VNCS::Generators::TetraMeshIntersection::tetraAABBTree() const
{
    return m_tetraAABBTree;
}
