#ifndef VNCS_GENERATORS_CDT_H
#define VNCS_GENERATORS_CDT_H

#include <VNCS/Spaces.h>
#include <VNCS/Generators/AdaptiveMeshCriteria.h>

#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

namespace VNCS
{
namespace Generators
{
namespace detail_impl
{
using K = VNCS::Space2D::K;
using Vb = CGAL::Triangulation_vertex_base_2<K>;
using Fb = CGAL::Delaunay_mesh_face_base_2<K>;
using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
}  // namespace detail_impl

using CDT = CGAL::Constrained_Delaunay_triangulation_2<detail_impl::K, detail_impl::Tds>;

VNCS::Space2D::Mesh cdt(const std::vector<VNCS::Space2D::Point> &points,
                        const std::vector<std::array<std::size_t, 2>> &edges,
                        VNCS::Generators::AdaptiveMeshCriteria<CDT> criteria);

VNCS::Space2D::Mesh cdt(const VNCS::Space2D::Mesh &mesh, VNCS::Generators::AdaptiveMeshCriteria<CDT> criteria);

template <typename Tr>
VNCS::Space2D::Mesh meshFromTriangulation(const Tr &tr)
{
    VNCS::Space2D::Mesh mesh;

    const auto vIt = tr.vertices_begin();
    const auto vEnd = tr.vertices_end();
    std::unordered_map<typename Tr::Vertex_handle, VNCS::Space2D::Mesh::Vertex_index> vertexHandlesMap;
    size_t index = 0;
    for (auto it = vIt; it != vEnd; ++it) {
        const auto handle = mesh.add_vertex({it->point()[0], it->point()[1]});
        vertexHandlesMap[it->handle()] = handle;
        index++;
    }

    const auto fIt = tr.faces_begin();
    const auto fEnd = tr.faces_end();
    std::vector<size_t> refineFacets;
    for (auto it = fIt; it != fEnd; ++it) {
        mesh.add_face(vertexHandlesMap[it->vertex(0)->handle()],
                      vertexHandlesMap[it->vertex(1)->handle()],
                      vertexHandlesMap[it->vertex(2)->handle()]);
    }

    return mesh;
}
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_GENERATORS_CDT_H
