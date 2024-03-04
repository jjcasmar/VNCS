#include "CDT.h"

VNCS::Space2D::Mesh VNCS::Generators::cdt(const std::vector<VNCS::Space2D::Point> &points,
                                          const std::vector<std::array<std::size_t, 2>> &edges,
                                          VNCS::Generators::AdaptiveMeshCriteria<VNCS::Generators::CDT> criteria)
{
    VNCS::Generators::CDT cdt;
    std::vector<VNCS::Generators::CDT::Vertex_handle> vertexHandles;
    std::transform(std::begin(points), std::end(points), std::back_inserter(vertexHandles), [&cdt](const auto &point) {
        return cdt.insert(VNCS::Generators::CDT::Point(point[0], point[1]));
    });

    for (const auto &edge : edges) {
        const auto &va = vertexHandles[edge[0]];
        const auto &vb = vertexHandles[edge[1]];
        cdt.insert_constraint(va, vb);
    }

    CGAL::refine_Delaunay_mesh_2(cdt, criteria);
    //    std::cout << cdt.number_of_vertices() << "\n";
    //    std::cout << cdt.number_of_faces() << "\n";

    namespace PMP = CGAL::Polygon_mesh_processing;
    VNCS::Space2D::Mesh mesh;

    const auto vIt = cdt.vertices_begin();
    const auto vEnd = cdt.vertices_end();
    std::unordered_map<VNCS::Generators::CDT::Vertex_handle, VNCS::Space2D::Mesh::Vertex_index> vertexHandlesMap;
    size_t index = 0;
    for (auto it = vIt; it != vEnd; ++it) {
        const auto handle = mesh.add_vertex({it->point()[0], it->point()[1]});
        vertexHandlesMap[it->handle()] = handle;
        index++;
    }

    const auto fIt = cdt.faces_begin();
    const auto fEnd = cdt.faces_end();
    std::vector<size_t> refineFacets;
    for (auto it = fIt; it != fEnd; ++it) {
        if (it->is_in_domain()) {
            mesh.add_face(vertexHandlesMap[it->vertex(0)->handle()],
                          vertexHandlesMap[it->vertex(1)->handle()],
                          vertexHandlesMap[it->vertex(2)->handle()]);
        }
    }

    return mesh;
}

VNCS::Space2D::Mesh VNCS::Generators::cdt(const VNCS::Space2D::Mesh &mesh,
                                          VNCS::Generators::AdaptiveMeshCriteria<VNCS::Generators::CDT> criteria)
{
    // Extract boundaries that will guide the CDT refinement
    std::vector<VNCS::Space2D::Point> boundaryVertices;
    std::unordered_map<int, std::size_t>
        boundaryVerticesIndices;  // Given the index in the original mesh, returns the index in the list of boundaries
    for (const auto vertexHandle : mesh.vertices())
        if (mesh.is_border(vertexHandle)) {
            boundaryVerticesIndices[vertexHandle.idx()] = boundaryVertices.size();
            boundaryVertices.push_back(mesh.point(vertexHandle));
        }

    std::vector<std::array<std::size_t, 2>> boundaryEdges;
    for (const auto edgeHandle : mesh.edges())
        if (mesh.is_border(edgeHandle))
            boundaryEdges.push_back({boundaryVerticesIndices[mesh.vertex(edgeHandle, 0)],
                                     boundaryVerticesIndices[mesh.vertex(edgeHandle, 1)]});

    std::cout << "Vertices: " << boundaryVertices.size() << "\n";
    std::cout << "Edges: " << boundaryEdges.size() << "\n";
    return cdt(boundaryVertices, boundaryEdges, criteria);
}
