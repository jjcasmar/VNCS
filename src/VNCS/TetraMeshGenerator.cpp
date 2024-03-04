#include "TetraMeshGenerator.h"

#include <sofa/core/topology/Topology.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include <CPFSofaPlugin/DataExtensions.h>

using SurfaceMesh = CGAL::Surface_mesh<CPF::Point_3>;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;
#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

#include <sofa/core/ObjectFactory.h>

static int register_TetraMeshGenerator =
    sofa::core::RegisterObject("TetraMeshGenerator").add<CPF::TetraMeshGenerator>();

CPF::TetraMeshGenerator::TetraMeshGenerator()
    : m_triangleTopologyLink(initLink("topology", "Topology"))
    , m_outputPoints(initData(&m_outputPoints, "outputPoints", "Output points"))
    , m_outputTetras(initData(&m_outputTetras, "outputTetras", "Output tetras"))
{
}

void CPF::TetraMeshGenerator::init()
{
    // Create the input domain
    const ::sofa::component::topology::TriangleSetTopologyContainer &triangleSet = *m_triangleTopologyLink;
    SurfaceMesh inputMeshDomain;
    const auto &points = make_read_accessor(triangleSet.d_initPoints);
    for (const auto &p : points)
        inputMeshDomain.add_vertex(CPF::Point_3(p[0], p[1], p[2]));

    const auto &triangles = make_read_accessor(triangleSet.d_triangle);
    for (const auto &triangle : triangles)
        inputMeshDomain.add_face(SurfaceMesh::Vertex_index(triangle[0]),  //
                                 SurfaceMesh::Vertex_index(triangle[1]),
                                 SurfaceMesh::Vertex_index(triangle[2]));

    Polyhedron polyhedron;
    CGAL::copy_face_graph(inputMeshDomain, polyhedron);
    Mesh_domain domain(polyhedron);
    Mesh_criteria::Cell_criteria cellCriteria(1, 0.1);
    Mesh_criteria::Facet_criteria facetCriteria(30, 0.1, 0.25);
    Mesh_criteria criteria(facetCriteria, cellCriteria);

    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

    const auto &tr = c3t3.triangulation();
    std::unordered_map<Tr::Vertex_handle, CPF::SofaTypes::Tetra::value_type> vertexHandleToIndex;
    auto &outputPoints = make_write_accessor(m_outputPoints);
    int i = 0;
    for (const auto vertexHandle : tr.finite_vertex_handles()) {
        const auto &p = vertexHandle->point();
        outputPoints.push_back(CPF::SofaTypes::Point{p[0], p[1], p[2]});
        vertexHandleToIndex[vertexHandle] = i;
        ++i;
    }

    auto &tetras = make_write_accessor(m_outputTetras);
    for (const auto cellHandle : tr.finite_cell_handles()) {
        const auto vA_handle = cellHandle->vertex(0);
        const auto vB_handle = cellHandle->vertex(1);
        const auto vC_handle = cellHandle->vertex(2);
        const auto vD_handle = cellHandle->vertex(3);

        tetras.push_back(CPF::SofaTypes::Tetra{
            vertexHandleToIndex[vA_handle],
            vertexHandleToIndex[vB_handle],
            vertexHandleToIndex[vC_handle],
            vertexHandleToIndex[vD_handle],
        });
    }
}
