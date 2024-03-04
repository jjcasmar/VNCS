#include "Gen222.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/IO/print_wavefront.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <range/v3/to_container.hpp>
#include <range/v3/view/subrange.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/action/sort.hpp>
#include <range/v3/view/group_by.hpp>
#include <range/v3/numeric/accumulate.hpp>

#include <unsupported/Eigen/SparseExtra>

#include <fstream>

#include <VNCS/Generators/CDT.h>
#include <VNCS/Generators/PMPTriangleIntersection.h>
#include <VNCS/Generators/Gen222/SamplingPointsExporter.h>
#include <VNCS/Generators/Gen222/VisualSamplersExporter.h>

namespace
{
class My_point_property_map
{
    const std::vector<VNCS::Space2D::Point> &points;

public:
    typedef VNCS::Space2D::Point value_type;
    typedef const value_type &reference;
    typedef std::size_t key_type;
    typedef boost::lvalue_property_map_tag category;
    My_point_property_map(const std::vector<value_type> &pts)
        : points(pts)
    {
    }
    reference operator[](key_type k) const { return points[k]; }
    friend reference get(const My_point_property_map &ppmap, key_type i) { return ppmap[i]; }
};

typedef CGAL::Search_traits_2<VNCS::Space2D::K> TreeTraitsBase;
typedef CGAL::Search_traits_adapter<std::size_t, My_point_property_map, TreeTraitsBase> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;
typedef Tree::Splitter Splitter;
typedef Neighbor_search::Distance Distance;

using Mesh = VNCS::Space2D::Mesh;

using Point = VNCS::Space2D::Point;
using Triangle = VNCS::Space2D::Triangle;
using Real = VNCS::Space2D::Real;

}  // namespace

void VNCS::Generators::Gen222::Generator::operator()()
{
    // Load OBJ

    Mesh coarseInputMesh;
    Mesh fineInputMesh;
    Mesh visualInputMesh;

    {
        std::vector<VNCS::Space2D::Point> coarsePoints;
        std::vector<std::vector<std::size_t>> coarseFaces;
        std::ifstream coarseMeshIn(m_coarseMeshPath);
        CGAL::read_OBJ(coarseMeshIn, coarsePoints, coarseFaces);

        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::polygon_soup_to_polygon_mesh(coarsePoints, coarseFaces, coarseInputMesh);
    }

    {
        std::vector<VNCS::Space2D::Point> coarsePoints;
        std::vector<std::vector<std::size_t>> coarseFaces;
        std::ifstream coarseMeshIn(m_fineMeshPath);
        CGAL::read_OBJ(coarseMeshIn, coarsePoints, coarseFaces);

        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::polygon_soup_to_polygon_mesh(coarsePoints, coarseFaces, fineInputMesh);
    }

    {
        std::vector<VNCS::Space2D::Point> coarsePoints;
        std::vector<std::vector<std::size_t>> coarseFaces;
        std::ifstream coarseMeshIn(m_visualMeshPath);
        CGAL::read_OBJ(coarseMeshIn, coarsePoints, coarseFaces);

        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::polygon_soup_to_polygon_mesh(coarsePoints, coarseFaces, visualInputMesh);
    }

    const auto coarseMesh = VNCS::Generators::cdt(
        coarseInputMesh, VNCS::Generators::AdaptiveMeshCriteria<VNCS::Generators::CDT>{m_coarseCriteria});
    const auto fineMesh = VNCS::Generators::cdt(
        fineInputMesh, VNCS::Generators::AdaptiveMeshCriteria<VNCS::Generators::CDT>{m_fineCriteria});
    const auto visualMesh = VNCS::Generators::cdt(
        visualInputMesh, VNCS::Generators::AdaptiveMeshCriteria<VNCS::Generators::CDT>{m_visualCriteria});

    {
        CGAL::Polyhedron_3<VNCS::Space2D::K> visualOutputMesh;

        std::vector<VNCS::Space2D::K::Point_3> points = visualMesh.vertices() |
                                                        ranges::views::transform([visualMesh](const auto &v) {
                                                            const auto p = visualMesh.point(v);
                                                            return VNCS::Space2D::K::Point_3{p.x(), p.y(), 0};
                                                        }) |
                                                        ranges::to_vector;
        std::vector<std::vector<std::size_t>> triangles;
        for (const auto face : visualMesh.faces()) {
            std::vector<std::size_t> triangle;
            for (const auto vd : visualMesh.vertices_around_face(visualMesh.halfedge(face))) {
                triangle.push_back(vd);
            }
            triangles.push_back(triangle);
        }

        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::polygon_soup_to_polygon_mesh(points, triangles, visualOutputMesh);

        std::ofstream outVisual("visual.obj");
        CGAL::print_polyhedron_wavefront(outVisual, visualOutputMesh);
    }

    enum class XN_TYPE { CoarseRegion, Hybrid_Boundary, Hybrid_RegionBoundary, xC2xN };

    // Separate the vertices into three sets:
    // 1) xN vertices
    //    - Vertices that have a coarse value of 1 (this are vertices in the coarse region)
    //    - Vertices that are a boundary and have at least one adjacent vertex with blending value different than 0
    //    (in order to have a tetrahedron in the hybrid region)
    //    - Vertices that have a blending value of 0 and have at least one adjacent vertex with blending value
    //    different than 0 (in order to have a tetrahedron in the transition between hybrid and fine regions)
    // 2) xC vertices
    //    - Vertices that have a blending between 0 or 1 and are not boundary
    // 3) Rest of vertices
    //    - They dont participate in the simulation

    std::vector<unsigned int> xNIndices;
    std::vector<unsigned int> xCIndices;
    std::vector<VNCS::Space2D::Point> xC;

    std::vector<XN_TYPE> xNVerticesType;
    {
        for (const Mesh::Vertex_index vd : coarseMesh.vertices()) {
            const auto p = coarseMesh.point(vd);
            if (m_blendingField->blending(p) == 1.0) {
                xNIndices.push_back(vd.idx());
            } else {
                if (!m_allowBoundary) {
                    if (coarseMesh.is_border(vd)) {
                        // Circulate on the vertices of vd and check if there is at least one with blending != 0.0
                        for (const auto vdCirculator : coarseMesh.vertices_around_target(coarseMesh.halfedge(vd))) {
                            if (m_blendingField->blending(coarseMesh.point(Mesh::Vertex_index(vdCirculator))) != 0.0) {
                                xNIndices.push_back(vd.idx());
                                break;
                            }
                        }
                    } else {
                        if (m_blendingField->blending(p) == 0.0) {
                            for (const auto vdCirculator : coarseMesh.vertices_around_target(coarseMesh.halfedge(vd))) {
                                if (m_blendingField->blending(coarseMesh.point(vdCirculator)) != 0.0) {
                                    xNIndices.push_back(vd.idx());
                                    break;
                                }
                            }
                        } else {
                            xCIndices.push_back(vd.idx());
                            xC.push_back(p);
                        }
                    }
                } else {
                    if (m_blendingField->blending(p) == 0.0) {
                        for (const auto vdCirculator : coarseMesh.vertices_around_target(coarseMesh.halfedge(vd))) {
                            if (m_blendingField->blending(coarseMesh.point(vdCirculator)) != 0.0) {
                                xNIndices.push_back(vd.idx());
                                break;
                            }
                        }
                    } else {
                        xCIndices.push_back(vd.idx());
                        xC.push_back(p);
                    }
                }
            }
        }
    }

    // xC degrees of freedom will be initialized later as we need to compute the cluster point

    // Fine tetrahedra appears in topology if not all its vertices have blending == 1.0
    // The fine points can be in two sets, pS and d.
    // They belong to pS if the blending value
    // They belong to d if the blending value E (0,1]
    // They belong to pS if the blending value is 0
    // If the blending value == 255 -> Full -> pS
    // If the blending value E (0,1) -> Hybrid -> d
    // If the blending value == 0 -> Coarse
    // Coarse nodes are not simulated

    std::vector<unsigned int> dIndices;
    std::vector<unsigned int> pSIndices;

    for (const auto vd : fineMesh.vertices()) {
        const auto &p = fineMesh.point(vd);

        const auto blending = m_blendingField->blending(p);
        if (blending > 0.0 && blending < 1.0)
            dIndices.push_back(vd.idx());

        if (blending == 0.0)
            pSIndices.push_back(vd.idx());

        if (blending == 1.0) {
            for (const auto vdCirculator : fineMesh.vertices_around_target(fineMesh.halfedge(vd))) {
                if (m_blendingField->blending(fineMesh.point(vdCirculator)) != 1.0) {
                    pSIndices.push_back(vd.idx());
                    break;
                }
            }
        }
    }

    std::vector<VNCS::Space2D::Point> d;
    d.reserve(dIndices.size());
    for (auto dIndex : dIndices) {
        const auto p = fineMesh.point(Mesh::Vertex_index(dIndex));
        d.push_back(p);
    }

    std::vector<VNCS::Space2D::Point> pS;
    pS.reserve(pSIndices.size());
    for (auto pSIndex : pSIndices) {
        const auto p = fineMesh.point(Mesh::Vertex_index(pSIndex));
        pS.push_back(p);
    }

    // We have to find the closest xC vertex to each d vertex
    // If any xC vertex is not selected, we will crash later!
    My_point_property_map ppmap(xC);
    Tree searchTree(boost::counting_iterator<std::size_t>(0),
                    boost::counting_iterator<std::size_t>(xC.size()),
                    Splitter(),
                    TreeTraits(ppmap));
    Distance trDistance(ppmap);

    std::vector<unsigned int> assignedXCvertex(dIndices.size());
    for (int i = 0; i < dIndices.size(); ++i) {
        const auto &dPoint = fineMesh.point(Mesh::Vertex_index(dIndices[i]));
        Neighbor_search search(searchTree, dPoint, 1, 0.0001, true, trDistance);
        for (const auto &searchResult : search) {
            assignedXCvertex[i] = xCIndices[searchResult.first];
        }
    }

    std::vector<Eigen::Triplet<VNCS::Space2D::Real>> triplets;
    triplets.reserve(3 * dIndices.size());
    for (int i = 0; i < dIndices.size(); ++i) {
        auto it = std::find(std::begin(xCIndices), std::end(xCIndices), assignedXCvertex[i]);
        auto index = std::distance(std::begin(xCIndices), it);
        triplets.push_back(Eigen::Triplet<VNCS::Space2D::Real>(2 * i + 0, 2 * index + 0, 1.0));
        triplets.push_back(Eigen::Triplet<VNCS::Space2D::Real>(2 * i + 1, 2 * index + 1, 1.0));
    }

    // Create xN degrees of freedom
    std::vector<VNCS::Space2D::Point> xN;
    xN.reserve(xNIndices.size());
    for (auto xNVertexId : xNIndices) {
        const auto p = coarseMesh.point(Mesh::Vertex_index(xNVertexId));
        xN.push_back({p.x(), p.y()});
    }

    Eigen::SparseMatrix<VNCS::Space2D::Real> C;
    C.resize(2 * dIndices.size(),  //
             2 * xCIndices.size());
    C.setFromTriplets(std::begin(triplets), std::end(triplets), [](const auto &a, const auto &b) {
        spdlog::get("VNCS")->info("Repeated triplets shouldnt happen");
        return a;
    });

    Eigen::saveMarket(C, m_clusterMatrixFilePath);

    Eigen::SparseMatrix<VNCS::Space2D::Real> CtC = (C.transpose() * C);
    Eigen::DiagonalMatrix<VNCS::Space2D::Real, Eigen::Dynamic> CtCInverse;
    CtCInverse.resize(CtC.rows());
    for (int i = 0; i < CtC.rows(); ++i) {
        if (CtC.coeff(i, i) == VNCS::Space2D::Real{0.0})
            spdlog::get("VNCS")->error("Coarse node {} doesnt have associated d nodes", xN.size() + i / 2);
        CtCInverse.diagonal()[i] = 1.0 / CtC.coeff(i, i);
    }

    xC.resize(xCIndices.size());

    Eigen::Matrix<VNCS::Space2D::Real, Eigen::Dynamic, 1> dEig;
    dEig.resize(2 * d.size());
    for (const auto &[index, p] : ranges::views::enumerate(d)) {
        dEig[2 * index + 0] = p[0];
        dEig[2 * index + 1] = p[1];
    }

    Eigen::Matrix<VNCS::Space2D::Real, Eigen::Dynamic, 1> p = dEig;

    Eigen::Matrix<VNCS::Space2D::Real, Eigen::Dynamic, 1> xC_eig;
    xC_eig.resize(2 * xC.size());
    for (const auto &[index, p] : ranges::views::enumerate(xC)) {
        xC_eig[2 * index + 0] = p[0];
        xC_eig[2 * index + 1] = p[1];
    }

    xC_eig = CtCInverse * (C.transpose() * dEig);
    dEig = p - C * xC_eig;

    for (const auto &[index, p] : ranges::views::enumerate(xC)) {
        p = VNCS::Space2D::Point{xC_eig[2 * index + 0], xC_eig[2 * index + 1]};
    }

    const auto dUntouched = d;
    for (const auto &[index, p] : ranges::views::enumerate(d)) {
        p = VNCS::Space2D::Point{dEig[2 * index + 0], dEig[2 * index + 1]};
    }

    // Create topologies
    // Coarse positions will be mapped with an IdentityMultiMapping
    // coarsePositions = [xN; xC]
    // Iterate of the triangles of the coarse topology and create the triangles

    std::vector<std::array<int, 3>> coarseTriangles;
    for (const auto &face : coarseMesh.faces()) {
        bool insertTriangleInTopology = true;
        std::array<int, 3> vertices;
        auto vertexIt = vertices.begin();
        for (auto pointId : coarseMesh.vertices_around_face(coarseMesh.halfedge(face))) {
            // Check if this vertex participates in the simulation
            // TODO concatenate this in a single loop
            auto xN_iterator = std::find(std::begin(xNIndices), std::end(xNIndices), pointId.idx());
            auto xC_iterator = std::find(std::begin(xCIndices), std::end(xCIndices), pointId.idx());

            if (xN_iterator != std::end(xNIndices) && xC_iterator != std::end(xCIndices)) {
                spdlog::get("VNCS")->info("This shouldnt happen");
            }

            if (xN_iterator != std::end(xNIndices) || xC_iterator != std::end(xCIndices)) {
                // This vertex participates in the simulation
                if (xN_iterator != std::end(xNIndices)) {
                    *vertexIt = std::distance(std::begin(xNIndices), xN_iterator);
                    vertexIt++;
                }
                if (xC_iterator != std::end(xCIndices)) {
                    *vertexIt = std::distance(std::begin(xCIndices), xC_iterator) + xNIndices.size();
                    vertexIt++;
                }
            } else {
                // The vertex is not simulated and therefore the whole tetrahedron is ignored
                insertTriangleInTopology = false;
            }
        }

        if (insertTriangleInTopology) {
            coarseTriangles.push_back(vertices);
        }
    }

    std::vector<std::array<int, 3>> fineTriangles;
    fineTriangles.clear();
    for (const auto &face : fineMesh.faces()) {
        bool insertTetraInTopology = true;
        std::array<int, 3> vertices;
        auto vertexIt = vertices.begin();
        for (auto pointId : fineMesh.vertices_around_face(fineMesh.halfedge(face))) {
            // Check if this vertex participates in the simulation
            // TODO concatenate this in a single loop
            auto d_iterator = std::find(std::begin(dIndices), std::end(dIndices), pointId.idx());
            auto pS_iterator = std::find(std::begin(pSIndices), std::end(pSIndices), pointId.idx());
            if (d_iterator != std::end(dIndices) || pS_iterator != std::end(pSIndices)) {
                // This vertex participates in the simulation
                if (d_iterator != std::end(dIndices)) {
                    *vertexIt = std::distance(std::begin(dIndices), d_iterator);
                    vertexIt++;
                }
                if (pS_iterator != std::end(pSIndices)) {
                    *vertexIt = std::distance(std::begin(pSIndices), pS_iterator) + dIndices.size();
                    vertexIt++;
                }
            } else {
                // The vertex doesn't
                insertTetraInTopology = false;
            }
        }

        if (insertTetraInTopology) {
            fineTriangles.push_back(vertices);
        }
    }

    {
        json dofJson = json::object();
        json x0 = json::object();
        auto xNjson = json::array();
        for (const auto &xN : xN) {
            auto v = json::array();
            v.push_back(xN[0]);
            v.push_back(xN[1]);
            xNjson.push_back(v);
        }
        x0["xN"] = xNjson;

        auto xCjson = json::array();
        for (const auto &xC : xC) {
            auto v = json::array();
            v.push_back(xC[0]);
            v.push_back(xC[1]);
            xCjson.push_back(v);
        }
        x0["xC"] = xCjson;

        auto dJson = json::array();
        for (const auto &d : d) {
            auto v = json::array();
            v.push_back(d[0]);
            v.push_back(d[1]);
            dJson.push_back(v);
        }
        x0["d"] = dJson;

        auto pSJson = json::array();
        for (const auto &pS : pS) {
            auto v = json::array();
            v.push_back(pS[0]);
            v.push_back(pS[1]);
            pSJson.push_back(v);
        }
        x0["pS"] = pSJson;

        dofJson["x0"] = x0;

        std::ofstream dof(m_dofFilePath);
        dof << std::setw(4) << dofJson << std::endl;
    }

    Mesh coarseSimulationMesh;
    {
        std::vector<VNCS::Space2D::Point> points2D = ranges::views::concat(xN, xC) |
                                                     ranges::views::transform([](const auto &p) {
                                                         return VNCS::Space2D::Point{p.x(), p.y()};
                                                     }) |
                                                     ranges::to_vector;
        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::polygon_soup_to_polygon_mesh(points2D, coarseTriangles, coarseSimulationMesh);
    }

    Mesh fineSimulationMesh;
    {
        std::vector<VNCS::Space2D::Point> points2D = ranges::views::concat(dUntouched, pS) |
                                                     ranges::views::transform([](const auto &p) {
                                                         return VNCS::Space2D::Point{p.x(), p.y()};
                                                     }) |
                                                     ranges::to_vector;
        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::polygon_soup_to_polygon_mesh(points2D, fineTriangles, fineSimulationMesh);
    }

    auto [weightPMap, isCreated] = fineSimulationMesh.add_property_map<Mesh::Vertex_index, VNCS::Space2D::Real>("w");

    // Populate the weights with values
    for (const auto vertexId : fineSimulationMesh.vertices())
        weightPMap[vertexId] = m_blendingField->blending(fineSimulationMesh.point(vertexId));

    // We need to slightly modify the weights.
    // In case a fine triangle is partially covered by coarse triangle, then all the weights must be 0
    const PMPTriangleIntersection coarseIntersectionFunc(coarseSimulationMesh);
    for (const auto &[faceId, fineFace] : fineSimulationMesh.faces() | ranges::views::enumerate) {
        const auto vIndices =
            fineSimulationMesh.vertices_around_face(fineSimulationMesh.halfedge(fineFace)) | ranges::to_vector;
        const auto v =
            vIndices |
            ranges::views::transform([fineSimulationMesh](const auto &vd) { return fineSimulationMesh.point(vd); }) |
            ranges::to_vector;
        const Triangle fineTriangle{v[0], v[1], v[2]};
        const std::array<std::optional<Mesh::Face_index>, 3> coarseTrianglesForFineTriangle = {
            coarseIntersectionFunc.query(v[0]),  //
            coarseIntersectionFunc.query(v[1]),
            coarseIntersectionFunc.query(v[2])};

        if (ranges::any_of(coarseTrianglesForFineTriangle, [](auto opt) { return opt.has_value(); }) &&
            ranges::any_of(coarseTrianglesForFineTriangle, [](auto opt) { return !opt.has_value(); })) {
            weightPMap[vIndices[0]] = 0.0;
            weightPMap[vIndices[1]] = 0.0;
            weightPMap[vIndices[2]] = 0.0;
        }
    }

    std::vector<VNCS::Real> fineBlending(weightPMap.begin(), weightPMap.end());

    m_exportFunction(coarseSimulationMesh, {}, "coarse.vtu");
    m_exportFunction(fineSimulationMesh, fineBlending, "fine.vtu");

    SamplingPointsExporter samplers;
    samplers.setCoarseSamplersFilePath(m_coarseSamplersFilePath);
    samplers.setFineSamplersFilePath(m_fineSamplersFilePath);
    samplers.setFineNodeSamplersFilePath(m_fineNodeSamplersFilePath);
    samplers(coarseSimulationMesh, fineSimulationMesh, m_gridSize);

    VisualSamplersExporter visualSamplers;
    visualSamplers.setSamplersFilePath(m_visualSamplersFilePath);
    visualSamplers(visualMesh, coarseSimulationMesh, fineSimulationMesh);
}

void VNCS::Generators::Gen222::Generator::setBlendingField(
    std::shared_ptr<VNCS::Generators::BlendingField<VNCS::Space2D::Real, VNCS::Space2D::Point>> blendingField)
{
    m_blendingField = blendingField;
}

std::filesystem::path VNCS::Generators::Gen222::Generator::coarseSamplersFilePath() const
{
    return m_coarseSamplersFilePath;
}

void VNCS::Generators::Gen222::Generator::setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath)
{
    m_coarseSamplersFilePath = coarseSamplersFilePath;
}

std::filesystem::path VNCS::Generators::Gen222::Generator::fineSamplersFilePath() const
{
    return m_fineSamplersFilePath;
}

void VNCS::Generators::Gen222::Generator::setFineSamplersFilePath(const std::filesystem::path &fineSamplersFilePath)
{
    m_fineSamplersFilePath = fineSamplersFilePath;
}

std::filesystem::path VNCS::Generators::Gen222::Generator::fineNodeSamplersFilePath() const
{
    return m_fineNodeSamplersFilePath;
}

void VNCS::Generators::Gen222::Generator::setFineNodeSamplersFilePath(
    const std::filesystem::path &fineNodeSamplersFilePath)
{
    m_fineNodeSamplersFilePath = fineNodeSamplersFilePath;
}

std::filesystem::path VNCS::Generators::Gen222::Generator::visualSamplersFilePath() const
{
    return m_visualSamplersFilePath;
}

void VNCS::Generators::Gen222::Generator::setVisualSamplersFilePath(const std::filesystem::path &visualSamplersFilePath)
{
    m_visualSamplersFilePath = visualSamplersFilePath;
}

std::filesystem::path VNCS::Generators::Gen222::Generator::clusterMatrixFilePath() const
{
    return m_clusterMatrixFilePath;
}

void VNCS::Generators::Gen222::Generator::setClusterMatrixFilePath(const std::filesystem::path &clusterMatrixFilePath)
{
    m_clusterMatrixFilePath = clusterMatrixFilePath;
}

std::filesystem::path VNCS::Generators::Gen222::Generator::dofFilePath() const
{
    return m_dofFilePath;
}

void VNCS::Generators::Gen222::Generator::setDofFilePath(const std::filesystem::path &dofFilePath)
{
    m_dofFilePath = dofFilePath;
}

void VNCS::Generators::Gen222::Generator::setGridSize(const std::size_t gridSize)
{
    m_gridSize = gridSize;
}

std::filesystem::path VNCS::Generators::Gen222::Generator::coarseMeshPath() const
{
    return m_coarseMeshPath;
}

void VNCS::Generators::Gen222::Generator::setCoarseMeshPath(const std::filesystem::path &coarseMeshPath)
{
    m_coarseMeshPath = coarseMeshPath;
}

std::filesystem::path VNCS::Generators::Gen222::Generator::fineMeshPath() const
{
    return m_fineMeshPath;
}

void VNCS::Generators::Gen222::Generator::setFineMeshPath(const std::filesystem::path &fineMeshPath)
{
    m_fineMeshPath = fineMeshPath;
}

std::filesystem::path VNCS::Generators::Gen222::Generator::visualMeshPath() const
{
    return m_visualMeshPath;
}

void VNCS::Generators::Gen222::Generator::setVisualMeshPath(const std::filesystem::path &visualMeshPath)
{
    m_visualMeshPath = visualMeshPath;
}

std::shared_ptr<VNCS::Generators::AdaptiveMeshCriteriaFunctor> VNCS::Generators::Gen222::Generator::coarseCriteria()
    const
{
    return m_coarseCriteria;
}

void VNCS::Generators::Gen222::Generator::setCoarseCriteria(
    const std::shared_ptr<AdaptiveMeshCriteriaFunctor> &coarseCriteria)
{
    m_coarseCriteria = coarseCriteria;
}

std::shared_ptr<VNCS::Generators::AdaptiveMeshCriteriaFunctor> VNCS::Generators::Gen222::Generator::fineCriteria() const
{
    return m_fineCriteria;
}

void VNCS::Generators::Gen222::Generator::setFineCriteria(
    const std::shared_ptr<AdaptiveMeshCriteriaFunctor> &fineCriteria)
{
    m_fineCriteria = fineCriteria;
}

std::shared_ptr<VNCS::Generators::AdaptiveMeshCriteriaFunctor> VNCS::Generators::Gen222::Generator::visualCriteria()
    const
{
    return m_visualCriteria;
}

void VNCS::Generators::Gen222::Generator::setVisualCriteria(
    const std::shared_ptr<AdaptiveMeshCriteriaFunctor> &visualCriteria)
{
    m_visualCriteria = visualCriteria;
}

void VNCS::Generators::Gen222::Generator::setExportFunction(
    std::function<void(const Space2D::Mesh &, const std::vector<VNCS::Real> &, const std::string &)> exportFunction)
{
    m_exportFunction = exportFunction;
}

bool VNCS::Generators::Gen222::Generator::allowBoundary() const
{
    return m_allowBoundary;
}

void VNCS::Generators::Gen222::Generator::setAllowBoundary(bool allowBoundary)
{
    m_allowBoundary = allowBoundary;
}