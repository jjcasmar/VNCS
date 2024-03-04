#include "Gen333.h"

#include "SamplingPointsExporter.h"
#include "VisualSamplersExporter.h"

#include <VNCS/Generators/TetraMeshGenerator.h>
#include <VNCS/SamplingPoints.h>
#include <VNCS/Generators/TetraMeshIntersection.h>

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <boost/iterator/counting_iterator.hpp>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/IO/print_wavefront.h>
#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <spdlog/fmt/fmt.h>

#include <range/v3/view/zip.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/unique.hpp>
#include <range/v3/view/group_by.hpp>
#include <range/v3/view/adjacent_remove_if.hpp>
#include <range/v3/action/unique.hpp>
#include <range/v3/action/sort.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/for_each.hpp>
#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/numeric/accumulate.hpp>

#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/IO/output_to_vtu.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <cmath>
#include <unsupported/Eigen/SparseExtra>

#include <spdlog/spdlog.h>
#include <spdlog/formatter.h>
#include <spdlog/pattern_formatter.h>
#include <spdlog/common.h>
#include <spdlog/sinks/base_sink.h>
#include <spdlog/fmt/chrono.h>
#include <spdlog/fmt/ostr.h>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace
{
namespace PMP = CGAL::Polygon_mesh_processing;
class My_point_property_map
{
    const std::vector<VNCS::Space3D::Point> &points;

public:
    typedef VNCS::Space3D::Point value_type;
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

typedef CGAL::Search_traits_3<VNCS::Space3D::K> TreeTraitsBase;
typedef CGAL::Search_traits_adapter<std::size_t, My_point_property_map, TreeTraitsBase> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;
typedef Tree::Splitter Splitter;
typedef Neighbor_search::Distance Distance;

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
    bool operator()(const std::array<int, 3> &a, const std::array<int, 3> &b) const
    {
        auto a_ = a;
        auto b_ = b;
        std::sort(std::begin(a_), std::end(a_));
        std::sort(std::begin(b_), std::end(b_));
        return a_[0] == b_[0] &&  //
               a_[1] == b_[1] &&  //
               a_[2] == b_[2];
    }
};

}  // namespace

namespace VNCS
{
namespace Generators
{
namespace Gen333
{
void Generator::operator()()
{
    const auto logger = spdlog::get("Gen333") ? spdlog::get("Gen333") : []() {
        std::cout << "There is no logger!\n";
        auto logger = std::make_shared<spdlog::logger>("Gen333");
        spdlog::register_logger(logger);
        return logger;
    }();

    logger->info("Executing Gen333");
    // Create the coarse and fine tetrahedral meshes
    VNCS::Space3D::Mesh coarseInputMesh;
    VNCS::Space3D::Mesh fineInputMesh;
    VNCS::Space3D::Mesh visualInputMesh;

    {
        logger->info("Loading coarse mesh from {}", m_coarseMeshPath);
        std::vector<VNCS::Space3D::Point> coarsePoints;
        std::vector<std::vector<std::size_t>> coarseFaces;
        std::ifstream coarseMeshIn(m_coarseMeshPath);
        CGAL::read_OBJ(coarseMeshIn, coarsePoints, coarseFaces);

        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::polygon_soup_to_polygon_mesh(coarsePoints, coarseFaces, coarseInputMesh);
    }

    {
        logger->info("Loading fine mesh from {}", m_fineMeshPath);
        std::vector<VNCS::Space3D::Point> coarsePoints;
        std::vector<std::vector<std::size_t>> coarseFaces;
        std::ifstream coarseMeshIn(m_fineMeshPath);
        CGAL::read_OBJ(coarseMeshIn, coarsePoints, coarseFaces);

        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::polygon_soup_to_polygon_mesh(coarsePoints, coarseFaces, fineInputMesh);
    }

    {
        logger->info("Loading visual mesh from {}", m_visualMeshPath);
        std::vector<VNCS::Space3D::Point> coarsePoints;
        std::vector<std::vector<std::size_t>> coarseFaces;
        std::ifstream coarseMeshIn(m_visualMeshPath);
        CGAL::read_OBJ(coarseMeshIn, coarsePoints, coarseFaces);

        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::polygon_soup_to_polygon_mesh(coarsePoints, coarseFaces, visualInputMesh);
    }

    logger->info("Generating coarse tetrahedra mesh");
    const auto coarseMesh = VNCS::Generators::makeC3T3FromMesh(coarseInputMesh, m_coarseCriteria);
    const auto [coarseMeshPoints, coarseMeshTetra] = VNCS::Generators::c3t3ToTetrahedra(coarseMesh);
    const auto coarseMeshBlendingValues =
        coarseMeshPoints | ranges::views::transform([this](const auto &p) { return m_blendingField->blending(p); }) |
        ranges::to_vector;
    std::unordered_map<std::string, std::vector<VNCS::Space3D::Real>> attributes;
    attributes["b"] = coarseMeshBlendingValues;
    m_exportFunction({coarseMeshPoints, coarseMeshTetra}, "coarseTest.vtu", attributes);

    logger->info("Generating fine tetrahedra mesh");
    const auto fineMesh = VNCS::Generators::makeC3T3FromMesh(fineInputMesh, m_fineCriteria);
    const auto [fineMeshPoints, fineMeshTetra] = VNCS::Generators::c3t3ToTetrahedra(fineMesh);
    const auto fineMeshBlendingValues =
        fineMeshPoints | ranges::views::transform([this](const auto &p) { return m_blendingField->blending(p); }) |
        ranges::to_vector;

    attributes["b"] = fineMeshBlendingValues;
    m_exportFunction({fineMeshPoints, fineMeshTetra}, "fineTest.vtu", attributes);

    logger->info("Generating visual surface mesh");
    auto visualMesh = visualInputMesh;
    {
        VNCS::Space3D::Mesh::Property_map<VNCS::Space3D::Mesh::Edge_index, bool> is_constrained =
            visualMesh.add_property_map<VNCS::Space3D::Mesh::Edge_index, bool>("e:is_constrained", false).first;

        // detect sharp features
        for (const auto edge : visualMesh.edges()) {
            auto hd = visualMesh.halfedge(edge);
            if (!visualMesh.is_border(hd)) {
                double angle = CGAL::Mesh_3::dihedral_angle(
                    visualMesh.point(source(hd, visualMesh)),
                    visualMesh.point(target(hd, visualMesh)),
                    visualMesh.point(target(next(hd, visualMesh), visualMesh)),
                    visualMesh.point(target(next(opposite(hd, visualMesh), visualMesh), visualMesh)));
                if (CGAL::abs(angle) < m_visualCriteria.sharpAngle)
                    is_constrained[edge] = true;
            }
        }

        // remesh
        PMP::isotropic_remeshing(
            faces(visualMesh),
            m_visualCriteria.target_edge_length,
            visualMesh,
            PMP::parameters::number_of_iterations(m_visualCriteria.iterations).edge_is_constrained_map(is_constrained));

        visualMesh.collect_garbage();
    }

    {
        CGAL::Polyhedron_3<VNCS::Space3D::K> visualOutputMesh;

        std::vector<VNCS::Space3D::K::Point_3> points =
            visualMesh.vertices() |
            ranges::views::transform([visualMesh](const auto &v) { return visualMesh.point(v); }) | ranges::to_vector;
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
        logger->info("Visual output mesh has {} points and {} triangles",
                     visualOutputMesh.size_of_vertices(),
                     visualOutputMesh.size_of_facets());

        std::ofstream outVisual("visual.obj");
        CGAL::print_polyhedron_wavefront(outVisual, visualOutputMesh);
    }

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

    std::vector<int> boundaryPointsIndices;
    {
        std::vector<std::array<int, 3>> tetraTriangles;
        for (const std::array<int, 4> &tetra : coarseMeshTetra) {
            tetraTriangles.push_back({tetra[0], tetra[1], tetra[2]});
            tetraTriangles.push_back({tetra[1], tetra[2], tetra[3]});
            tetraTriangles.push_back({tetra[0], tetra[2], tetra[3]});
            tetraTriangles.push_back({tetra[0], tetra[1], tetra[3]});
        }

        // Order the triangles and remove repeated ones, no matter the order of the indices: (a,b,c) is the same
        // triangle as (b,a,c)
        std::sort(std::begin(tetraTriangles), std::end(tetraTriangles), STriangleLessComparator{});
        for (auto it = std::adjacent_find(tetraTriangles.begin(), tetraTriangles.end(), STriangleEqualComparator{});  //
             it != tetraTriangles.end();                                                                              //
             it = std::adjacent_find(tetraTriangles.begin(), tetraTriangles.end(), STriangleEqualComparator{}))       //
            tetraTriangles.erase(it, it + 2);

        // The remaining triangles are the boundary
        for (const std::array<int, 3> &boundaryTriangle : tetraTriangles) {
            boundaryPointsIndices.push_back(boundaryTriangle[0]);
            boundaryPointsIndices.push_back(boundaryTriangle[1]);
            boundaryPointsIndices.push_back(boundaryTriangle[2]);
        }

        // Some points are repeated, keep one single copy of boundary vertices
        std::sort(std::begin(boundaryPointsIndices), std::end(boundaryPointsIndices));
        for (auto it = std::adjacent_find(boundaryPointsIndices.begin(), boundaryPointsIndices.end());  //
             it != boundaryPointsIndices.end();                                                         //
             it = std::adjacent_find(boundaryPointsIndices.begin(), boundaryPointsIndices.end()))       //
            boundaryPointsIndices.erase(it + 1, it + 2);
    }

    std::vector<unsigned int> xNIndices;
    std::vector<unsigned int> xCIndices;

    {
        for (const auto &tetra : coarseMeshTetra) {
            std::array<VNCS::Space3D::Real, 4> tetraBlendingValues = {coarseMeshBlendingValues[tetra[0]],
                                                                      coarseMeshBlendingValues[tetra[1]],
                                                                      coarseMeshBlendingValues[tetra[2]],
                                                                      coarseMeshBlendingValues[tetra[3]]};

            auto isCoarseDomain = [&tetraBlendingValues]() {
                return std::all_of(std::begin(tetraBlendingValues),
                                   std::end(tetraBlendingValues),
                                   [](auto blendingValue) { return blendingValue == 1.0; });
            };
            auto isHybridDomain = [&tetraBlendingValues]() {
                return std::any_of(std::begin(tetraBlendingValues),
                                   std::end(tetraBlendingValues),
                                   [](auto blendingValue) { return blendingValue > 0.0 && blendingValue < 1.0; });
            };

            // Tetrahedron appears in topology if its in the coarse region or in the hybrid region
            // If tetrahedron is in the coarse region, all its vertices belong to the xN set
            // if tetrahedron is in the hybrid region, we must identify three cases:
            //    - vertex is in the mesh boundary -> vertex goes into the xN set
            //    - vertex has a blending factor of 0 or 1 -> vertex goes into the xN set
            //    - otherwise, it goes nito the xC set
            if (isCoarseDomain()) {
                xNIndices.push_back(tetra[0]);
                xNIndices.push_back(tetra[1]);
                xNIndices.push_back(tetra[2]);
                xNIndices.push_back(tetra[3]);
            }

            if (isHybridDomain()) {
                for (int i = 0; i < 4; ++i) {
                    const auto &pointId = tetra[i];
                    if (std::find(std::begin(boundaryPointsIndices), std::end(boundaryPointsIndices), pointId) !=
                        std::end(boundaryPointsIndices)) {
                        xNIndices.push_back(pointId);
                    } else {
                        const auto &blendingValue = tetraBlendingValues[i];
                        if (blendingValue >= 1.0 || blendingValue <= 0.0)
                            xNIndices.push_back(pointId);
                        else
                            xCIndices.push_back(pointId);
                    }
                }
            }
        }

        // Some vertices may be repeated, we need to remove them!
        ranges::actions::sort(xCIndices);
        ranges::actions::sort(xNIndices);
        xCIndices = ranges::views::adjacent_remove_if(xCIndices, [](const auto &a, const auto &b) { return a == b; }) |
                    ranges::to_vector;
        xNIndices = ranges::views::adjacent_remove_if(xNIndices, [](const auto &a, const auto &b) { return a == b; }) |
                    ranges::to_vector;
    }
    logger->info("xN has {} points", xNIndices.size());
    logger->info("xC has {} points", xCIndices.size());

    // Later we will modify xC, so cant be const
    std::vector<VNCS::Space3D::Point> xC =
        xCIndices | ranges::views::transform([&coarseMeshPoints](const auto vId) { return coarseMeshPoints[vId]; }) |
        ranges::to_vector;

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

    for (const auto &tetra : fineMeshTetra) {
        std::array<VNCS::Space3D::Real, 4> tetraBlendingValues = {fineMeshBlendingValues[tetra[0]],  //
                                                                  fineMeshBlendingValues[tetra[1]],
                                                                  fineMeshBlendingValues[tetra[2]],
                                                                  fineMeshBlendingValues[tetra[3]]};

        if (ranges::any_of(tetraBlendingValues, [](const auto v) { return v != 1.0; })) {
            for (const auto &[index, blending] : ranges::views::zip(tetra, tetraBlendingValues)) {
                if (blending != 0.0)
                    dIndices.push_back(index);
                else
                    pSIndices.push_back(index);
            }
        }

        if (ranges::all_of(tetraBlendingValues, [](const auto v) { return v == 0.0; })) {
            pSIndices.push_back(tetra[0]);
            pSIndices.push_back(tetra[1]);
            pSIndices.push_back(tetra[2]);
            pSIndices.push_back(tetra[3]);
        }
    }

    dIndices = ranges::actions::sort(dIndices) |
               ranges::views::adjacent_remove_if([](const auto &a, const auto &b) { return a == b; }) |
               ranges::to_vector;

    pSIndices = ranges::actions::sort(pSIndices) |
                ranges::views::adjacent_remove_if([](const auto &a, const auto &b) { return a == b; }) |
                ranges::to_vector;
    logger->info("d has {} points", dIndices.size());
    logger->info("pS has {} points", pSIndices.size());

    std::vector<VNCS::Space3D::Point> d;
    d.reserve(dIndices.size());
    for (auto dIndex : dIndices) {
        const auto p = fineMeshPoints[dIndex];
        d.push_back(p);
    }

    std::vector<VNCS::Space3D::Point> pS;
    pS.reserve(pSIndices.size());
    for (auto pSIndex : pSIndices) {
        const auto p = fineMeshPoints[pSIndex];
        pS.push_back(p);
    }

    // We have to find the closest xC vertex to each d vertex
    // If any xC vertex is not selected, we will crash later!
    logger->info("Computing closest cluster points");
    My_point_property_map ppmap(xC);
    Tree searchTree(boost::counting_iterator<std::size_t>(0),
                    boost::counting_iterator<std::size_t>(xC.size()),
                    Splitter(),
                    TreeTraits(ppmap));
    Distance trDistance(ppmap);

    std::vector<unsigned int> assignedXCvertex(dIndices.size());
    for (int i = 0; i < dIndices.size(); ++i) {
        const auto &dPoint = fineMeshPoints[dIndices[i]];
        Neighbor_search search(searchTree, dPoint, 1, 0.0001, true, trDistance);
        for (const auto &searchResult : search) {
            assignedXCvertex[i] = xCIndices[searchResult.first];
        }
    }

    std::vector<Eigen::Triplet<VNCS::Space3D::Real>> triplets;
    triplets.reserve(3 * dIndices.size());
    for (int i = 0; i < dIndices.size(); ++i) {
        auto it = std::find(std::begin(xCIndices), std::end(xCIndices), assignedXCvertex[i]);
        auto index = std::distance(std::begin(xCIndices), it);
        triplets.push_back(Eigen::Triplet<VNCS::Space2D::Real>(3 * i + 0, 3 * index + 0, 1.0));
        triplets.push_back(Eigen::Triplet<VNCS::Space2D::Real>(3 * i + 1, 3 * index + 1, 1.0));
        triplets.push_back(Eigen::Triplet<VNCS::Space2D::Real>(3 * i + 2, 3 * index + 2, 1.0));
    }

    // We have to find the closest xC vertex to each d vertex
    // If any xC vertex is not selected, we will crash later!

    // Create xN degrees of freedom
    std::vector<VNCS::Space3D::Point> xN;
    xN.reserve(xNIndices.size());
    for (auto xNVertexId : xNIndices) {
        const auto p = coarseMeshPoints[xNVertexId];
        xN.push_back(p);
    }

    Eigen::SparseMatrix<VNCS::Space3D::Real> C;
    C.resize(3 * dIndices.size(),  //
             3 * xCIndices.size());
    C.setFromTriplets(std::begin(triplets), std::end(triplets), [](const auto &a, const auto &b) { return a; });

    Eigen::saveMarket(C, m_clusterMatrixFilePath);

    Eigen::SparseMatrix<VNCS::Space3D::Real> CtC = (C.transpose() * C);
    Eigen::DiagonalMatrix<VNCS::Space3D::Real, Eigen::Dynamic> CtCInverse;
    CtCInverse.resize(CtC.rows());
    for (int i = 0; i < CtC.rows(); ++i) {
        if (CtC.coeff(i, i) == VNCS::Space2D::Real{0.0})
            logger->error("d vertex {} doesn't have an associated xC node", i / 3);
        CtCInverse.diagonal()[i] = 1.0 / CtC.coeff(i, i);
    }

    Eigen::Matrix<VNCS::Space3D::Real, Eigen::Dynamic, 1> dEig;
    dEig.resize(3 * d.size());
    for (const auto &[index, p] : ranges::views::enumerate(d)) {
        dEig[3 * index + 0] = p[0];
        dEig[3 * index + 1] = p[1];
        dEig[3 * index + 2] = p[2];
    }

    Eigen::Matrix<VNCS::Space3D::Real, Eigen::Dynamic, 1> p = dEig;

    Eigen::Matrix<VNCS::Space3D::Real, Eigen::Dynamic, 1> xC_eig;
    xC_eig.resize(3 * xC.size());
    for (const auto &[index, p] : ranges::views::enumerate(xC)) {
        xC_eig[3 * index + 0] = p[0];
        xC_eig[3 * index + 1] = p[1];
        xC_eig[3 * index + 2] = p[2];
    }

    xC_eig = CtCInverse * (C.transpose() * dEig);
    dEig = p - C * xC_eig;

    for (const auto &[index, p] : ranges::views::enumerate(xC)) {
        p = VNCS::Space3D::Point{xC_eig[3 * index + 0], xC_eig[3 * index + 1], xC_eig[3 * index + 2]};
    }

    const auto dUntouched = d;
    for (const auto &[index, p] : ranges::views::enumerate(d)) {
        p = VNCS::Space3D::Point{dEig[3 * index + 0], dEig[3 * index + 1], dEig[3 * index + 2]};
    }

    // Create topologies
    // Coarse positions will be mapped with an IdentityMultiMapping
    // coarsePositions = [xN; xC]
    // Iterate of the triangles of the coarse topology and create the triangles

    logger->info("Creating simulation topologies");

    std::vector<std::array<int, 4>> coarseTetrahedra;
    for (const auto &tetra : coarseMeshTetra) {
        bool insertTetraInTopology = true;
        std::array<int, 4> vertices;
        auto vertexIt = vertices.begin();
        for (auto pointId : tetra) {
            // Check if this vertex participates in the simulation
            // TODO concatenate this in a single loop
            auto xN_iterator = std::find(std::begin(xNIndices), std::end(xNIndices), pointId);
            auto xC_iterator = std::find(std::begin(xCIndices), std::end(xCIndices), pointId);

            if (xN_iterator != std::end(xNIndices) && xC_iterator != std::end(xCIndices)) {
                throw std::runtime_error("");
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
                insertTetraInTopology = false;
            }
        }

        if (insertTetraInTopology) {
            coarseTetrahedra.push_back(vertices);
        }
    }

    std::vector<std::array<int, 4>> fineTetrahedra;
    for (const auto &tetra : fineMeshTetra) {
        bool insertTetraInTopology = true;
        std::array<int, 4> vertices;
        auto vertexIt = vertices.begin();
        for (auto pointId : tetra) {
            // Check if this vertex participates in the simulation
            // TODO concatenate this in a single loop
            auto d_iterator = std::find(std::begin(dIndices), std::end(dIndices), pointId);
            auto pS_iterator = std::find(std::begin(pSIndices), std::end(pSIndices), pointId);
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
            fineTetrahedra.push_back(vertices);
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
            v.push_back(xN[2]);
            xNjson.push_back(v);
        }
        x0["xN"] = xNjson;

        auto xCjson = json::array();
        for (const auto &xC : xC) {
            auto v = json::array();
            v.push_back(xC[0]);
            v.push_back(xC[1]);
            v.push_back(xC[2]);
            xCjson.push_back(v);
        }
        x0["xC"] = xCjson;

        auto dJson = json::array();
        for (const auto &d : d) {
            auto v = json::array();
            v.push_back(d[0]);
            v.push_back(d[1]);
            v.push_back(d[2]);
            dJson.push_back(v);
        }
        x0["d"] = dJson;

        auto pSJson = json::array();
        for (const auto &pS : pS) {
            auto v = json::array();
            v.push_back(pS[0]);
            v.push_back(pS[1]);
            v.push_back(pS[2]);
            pSJson.push_back(v);
        }
        x0["pS"] = pSJson;

        dofJson["x0"] = x0;

        std::ofstream dof(m_dofFilePath);
        dof << std::setw(4) << dofJson << std::endl;
    }

    VNCS::Space3D::TetraMesh coarseSimulationMesh = {ranges::views::concat(xN, xC) | ranges::to_vector,
                                                     coarseTetrahedra};
    if (coarseSimulationMesh.tetras.size())
        m_exportFunction(coarseSimulationMesh, "coarse.vtu", {});

    VNCS::Space3D::TetraMesh fineSimulationMesh = {ranges::views::concat(dUntouched, pS) | ranges::to_vector,
                                                   fineTetrahedra};
    if (fineSimulationMesh.tetras.size())
        m_exportFunction(fineSimulationMesh, "fine.vtu", {});

    logger->info("Coarse tetrahedra mesh has {} points and {} tetrahedra",
                 coarseSimulationMesh.points.size(),
                 coarseSimulationMesh.tetras.size());
    logger->info("Fine tetrahedra mesh has {} points and {} tetrahedra",
                 fineSimulationMesh.points.size(),
                 fineSimulationMesh.tetras.size());

    SamplingPointsExporter samplers;
    samplers.setCoarseSamplersFilePath(m_coarseSamplersFilePath);
    samplers.setFineSamplersFilePath(m_fineSamplersFilePath);
    samplers(coarseSimulationMesh, fineSimulationMesh);

    //VisualSamplersExporter visualSamplers;
    //visualSamplers.setSamplersFilePath(m_visualSamplersFilePath);
    //visualSamplers(visualMesh, coarseSimulationMesh, fineSimulationMesh, weights, m_blendingField);
}

void Generator::setBlendingField(std::shared_ptr<BlendingField<Space3D::Real, Space3D::Point>> blendingField)
{
    m_blendingField = blendingField;
}

std::filesystem::path Generator::coarseMeshPath() const
{
    return m_coarseMeshPath;
}
void Generator::setCoarseMeshPath(const std::filesystem::path &coarseMeshPath)
{
    m_coarseMeshPath = coarseMeshPath;
}

std::filesystem::path Generator::fineMeshPath() const
{
    return m_fineMeshPath;
}
void Generator::setFineMeshPath(const std::filesystem::path &fineMeshPath)
{
    m_fineMeshPath = fineMeshPath;
}

std::filesystem::path Generator::visualMeshPath() const
{
    return m_visualMeshPath;
}
void Generator::setVisualMeshPath(const std::filesystem::path &visualMeshPath)
{
    m_visualMeshPath = visualMeshPath;
}

VNCS::Generators::C3t3Criteria Generator::coarseCriteria() const
{
    return m_coarseCriteria;
}

void Generator::setCoarseCriteria(const VNCS::Generators::C3t3Criteria &coarseCriteria)
{
    m_coarseCriteria = coarseCriteria;
}

VNCS::Generators::C3t3Criteria Generator::fineCriteria() const
{
    return m_fineCriteria;
}

void Generator::setFineCriteria(const VNCS::Generators::C3t3Criteria &fineCriteria)
{
    m_fineCriteria = fineCriteria;
}

VNCS::Generators::RemeshCriteria Generator::visualCriteria() const
{
    return m_visualCriteria;
}

void Generator::setVisualCriteria(const VNCS::Generators::RemeshCriteria &visualCriteria)
{
    m_visualCriteria = visualCriteria;
}

void Generator::setExportFunction(
    std::function<void(const Space3D::TetraMesh &,
                       const std::string &,
                       const std::unordered_map<std::string, std::vector<VNCS::Space3D::Real>> &)> exportFunction)
{
    m_exportFunction = exportFunction;
}

std::filesystem::path Generator::coarseSamplersFilePath() const
{
    return m_coarseSamplersFilePath;
}

void Generator::setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath)
{
    m_coarseSamplersFilePath = coarseSamplersFilePath;
}

std::filesystem::path Generator::fineSamplersFilePath() const
{
    return m_fineSamplersFilePath;
}

void Generator::setFineSamplersFilePath(const std::filesystem::path &fineSamplersFilePath)
{
    m_fineSamplersFilePath = fineSamplersFilePath;
}

std::filesystem::path Generator::visualSamplersFilePath() const
{
    return m_visualSamplersFilePath;
}

void Generator::setVisualSamplersFilePath(const std::filesystem::path &visualSamplersFilePath)
{
    m_visualSamplersFilePath = visualSamplersFilePath;
}

std::filesystem::path Generator::clusterMatrixFilePath() const
{
    return m_clusterMatrixFilePath;
}

void Generator::setClusterMatrixFilePath(const std::filesystem::path &clusterMatrixFilePath)
{
    m_clusterMatrixFilePath = clusterMatrixFilePath;
}

std::filesystem::path Generator::dofFilePath() const
{
    return m_dofFilePath;
}

void Generator::setDofFilePath(const std::filesystem::path &dofFilePath)
{
    m_dofFilePath = dofFilePath;
}

}  // namespace Gen333
}  // namespace Generators
}  // namespace VNCS
