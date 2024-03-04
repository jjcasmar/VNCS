#include "Gen223.h"

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <VNCS/SamplingPoints.h>
#include <VNCS/Generators/IntegrationScheme.h>
#include <VNCS/Generators/Gen223/SamplingPointsExporter.h>

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
#include <range/v3/view/group_by.hpp>
#include <range/v3/action/adjacent_remove_if.hpp>
#include <range/v3/view/adjacent_remove_if.hpp>
#include <range/v3/action/sort.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/for_each.hpp>
#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/algorithm/contains.hpp>

#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/IO/output_to_vtu.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <cmath>
#include <unsupported/Eigen/SparseExtra>

#include <nlohmann/json.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/formatter.h>
#include <spdlog/pattern_formatter.h>
#include <spdlog/common.h>
#include <spdlog/sinks/base_sink.h>
#include <spdlog/fmt/chrono.h>
#include <spdlog/fmt/ostr.h>

#include <chrono>

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
namespace Gen223
{
void Generator::operator()()
{
    const auto logger = spdlog::get("Gen223") ? spdlog::get("Gen223") : []() {
        std::cout << "There is no logger!\n";
        auto logger = std::make_shared<spdlog::logger>("Gen223");
        spdlog::register_logger(logger);
        return logger;
    }();

    logger->info("Executing Gen223");
    // Create the coarse and fine tetrahedral meshes
    VNCS::Space3D::Mesh coarseInputMesh;
    VNCS::Space3D::Mesh fineInputMesh;

    {
        logger->info("Loading coarse mesh from {}", m_coarseMeshPath);
        std::vector<VNCS::Space3D::Point> coarsePoints;
        std::vector<std::vector<std::size_t>> coarseFaces;
        std::ifstream coarseMeshIn(m_coarseMeshPath);
        CGAL::read_OBJ(coarseMeshIn, coarsePoints, coarseFaces);

        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::orient_polygon_soup(coarsePoints, coarseFaces);
        PMP::polygon_soup_to_polygon_mesh(coarsePoints, coarseFaces, coarseInputMesh);
    }

    {
        logger->info("Loading fine mesh from {}", m_fineMeshPath);
        std::vector<VNCS::Space3D::Point> coarsePoints;
        std::vector<std::vector<std::size_t>> coarseFaces;
        std::ifstream coarseMeshIn(m_fineMeshPath);
        CGAL::read_OBJ(coarseMeshIn, coarsePoints, coarseFaces);

        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::orient_polygon_soup(coarsePoints, coarseFaces);
        PMP::polygon_soup_to_polygon_mesh(coarsePoints, coarseFaces, fineInputMesh);
    }

    auto remesh = [](const VNCS::Space3D::Mesh &inputMesh, const auto &criteria) {
        auto mesh = inputMesh;
        {
            VNCS::Space3D::Mesh::Property_map<VNCS::Space3D::Mesh::Edge_index, bool> is_constrained =
                mesh.add_property_map<VNCS::Space3D::Mesh::Edge_index, bool>("e:is_constrained", false).first;

            VNCS::Space3D::Mesh::Property_map<VNCS::Space3D::Mesh::Vertex_index, bool> vertex_is_constrained =
                mesh.add_property_map<VNCS::Space3D::Mesh::Vertex_index, bool>("v:is_constrained", false).first;

            // detect sharp features
            for (const auto edge : mesh.edges()) {
                auto hd = mesh.halfedge(edge);
                if (!mesh.is_border(hd)) {
                    double angle =
                        CGAL::Mesh_3::dihedral_angle(mesh.point(source(hd, mesh)),
                                                     mesh.point(target(hd, mesh)),
                                                     mesh.point(target(next(hd, mesh), mesh)),
                                                     mesh.point(target(next(opposite(hd, mesh), mesh), mesh)));
                    if (CGAL::abs(angle) < criteria.sharpAngle)
                        is_constrained[edge] = true;
                }
            }

            // remesh
            PMP::isotropic_remeshing(faces(mesh),
                                     criteria.target_edge_length,
                                     mesh,
                                     PMP::parameters::number_of_iterations(criteria.iterations)
                                         .edge_is_constrained_map(is_constrained)
                                         .vertex_is_constrained_map(vertex_is_constrained));

            mesh.collect_garbage();
            return mesh;
        }
    };

    auto smooth = [](const VNCS::Space3D::Mesh &inputMesh) {
        auto mesh = inputMesh;
        typedef boost::property_map<VNCS::Space3D::Mesh, CGAL::edge_is_feature_t>::type EIFMap;
        EIFMap eif = get(CGAL::edge_is_feature, mesh);
        PMP::detect_sharp_edges(mesh, 60, eif);

        // Smooth with both angle and area criteria + Delaunay flips
        PMP::smooth_mesh(mesh,
                         PMP::parameters::number_of_iterations(20)
                             .use_safety_constraints(false)
                             .do_project(false)
                             .use_area_smoothing(false));

        return mesh;
    };

    auto coarseRemeshMesh = m_coarseCriteria ? remesh(coarseInputMesh, m_coarseCriteria.value()) : coarseInputMesh;
    auto fineRemeshMesh = m_fineCriteria ? remesh(fineInputMesh, m_fineCriteria.value()) : fineInputMesh;

    auto coarseMesh = coarseRemeshMesh;
    auto fineMesh = smooth(fineRemeshMesh);

    // Get the indices of the possible cluster points. This are the non boundary vertices
    const auto possibleClusterPointsIndices = coarseMesh.vertices() | ranges::to_vector;

    const auto possibleClusterPoints =
        possibleClusterPointsIndices |
        ranges::views::transform([&coarseMesh](const auto vId) { return coarseMesh.point(vId); }) | ranges::to_vector;

    ::My_point_property_map ppmap(possibleClusterPoints);
    Tree searchTree(boost::counting_iterator<std::size_t>(0),
                    boost::counting_iterator<std::size_t>(possibleClusterPoints.size()),
                    Splitter(),
                    TreeTraits(ppmap));
    Distance trDistance(ppmap);

    // xC points are only those that are required because they are the closest to a fine point
    auto xCIndices =
        fineMesh.vertices() |
        ranges::views::transform([&fineMesh](const auto vertex) { return fineMesh.point(vertex); }) |
        ranges::views::transform([&possibleClusterPointsIndices, &searchTree, &trDistance](const auto &point) {
            Neighbor_search search(searchTree, point, 1, 0.0001, true, trDistance);
            return possibleClusterPointsIndices[search.begin()->first];
        }) |
        ranges::to_vector | ranges::actions::sort |
        ranges::actions::adjacent_remove_if([](const auto a, const auto b) { return a == b; });

    auto xC = xCIndices | ranges::views::transform([&coarseMesh](const auto vId) { return coarseMesh.point(vId); }) |
              ranges::to_vector;

    My_point_property_map ppmap2(xC);
    Tree searchTree2(boost::counting_iterator<std::size_t>(0),
                     boost::counting_iterator<std::size_t>(xC.size()),
                     Splitter(),
                     TreeTraits(ppmap2));
    Distance trDistance2(ppmap2);

    std::vector<decltype(xCIndices)::value_type> assignedXCvertices(fineMesh.num_vertices());
    for (const auto [i, vertexHandle] : fineMesh.vertices() | ranges::views::enumerate) {
        const auto &dPoint = fineMesh.point(vertexHandle);
        Neighbor_search search2(searchTree2, dPoint, 1, 0.0001, true, trDistance2);
        for (const auto searchResult : search2) {
            assignedXCvertices[i] = xCIndices[searchResult.first];
        }
    }

    const auto xNIndices =
        coarseMesh.vertices() |
        ranges::views::filter([&xCIndices](const auto vId) { return !ranges::contains(xCIndices, vId); }) |
        ranges::to_vector;

    std::vector<Eigen::Triplet<VNCS::Real>> triplets;
    triplets.reserve(3 * fineMesh.num_vertices());
    for (const auto [i, vertexHandle] : fineMesh.vertices() | ranges::views::enumerate) {
        const auto point = fineMesh.point(vertexHandle);
        auto it = std::find(std::begin(xCIndices), std::end(xCIndices), assignedXCvertices[i]);
        auto index = std::distance(std::begin(xCIndices), it);
        triplets.push_back(Eigen::Triplet<VNCS::Space2D::Real>(3 * i + 0, 3 * index + 0, 1.0));
        triplets.push_back(Eigen::Triplet<VNCS::Space2D::Real>(3 * i + 1, 3 * index + 1, 1.0));
        triplets.push_back(Eigen::Triplet<VNCS::Space2D::Real>(3 * i + 2, 3 * index + 2, 1.0));
    }

    // Create xN degrees of freedom
    std::vector<VNCS::Space3D::Point> xN;
    xN.reserve(xNIndices.size());
    for (auto xNVertexId : xNIndices) {
        const auto p = coarseMesh.point(xNVertexId);
        xN.push_back({p.x(), p.y(), p.z()});
    }

    Eigen::SparseMatrix<VNCS::Real> C;
    C.resize(3 * fineMesh.num_vertices(),  //
             3 * xCIndices.size());
    C.setFromTriplets(std::begin(triplets), std::end(triplets), [](const auto &a, const auto &b) {
        spdlog::get("VNCS")->warn("Repeated triplets shouldnt happen");
        return a;
    });

    Eigen::saveMarket(C, m_clusterMatrixFilePath);

    Eigen::SparseMatrix<VNCS::Real> CtC = (C.transpose() * C);
    Eigen::DiagonalMatrix<VNCS::Real, Eigen::Dynamic> CtCInverse;
    CtCInverse.resize(CtC.rows());
    for (int i = 0; i < CtC.rows(); ++i) {
        if (CtC.coeff(i, i) == VNCS::Real{0.0})
            spdlog::get("VNCS")->error("Coarse node {} doesnt have associated d nodes", xN.size() + i / 3);
        CtCInverse.diagonal()[i] = 1.0 / CtC.coeff(i, i);
    }

    Eigen::Matrix<VNCS::Real, Eigen::Dynamic, 1> dEig;
    dEig.resize(3 * boost::num_vertices(fineMesh));
    for (const auto &[index, vertexHandle] : ranges::views::enumerate(fineMesh.vertices())) {
        const auto &p = fineMesh.point(vertexHandle);
        dEig[3 * index + 0] = p[0];
        dEig[3 * index + 1] = p[1];
        dEig[3 * index + 2] = p[2];
    }

    Eigen::Matrix<VNCS::Real, Eigen::Dynamic, 1> p = dEig;

    Eigen::Matrix<VNCS::Real, Eigen::Dynamic, 1> xC_eig;
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

    std::vector<VNCS::Space3D::Point> dValues;
    for (const auto &[index, vertexHandle] : ranges::views::enumerate(fineMesh.vertices()))
        dValues.push_back(VNCS::Space3D::Point{dEig[3 * index + 0], dEig[3 * index + 1], dEig[3 * index + 2]});

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
        for (const auto &p : dValues) {
            auto v = json::array();
            v.push_back(p[0]);
            v.push_back(p[1]);
            v.push_back(p[2]);
            dJson.push_back(v);
        }
        x0["d"] = dJson;

        dofJson["x0"] = x0;

        std::ofstream dof(m_dofFilePath);
        dof << std::setw(4) << dofJson << std::endl;
    }

    // Create the simulation tetra mesh
    CGAL::Polyhedron_3<VNCS::Space3D::K> coarseMeshExport;
    std::vector<std::array<int, 3>> coarseTriangles;
    for (const auto &face : coarseMesh.faces()) {
        bool insertTetraInTopology = true;
        std::array<int, 3> vertices;
        auto vertexIt = vertices.begin();
        for (const auto vId : coarseMesh.vertices_around_face(coarseMesh.halfedge(face))) {
            // Check if this vertex participates in the simulation
            // TODO concatenate this in a single loop
            auto xN_iterator = std::find(std::begin(xNIndices), std::end(xNIndices), vId);
            auto xC_iterator = std::find(std::begin(xCIndices), std::end(xCIndices), vId);

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
            coarseTriangles.push_back(vertices);
        }
    }

    VNCS::Space3D::Mesh coarseSimulationMesh;
    std::vector<VNCS::Space3D::Mesh::Vertex_index> vertexIndices;
    for (const auto xNPoint : xN) {
        vertexIndices.push_back(coarseSimulationMesh.add_vertex(xNPoint));
    }

    for (const auto xCPoint : xC) {
        vertexIndices.push_back(coarseSimulationMesh.add_vertex(xCPoint));
    }

    for (const auto triangle : coarseTriangles) {
        const auto index0 = vertexIndices[triangle[0]];
        const auto index1 = vertexIndices[triangle[1]];
        const auto index2 = vertexIndices[triangle[2]];

        coarseSimulationMesh.add_face(index0, index1, index2);
    }

    {
        CGAL::Polyhedron_3<VNCS::Space3D::K> coarseMeshExport;

        std::vector<VNCS::Space2D::K::Point_3> points =
            coarseSimulationMesh.vertices() | ranges::views::transform([&coarseSimulationMesh](const auto vId) {
                const auto &p = coarseSimulationMesh.point(vId);
                return VNCS::Space2D::K::Point_3{p.x(), p.y(), p.z()};
            }) |
            ranges::to_vector;

        const auto fineTriangles =
            coarseSimulationMesh.faces() | ranges::views::transform([&coarseSimulationMesh](const auto fId) {
                const auto triangle = coarseSimulationMesh.vertices_around_face(coarseSimulationMesh.halfedge(fId)) |
                                      ranges::views::transform([](const auto vId) { return vId.idx(); }) |
                                      ranges::to_vector;
                return triangle;
            }) |
            ranges::to_vector;

        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::polygon_soup_to_polygon_mesh(points, fineTriangles, coarseMeshExport);

        std::ofstream out("coarse.obj");
        CGAL::print_polyhedron_wavefront(out, coarseMeshExport);
    }

    {
        CGAL::Polyhedron_3<VNCS::Space3D::K> fineMeshExport;

        std::vector<VNCS::Space2D::K::Point_3> points = fineMesh.vertices() |
                                                        ranges::views::transform([&fineMesh](const auto vId) {
                                                            const auto &p = fineMesh.point(vId);
                                                            return VNCS::Space2D::K::Point_3{p.x(), p.y(), p.z()};
                                                        }) |
                                                        ranges::to_vector;

        const auto fineTriangles = fineMesh.faces() | ranges::views::transform([&fineMesh](const auto fId) {
                                       const auto triangle =
                                           fineMesh.vertices_around_face(fineMesh.halfedge(fId)) |
                                           ranges::views::transform([](const auto vId) { return vId.idx(); }) |
                                           ranges::to_vector;
                                       return triangle;
                                   }) |
                                   ranges::to_vector;

        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::polygon_soup_to_polygon_mesh(points, fineTriangles, fineMeshExport);

        std::ofstream out("fine.obj");
        CGAL::print_polyhedron_wavefront(out, fineMeshExport);
    }

    // Generate interpolation data from coarseSimulationMesh to coarseMesh
    // This is needed to be able to animate the simulationMesh given an animation of the coarseMesh

    // for (const auto simVertexId : coarseSimulationMesh.vertices()) {
    //     // Find closest triangle... brute force, Im lazy
    //     auto coarseSimulationPoint = coarseSimulationMesh.point(simVertexId);
    //     std::vector<double> coordinates;

    //     CGAL::Barycentric_coordinates::mean_value_coordinates_3(coarseMesh, coarseSimulationPoint, coordinates);
    //     json interpolationForPoint = json::object();
    //     interpolationForPoint["coordinates"] = coordinates;
    //     interpolationForPoint["pointId"] = simVertexId.idx();

    //     interpolationJson.push_back(interpolationForPoint);
    // }

    json interpolationJson = json::array();
    for (const auto simVertexId : coarseSimulationMesh.vertices()) {
        // Find closest triangle... brute force, Im lazy
        auto distance = std::numeric_limits<double>::infinity();
        auto closestId = VNCS::Space3D::Mesh::Face_index{};
        auto coarseSimulationPoint = coarseSimulationMesh.point(simVertexId);
        for (const auto cTriangleId : coarseMesh.faces()) {
            auto trianglePoints =
                coarseMesh.vertices_around_face(coarseMesh.halfedge(cTriangleId)) |
                ranges::views::transform([&coarseMesh](const auto vId) { return coarseMesh.point(vId); }) |
                ranges::to_vector;
            VNCS::Space3D::Triangle triangle(trianglePoints[0], trianglePoints[1], trianglePoints[2]);

            const auto d = CGAL::squared_distance(coarseSimulationPoint, triangle);
            if (d < distance) {
                distance = d;
                closestId = cTriangleId;
            }
        }

        auto trianglePointsIds = coarseMesh.vertices_around_face(coarseMesh.halfedge(closestId)) | ranges::to_vector;

        auto trianglePoints =
            coarseMesh.vertices_around_face(coarseMesh.halfedge(closestId)) |
            ranges::views::transform([&coarseMesh](const auto vId) { return coarseMesh.point(vId); }) |
            ranges::to_vector;
        VNCS::Space3D::Triangle triangle(trianglePoints[0], trianglePoints[1], trianglePoints[2]);

        // This point is the point we are trying to animated before centroidal movement
        // We should compute the closest point from simVertexId to the triangle but lets try this first
        const auto coarsePoint = coarseMesh.point(simVertexId);

        // Compute barycentric coordinates of the projected point wrt the triangle
        std::vector<Eigen::Vector3d> trianglePointsEigen = trianglePoints | ranges::views::transform([](const auto p) {
                                                               return Eigen::Vector3d{p[0], p[1], p[2]};
                                                           }) |
                                                           ranges::to_vector;

        Eigen::Matrix3d m = Eigen::Matrix3d{{trianglePoints[0][0], trianglePoints[1][0], trianglePoints[2][0]},
                                            {trianglePoints[0][1], trianglePoints[1][1], trianglePoints[2][1]},
                                            {trianglePoints[0][2], trianglePoints[1][2], trianglePoints[2][2]}};

        const auto toEigen = [](const auto v) { return Eigen::Vector3d{v[0], v[1], v[2]}; };

        Eigen::Vector3d coordinates = m.inverse() * toEigen(coarsePoint);

        Eigen::Vector3d v0 = (trianglePointsEigen[1] - trianglePointsEigen[0]).normalized();
        Eigen::Vector3d v1 = toEigen(triangle.supporting_plane().orthogonal_vector()).normalized();
        Eigen::Vector3d v2 = v0.cross(v1);

        Eigen::Vector3d offset = toEigen(coarseSimulationPoint) - toEigen(coarsePoint);

        // Project the offset into the new base
        VNCS::Real offsetV0 = offset.dot(v0);
        VNCS::Real offsetV1 = offset.dot(v1);
        VNCS::Real offsetV2 = offset.dot(v2);

        json interpolationForPoint = json::object();
        interpolationForPoint["a"] = coordinates[0];
        interpolationForPoint["aId"] = trianglePointsIds[0].idx();
        interpolationForPoint["b"] = coordinates[1];
        interpolationForPoint["bId"] = trianglePointsIds[1].idx();
        interpolationForPoint["c"] = coordinates[2];
        interpolationForPoint["cId"] = trianglePointsIds[2].idx();
        interpolationForPoint["closesTriangleId"] = closestId.idx();
        interpolationForPoint["pointId"] = simVertexId.idx();
        interpolationForPoint["Objective_Frame0"] =
            std::array<VNCS::Real, 3>{coarseSimulationPoint[0], coarseSimulationPoint[1], coarseSimulationPoint[2]};
        interpolationForPoint["v0_axis"] = v0;
        interpolationForPoint["v1_axis"] = v1;
        interpolationForPoint["v2_axis"] = v2;

        interpolationForPoint["v0"] = offsetV0;
        interpolationForPoint["v1"] = offsetV1;
        interpolationForPoint["v2"] = offsetV2;

        interpolationJson.push_back(interpolationForPoint);
    }

    std::ofstream oNodalFine("interpolation.json");
    oNodalFine << std::setw(4) << interpolationJson << std::endl;

    SamplingPointsExporter samplers;
    samplers.setCoarseSamplersFilePath(m_coarseSamplersFilePath);
    samplers.setFineSamplersFilePath(m_fineSamplersFilePath);
    samplers(coarseSimulationMesh, fineMesh);
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

std::optional<VNCS::Generators::RemeshCriteria> Generator::coarseCriteria() const
{
    return m_coarseCriteria;
}

void Generator::setCoarseCriteria(const VNCS::Generators::RemeshCriteria &coarseCriteria)
{
    m_coarseCriteria = coarseCriteria;
}

std::optional<VNCS::Generators::RemeshCriteria> Generator::fineCriteria() const
{
    return m_fineCriteria;
}

void Generator::setFineCriteria(const VNCS::Generators::RemeshCriteria &fineCriteria)
{
    m_fineCriteria = fineCriteria;
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

}  // namespace Gen223
}  // namespace Generators
}  // namespace VNCS
