#include "Gen212.h"

#include <VNCS/Spaces.h>
#include <VNCS/Generators/CDT.h>
#include <VNCS/Generators/AdaptiveMeshCriteria.h>

#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/polygon_mesh_processing.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/print_wavefront.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <range/v3/view/subrange.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/action/sort.hpp>
#include <range/v3/action/adjacent_remove_if.hpp>
#include <range/v3/algorithm/find.hpp>
#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/algorithm/contains.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/group_by.hpp>

#include <algorithm>

#include <spdlog/spdlog.h>
#include <unsupported/Eigen/SparseExtra>
#include <nlohmann/json.hpp>

#include <VNCS/Generators/Gen212/SamplingPointsExporter.h>
#include <VNCS/Logger.h>
#include <VNCS/EdgeMesh.h>

using json = nlohmann::json;
namespace PMP = CGAL::Polygon_mesh_processing;

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

}  // namespace

namespace VNCS
{
namespace Generators
{
namespace Gen212
{
void Generator::operator()()
{
    // Load OBJ

    VNCS::Space2D::Mesh coarseInputMesh;

    {
        std::vector<VNCS::Space2D::Point> coarsePoints;
        std::vector<std::vector<std::size_t>> coarseFaces;
        std::ifstream coarseMeshIn(m_coarseMeshPath);
        CGAL::read_OBJ(coarseMeshIn, coarsePoints, coarseFaces);

        PMP::polygon_soup_to_polygon_mesh(coarsePoints, coarseFaces, coarseInputMesh);
    }

    const auto coarseMesh = [this, &coarseInputMesh]() {
        if (m_remeshCoarseMesh) {
            std::cout << "remeshing mesh\n";
            return VNCS::Generators::cdt(coarseInputMesh, AdaptiveMeshCriteria<CDT>{m_coarseCriteria});
        } else
            return coarseInputMesh;
    }();

    const VNCS::EdgeMesh<VNCS::Space2D> fineMesh = [this]() {
        std::vector<VNCS::Space2D::Point> fineMeshPoints;
        std::vector<std::vector<std::size_t>> coarseLines;
        std::ifstream fineMeshIn(m_edgeMeshPath);
        VNCS::read_OBJ(fineMeshIn, fineMeshPoints, coarseLines);

        return VNCS::lines_soup_to_edge_mesh<VNCS::Space2D>(fineMeshPoints, coarseLines);
    }();

    const auto &fineMeshPointMap = boost::get(PointVertexTag{}, fineMesh);
    const auto &[fineVerticesBegin, fineVerticesEnd] = boost::vertices(fineMesh);
    const auto fineVerticesSubrange = ranges::subrange(fineVerticesBegin, fineVerticesEnd);

    // Get the indices of the possible cluster points. This are the non boundary vertices
    const auto possibleClusterPointsIndices =
        coarseMesh.vertices() |                                                                    //
        ranges::views::filter([&coarseMesh](const auto v) { return !coarseMesh.is_border(v); }) |  //
        ranges::to_vector;

    const auto possibleClusterPoints =
        possibleClusterPointsIndices |                                                               //
        ranges::views::transform([&coarseMesh](const auto vId) { return coarseMesh.point(vId); }) |  //
        ranges::to_vector;

    ::My_point_property_map ppmap(possibleClusterPoints);
    Tree searchTree(boost::counting_iterator<std::size_t>(0),
                    boost::counting_iterator<std::size_t>(possibleClusterPoints.size()),
                    Splitter(),
                    TreeTraits(ppmap));
    Distance trDistance(ppmap);

    // xC points are only those that are required because they are the closest to a fine point
    auto xCIndices =
        fineVerticesSubrange |
        ranges::views::transform([&fineMeshPointMap](const auto vertex) { return fineMeshPointMap[vertex]; }) |
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

    std::vector<decltype(xCIndices)::value_type> assignedXCvertices(boost::num_vertices(fineMesh));
    for (int i = 0; i < boost::num_vertices(fineMesh); ++i) {
        const auto &dPoint = fineMeshPointMap[i];
        Neighbor_search search2(searchTree2, dPoint, 1, 0.0001, true, trDistance2);
        for (const auto searchResult : search2) {
            assignedXCvertices[i] = xCIndices[searchResult.first];
        }
    }

    // xNIndices are all the other vertices
    const auto xNIndices =
        coarseMesh.vertices() |
        ranges::views::filter([&xCIndices](const auto vId) { return !ranges::contains(xCIndices, vId); }) |
        ranges::to_vector;

    std::vector<Eigen::Triplet<VNCS::Space2D::Real>> triplets;
    triplets.reserve(2 * boost::num_vertices(fineMesh));
    for (int i = 0; i < boost::num_vertices(fineMesh); ++i) {
        const auto point = fineMeshPointMap[i];
        auto it = std::find(std::begin(xCIndices), std::end(xCIndices), assignedXCvertices[i]);
        auto index = std::distance(std::begin(xCIndices), it);
        triplets.push_back(Eigen::Triplet<VNCS::Space2D::Real>(2 * i + 0, 2 * index + 0, 1.0));
        triplets.push_back(Eigen::Triplet<VNCS::Space2D::Real>(2 * i + 1, 2 * index + 1, 1.0));
    }

    // Create xN degrees of freedom
    std::vector<VNCS::Space2D::Point> xN;
    xN.reserve(xNIndices.size());
    for (auto xNVertexId : xNIndices) {
        const auto p = coarseMesh.point(xNVertexId);
        xN.push_back({p.x(), p.y()});
    }

    {
        CGAL::Delaunay_triangulation_2<VNCS::Space2D::K> triangulation;
        for (const auto &point : ranges::views::concat(xN, xC))
            triangulation.insert(point);

        const VNCS::Space2D::Mesh coarseSimulationMesh = VNCS::Generators::meshFromTriangulation(triangulation);
        CGAL::Polyhedron_3<VNCS::Space2D::K> coarseSimulationMeshExport;

        std::vector<VNCS::Space2D::K::Point_3> points =
            coarseSimulationMesh.vertices() | ranges::views::transform([&coarseSimulationMesh](const auto vId) {
                const auto &p = coarseSimulationMesh.point(vId);
                return VNCS::Space2D::K::Point_3{p.x(), p.y(), 0.0};
            }) |
            ranges::to_vector;

        const auto coarseTriangles =
            coarseSimulationMesh.faces() | ranges::views::transform([&coarseSimulationMesh](const auto fId) {
                const auto triangle = coarseSimulationMesh.vertices_around_face(coarseSimulationMesh.halfedge(fId)) |
                                      ranges::views::transform([](const auto vId) { return vId.idx(); }) |
                                      ranges::to_vector;
                return triangle;
            }) |
            ranges::to_vector;

        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::polygon_soup_to_polygon_mesh(points, coarseTriangles, coarseSimulationMeshExport);

        std::ofstream out("coarse_init.obj");
        CGAL::print_polyhedron_wavefront(out, coarseSimulationMeshExport);
    }

    Eigen::SparseMatrix<VNCS::Space2D::Real> C;
    C.resize(2 * boost::num_vertices(fineMesh),  //
             2 * xCIndices.size());
    C.setFromTriplets(std::begin(triplets), std::end(triplets), [](const auto &a, const auto &b) {
        spdlog::get("VNCS")->warn("Repeated triplets shouldnt happen");
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

    Eigen::Matrix<VNCS::Space2D::Real, Eigen::Dynamic, 1> dEig;
    auto d = fineVerticesSubrange;
    dEig.resize(2 * boost::num_vertices(fineMesh));
    for (const auto &[index, vertexHandle] : ranges::views::enumerate(d)) {
        const auto &p = fineMeshPointMap[vertexHandle];
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

    std::vector<VNCS::Space2D::Point> dValues;
    for (const auto &[index, vertexHandle] : ranges::views::enumerate(fineVerticesSubrange))
        dValues.push_back(VNCS::Space2D::Point{dEig[2 * index + 0], dEig[2 * index + 1]});

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
        for (const auto &p : dValues) {
            auto v = json::array();
            v.push_back(p[0]);
            v.push_back(p[1]);
            dJson.push_back(v);
        }
        x0["d"] = dJson;

        dofJson["x0"] = x0;

        std::ofstream dof(m_dofFilePath);
        dof << std::setw(4) << dofJson << std::endl;
    }

    CGAL::Delaunay_triangulation_2<VNCS::Space2D::K> triangulation;
    for (const auto &point : ranges::views::concat(xN, xC))
        triangulation.insert(point);

    const VNCS::Space2D::Mesh coarseSimulationMesh = VNCS::Generators::meshFromTriangulation(triangulation);

    {
        CGAL::Polyhedron_3<VNCS::Space2D::K> coarseSimulationMeshExport;

        std::vector<VNCS::Space2D::K::Point_3> points =
            coarseSimulationMesh.vertices() | ranges::views::transform([&coarseSimulationMesh](const auto vId) {
                const auto &p = coarseSimulationMesh.point(vId);
                return VNCS::Space2D::K::Point_3{p.x(), p.y(), 0.0};
            }) |
            ranges::to_vector;

        const auto coarseTriangles =
            coarseSimulationMesh.faces() | ranges::views::transform([&coarseSimulationMesh](const auto fId) {
                const auto triangle = coarseSimulationMesh.vertices_around_face(coarseSimulationMesh.halfedge(fId)) |
                                      ranges::views::transform([](const auto vId) { return vId.idx(); }) |
                                      ranges::to_vector;
                return triangle;
            }) |
            ranges::to_vector;

        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::polygon_soup_to_polygon_mesh(points, coarseTriangles, coarseSimulationMeshExport);

        std::ofstream out("coarse.obj");
        CGAL::print_polyhedron_wavefront(out, coarseSimulationMeshExport);
    }

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

std::filesystem::path Generator::edgeMeshPath() const
{
    return m_edgeMeshPath;
}

void Generator::setEdgeMeshPath(const std::filesystem::path &edgeMeshPath)
{
    m_edgeMeshPath = edgeMeshPath;
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

std::shared_ptr<AdaptiveMeshCriteriaFunctor> Generator::coarseCriteria() const
{
    return m_coarseCriteria;
}

void Generator::setCoarseCriteria(const std::shared_ptr<AdaptiveMeshCriteriaFunctor> &coarseCriteria)
{
    m_coarseCriteria = coarseCriteria;
}

bool Generator::remeshCoarseMesh() const
{
    return m_remeshCoarseMesh;
}

void Generator::setRemeshCoarseMesh(bool newRemeshCoarseMesh)
{
    m_remeshCoarseMesh = newRemeshCoarseMesh;
}
}  // namespace Gen212
}  // namespace Generators
}  // namespace VNCS
