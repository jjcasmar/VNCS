#include "Gen202.h"

#include <VNCS/Spaces.h>
#include <VNCS/Generators/CDT.h>
#include <VNCS/Generators/AdaptiveMeshCriteria.h>

#include <VNCS/Generators/Gen202/SamplingPointsExporter.h>

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

using json = nlohmann::json;
namespace PMP = CGAL::Polygon_mesh_processing;

namespace VNCS
{
namespace Generators
{
namespace Gen202
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

    // xNIndices are all the other vertices
    const auto xNIndices = coarseMesh.vertices() | ranges::to_vector;

    // Create xN degrees of freedom
    std::vector<VNCS::Space2D::Point> xN;
    xN.reserve(xNIndices.size());
    for (auto xNVertexId : xNIndices) {
        const auto p = coarseMesh.point(xNVertexId);
        xN.push_back({p.x(), p.y()});
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

        dofJson["x0"] = x0;

        std::ofstream dof(m_dofFilePath);
        dof << std::setw(4) << dofJson << std::endl;
    }

    const VNCS::Space2D::Mesh coarseSimulationMesh = coarseMesh;

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
    samplers(coarseSimulationMesh);
}

std::filesystem::path Generator::coarseMeshPath() const
{
    return m_coarseMeshPath;
}

void Generator::setCoarseMeshPath(const std::filesystem::path &coarseMeshPath)
{
    m_coarseMeshPath = coarseMeshPath;
}

std::filesystem::path Generator::coarseSamplersFilePath() const
{
    return m_coarseSamplersFilePath;
}

void Generator::setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath)
{
    m_coarseSamplersFilePath = coarseSamplersFilePath;
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
}  // namespace Gen202
}  // namespace Generators
}  // namespace VNCS
