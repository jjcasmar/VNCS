#include "VisualSamplersExporter.h"

#include <VNCS/Generators/TetraMeshIntersection.h>
#include <VNCS/SamplingPoints.h>

#include <spdlog/spdlog.h>
#include <nlohmann/json.hpp>

#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>

#include <Eigen/Sparse>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

using json = nlohmann::json;

std::pair<Eigen::SparseMatrix<VNCS::Real>, Eigen::SparseMatrix<VNCS::Real>> VNCS::Generators::Gen333::VisualSamplersExporter::matrices() const
{
    auto logger = spdlog::get("Gen333");
    VNCS::Generators::TetraMeshIntersection fineIntersectionFunc(m_fineMesh);
    VNCS::Generators::TetraMeshIntersection coarseIntersectionFunc(m_coarseMesh);

    json visualSamplers = json::array();

    std::vector<VNCS::Real> coarseWeights = m_coarseMesh.points | ranges::views::transform([this](const auto &p) {
        return m_blendingField->blending(p);
    }) | ranges::to_vector;

    std::vector<VNCS::Real> fineWeights = m_fineMesh.points | ranges::views::transform([this](const auto &p) {
        return m_blendingField->blending(p);
    }) | ranges::to_vector;

    std::vector<Eigen::Triplet<VNCS::Real>> coarseTriplets;
    std::vector<Eigen::Triplet<VNCS::Real>> fineTriplets;

    for (const auto &[vdId, vd] : m_visualMesh.vertices() | ranges::views::enumerate) {
        const auto &point = m_visualMesh.point(vd);

        auto blendingValue = m_blendingField->blending(point);

        // Find in which fine tetra this visual vertex is, if any
        // We only look for a fine tetrahedron if we are in the hybrid or fine region!
        const auto fineTetraId = blendingValue < 1.0 ? fineIntersectionFunc.query(point, true) : std::optional<int>{};

        const auto alphaInTriangle = [this, vdId, fineTetraId, &fineWeights](const auto &p) {
            const auto vIndices = m_fineMesh.tetras[fineTetraId.value()];
            const auto v = vIndices |
                           ranges::views::transform([this](const auto vd) { return m_fineMesh.points[vd]; }) |
                           ranges::to_vector;
            const auto blendingFactors = vIndices |
                                         ranges::views::transform([&fineWeights](const auto &vd) { return fineWeights[vd]; }) |
                                         ranges::to_vector;

            std::array<Real, 4> coordinates;

            const VNCS::Space3D::Tetra tetrahedron(v[0],  //
                                                   v[1],
                                                   v[2],
                                                   v[3]);

            Eigen::Matrix4d barMat;
            barMat << v[0].x(), v[1].x(), v[2].x(), v[3].x(),  //
                v[0].y(), v[1].y(), v[2].y(), v[3].y(),        //
                v[0].z(), v[1].z(), v[2].z(), v[3].z(),        //
                1, 1, 1, 1;

            Eigen::Map<Eigen::Vector4d> coordinatesMap(coordinates.data());
            coordinatesMap = barMat.inverse() * Eigen::Vector4d(p.x(), p.y(), p.z(), 1.0);

            const auto blendingValue = std::inner_product(std::begin(blendingFactors),  //
                                                          std::end(blendingFactors),
                                                          std::begin(coordinates),
                                                          0.0);
            return blendingValue;
        };

        // Find in which coarse tetra this visual vertex is, if any
        // We only look for a coarse tetrahedron if we are in the hybrid or coarse region!
        const auto alpha = fineTetraId ? alphaInTriangle(point) : 1.0;
        const auto coarseTetraId = alpha > 0.0 ? coarseIntersectionFunc.query(point, true) : std::optional<int>{};

        if (coarseTetraId) {
            std::vector<VNCS::ShapeFunction<VNCS::Space3D>> coarseShapeFunctions;
            const auto vIndices = m_coarseMesh.tetras[coarseTetraId.value()];
            const auto v = vIndices |
                           ranges::views::transform([this](const auto vd) { return m_coarseMesh.points[vd]; }) |
                           ranges::to_vector;

            std::array<Real, 4> coordinates;

            const VNCS::Space3D::Tetra tetrahedron(v[0],  //
                                                   v[1],
                                                   v[2],
                                                   v[3]);

            Eigen::Matrix4d barMat;
            barMat << v[0].x(), v[1].x(), v[2].x(), v[3].x(),  //
                v[0].y(), v[1].y(), v[2].y(), v[3].y(),        //
                v[0].z(), v[1].z(), v[2].z(), v[3].z(),        //
                1, 1, 1, 1;

            Eigen::Map<Eigen::Vector4d> coordinatesMap(coordinates.data());
            coordinatesMap = barMat.inverse() * Eigen::Vector4d(point.x(), point.y(), point.z(), 1.0);

            coarseTriplets.emplace_back(3*vdId + 0, 3*vIndices[0] + 0, alpha * coordinates[0]);
            coarseTriplets.emplace_back(3*vdId + 1, 3*vIndices[0] + 1, alpha * coordinates[0]);
            coarseTriplets.emplace_back(3*vdId + 2, 3*vIndices[0] + 2, alpha * coordinates[0]);

            coarseTriplets.emplace_back(3*vdId + 0, 3*vIndices[1] + 0, alpha * coordinates[1]);
            coarseTriplets.emplace_back(3*vdId + 1, 3*vIndices[1] + 1, alpha * coordinates[1]);
            coarseTriplets.emplace_back(3*vdId + 2, 3*vIndices[1] + 2, alpha * coordinates[1]);

            coarseTriplets.emplace_back(3*vdId + 0, 3*vIndices[2] + 0, alpha * coordinates[2]);
            coarseTriplets.emplace_back(3*vdId + 1, 3*vIndices[2] + 1, alpha * coordinates[2]);
            coarseTriplets.emplace_back(3*vdId + 2, 3*vIndices[2] + 2, alpha * coordinates[2]);

            coarseTriplets.emplace_back(3*vdId + 0, 3*vIndices[3] + 0, alpha * coordinates[3]);
            coarseTriplets.emplace_back(3*vdId + 1, 3*vIndices[3] + 1, alpha * coordinates[3]);
            coarseTriplets.emplace_back(3*vdId + 2, 3*vIndices[3] + 2, alpha * coordinates[3]);
        }

        if (fineTetraId) {
            std::vector<VNCS::ShapeFunction<VNCS::Space3D>> fineShapeFunctions;
            const auto vIndices = m_fineMesh.tetras[fineTetraId.value()];
            const auto v = vIndices |
                           ranges::views::transform([this](const auto vd) { return m_fineMesh.points[vd]; }) |
                           ranges::to_vector;

            std::array<Real, 4> coordinates;

            const VNCS::Space3D::Tetra tetrahedron(v[0],  //
                                                   v[1],
                                                   v[2],
                                                   v[3]);

            Eigen::Matrix4d barMat;
            barMat << v[0].x(), v[1].x(), v[2].x(), v[3].x(),  //
                v[0].y(), v[1].y(), v[2].y(), v[3].y(),        //
                v[0].z(), v[1].z(), v[2].z(), v[3].z(),        //
                1, 1, 1, 1;

            Eigen::Map<Eigen::Vector4d> coordinatesMap(coordinates.data());
            coordinatesMap = barMat.inverse() * Eigen::Vector4d(point.x(), point.y(), point.z(), 1.0);

            fineTriplets.emplace_back(3*vdId + 0, 3*vIndices[0] + 0, (1-alpha) * coordinates[0]);
            fineTriplets.emplace_back(3*vdId + 1, 3*vIndices[0] + 1, (1-alpha) * coordinates[0]);
            fineTriplets.emplace_back(3*vdId + 2, 3*vIndices[0] + 2, (1-alpha) * coordinates[0]);

            fineTriplets.emplace_back(3*vdId + 0, 3*vIndices[1] + 0, (1-alpha) * coordinates[1]);
            fineTriplets.emplace_back(3*vdId + 1, 3*vIndices[1] + 1, (1-alpha) * coordinates[1]);
            fineTriplets.emplace_back(3*vdId + 2, 3*vIndices[1] + 2, (1-alpha) * coordinates[1]);

            fineTriplets.emplace_back(3*vdId + 0, 3*vIndices[2] + 0, (1-alpha) * coordinates[2]);
            fineTriplets.emplace_back(3*vdId + 1, 3*vIndices[2] + 1, (1-alpha) * coordinates[2]);
            fineTriplets.emplace_back(3*vdId + 2, 3*vIndices[2] + 2, (1-alpha) * coordinates[2]);

            fineTriplets.emplace_back(3*vdId + 0, 3*vIndices[3] + 0, (1-alpha) * coordinates[3]);
            fineTriplets.emplace_back(3*vdId + 1, 3*vIndices[3] + 1, (1-alpha) * coordinates[3]);
            fineTriplets.emplace_back(3*vdId + 2, 3*vIndices[3] + 2, (1-alpha) * coordinates[3]);

        }
    }

    Eigen::SparseMatrix<VNCS::Real> coarseMatrix;
    coarseMatrix.resize(3*m_visualMesh.num_vertices(), 3*m_coarseMesh.points.size());
    coarseMatrix.setFromTriplets(coarseTriplets.begin(), coarseTriplets.end());

    Eigen::SparseMatrix<VNCS::Real> fineMatrix;
    fineMatrix.resize(3*m_visualMesh.num_vertices(), 3*m_fineMesh.points.size());
    fineMatrix.setFromTriplets(fineTriplets.begin(), fineTriplets.end());

    return std::make_pair(coarseMatrix, fineMatrix);
}

void VNCS::Generators::Gen333::VisualSamplersExporter::setBlendingField(
    std::shared_ptr<BlendingField<VNCS::Space3D::Real, VNCS::Space3D::Point>> blendingField)
{
    m_blendingField = blendingField;
}
void VNCS::Generators::Gen333::VisualSamplersExporter::setVisualMeshPath(const std::filesystem::path &visualPath)
{
    {
        std::vector<VNCS::Space3D::Point> points;
        std::vector<std::vector<std::size_t>> faces;
        std::ifstream coarseMeshIn(visualPath);
        CGAL::read_OBJ(coarseMeshIn, points, faces);

        namespace PMP = CGAL::Polygon_mesh_processing;
        PMP::polygon_soup_to_polygon_mesh(points, faces, m_visualMesh);
    }
}
void VNCS::Generators::Gen333::VisualSamplersExporter::setFineTetraMesh(const VNCS::Space3D::TetraMesh &mesh)
{
    m_fineMesh = mesh;
}
void VNCS::Generators::Gen333::VisualSamplersExporter::setCoarseTetraMesh(const VNCS::Space3D::TetraMesh &mesh)
{
    m_coarseMesh = mesh;
}
