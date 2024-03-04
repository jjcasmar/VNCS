#include "VisualSamplersExporter.h"
#include <VNCS/DataExtensions.h>

#include <VNCS/IntegrationScheme.h>

#include <nlohmann/json.hpp>
#include <fstream>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>

#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/transform.hpp>

#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>
#include "PMPTriangleIntersection.h"

#include <VNCS/SamplingPoints.h>

using json = nlohmann::json;

namespace
{
using Mesh = CGAL::Surface_mesh<VNCS::Sim2D::Kernel::Point>;
using Point = VNCS::Sim2D::Kernel::Point;
using Triangle = VNCS::Sim2D::Kernel::K::Triangle_2;
using Real = VNCS::Sim2D::Kernel::Real;
using BarycentricCoordinates = CGAL::Barycentric_coordinates::Triangle_coordinates_2<VNCS::Sim2D::Kernel::K>;

}  // namespace

namespace VNCS
{
namespace Sim2D
{
VisualSamplersExporter::VisualSamplersExporter()
{
}

void VisualSamplersExporter::operator()(const Mesh &visualMesh, const Mesh &coarseMesh, const Mesh &fineMesh)
{
    {
        VNCS::Sim2D::PMPTriangleIntersection fineIntersectionFunc(fineMesh);
        VNCS::Sim2D::PMPTriangleIntersection coarseIntersectionFunc(coarseMesh);

        json visualSamplers = json::array();

        for (const auto &[vdId, vd] : visualMesh.vertices() | ranges::views::enumerate) {
            const auto &point = visualMesh.point(vd);

            // Find in which coarse triangle this visual vertex is, if any
            const auto coarseTriangleId = coarseIntersectionFunc.query(visualMesh.point(vd), true);

            // Find in which fine triangle this visual vertex is, if any
            const auto fineTriangleId = fineIntersectionFunc.query(visualMesh.point(vd));

            const auto alphaInTriangle = [this, &fineMesh, fineTriangleId](const auto &p) {
                const auto vIndices = fineMesh.vertices_around_face(fineMesh.halfedge(fineTriangleId.value()));
                const auto v = vIndices |
                               ranges::views::transform([&fineMesh](const auto vd) { return fineMesh.point(vd); }) |
                               ranges::to_vector;
                std::array<Real, 3> blendingFactors = {
                    m_blendingField->blending(v[0]), m_blendingField->blending(v[1]), m_blendingField->blending(v[2])};

                std::array<Real, 3> coordinates;

                BarycentricCoordinates coordinatesComputation(v[0],  //
                                                              v[1],
                                                              v[2]);

                coordinatesComputation(p, std::begin(coordinates));

                const auto blendingValue = std::inner_product(std::begin(blendingFactors),  //
                                                              std::end(blendingFactors),
                                                              std::begin(coordinates),
                                                              0.0);
                return blendingValue;
            };

            VNCS::Sampler<VNCS::Space2D> sampler;
            sampler.x = {point[0], point[1]};

            sampler.a = fineTriangleId ? alphaInTriangle(point) : 1.0;
            sampler.coarseFaceIdx = coarseTriangleId;
            sampler.fineFaceIdx = fineTriangleId;

            if (coarseTriangleId) {
                std::vector<ShapeFunction<VNCS::Space2D>> coarseShapeFunctions;
                const auto vIndices =
                    coarseMesh.vertices_around_face(coarseMesh.halfedge(coarseTriangleId.value())) | ranges::to_vector;
                const auto v = vIndices |
                               ranges::views::transform([&coarseMesh](const auto vd) { return coarseMesh.point(vd); }) |
                               ranges::to_vector;

                std::array<Real, 3> coordinates;
                BarycentricCoordinates coordinatesComputation(v[0],  //
                                                              v[1],
                                                              v[2]);

                coordinatesComputation(point, std::begin(coordinates));

                for (int i = 0; i < 3; ++i) {
                    ShapeFunction<VNCS::Space2D> coarseShapeFunction;

                    coarseShapeFunction.nodeIndex = vIndices[i].idx();
                    coarseShapeFunction.v = coordinates[i];
                    coarseShapeFunctions.push_back(coarseShapeFunction);
                }
                sampler.coarseShapeFunctions = coarseShapeFunctions;
            }

            if (fineTriangleId) {
                std::vector<ShapeFunction<VNCS::Space2D>> fineShapeFunctions;
                const auto vIndices =
                    fineMesh.vertices_around_face(fineMesh.halfedge(fineTriangleId.value())) | ranges::to_vector;
                const auto v = vIndices |
                               ranges::views::transform([&fineMesh](const auto vd) { return fineMesh.point(vd); }) |
                               ranges::to_vector;

                std::array<Real, 3> coordinates;
                BarycentricCoordinates coordinatesComputation(v[0],  //
                                                              v[1],
                                                              v[2]);

                coordinatesComputation(point, std::begin(coordinates));

                for (int i = 0; i < 3; ++i) {
                    ShapeFunction<VNCS::Space2D> fineShapeFunction;

                    fineShapeFunction.nodeIndex = vIndices[i].idx();
                    fineShapeFunction.v = coordinates[i];
                    fineShapeFunctions.push_back(fineShapeFunction);
                }
                sampler.fineShapeFunctions = fineShapeFunctions;
            }

            visualSamplers.push_back(sampler);
        }

        json visualSamplersJson = visualSamplers;

        std::ofstream oCoarse(m_samplersFilePath);
        oCoarse << std::setw(4) << visualSamplersJson << std::endl;
    }
}

const std::filesystem::path &VisualSamplersExporter::samplersFilePath() const
{
    return m_samplersFilePath;
}

void VisualSamplersExporter::setSamplersFilePath(const std::filesystem::path &filePath)
{
    m_samplersFilePath = filePath;
}

void VisualSamplersExporter::setBlendingField(std::shared_ptr<VNCS::BlendingField<Real, K::Point_2>> blendingField)
{
    m_blendingField = blendingField;
}

}  // namespace Sim2D
}  // namespace VNCS
