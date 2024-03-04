#include "SamplingPointsExporter.h"

#include <VNCS/Generators/Gen222/GenerateSampler.h>
#include <VNCS/Generators/PMPTriangleIntersection.h>
#include <VNCS/Generators/IntegrationScheme.h>

#include <range/v3/view/transform.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/to_container.hpp>

using Triangle = VNCS::Space2D::Triangle;

void VNCS::Generators::Gen222::SamplingPointsExporter::operator()(const Mesh &coarseMesh,
                                                                  const Mesh &fineMesh,
                                                                  std::size_t gridSize)
{
    const VNCS::Generators::PMPTriangleIntersection coarseIntersectionFunc(coarseMesh);
    const VNCS::Generators::PMPTriangleIntersection fineIntersectionFunc(fineMesh);
    VNCS::Generators::Gen222::GenerateSampler generateSampler(coarseIntersectionFunc, fineIntersectionFunc);

    for (const auto &[coarseFaceId, coarseFace] : coarseMesh.faces() | ranges::views::enumerate) {
        // Get the indices and points of the face
        const auto vIndices = coarseMesh.vertices_around_face(coarseMesh.halfedge(coarseFace)) | ranges::to_vector;
        const auto v = vIndices |
                       ranges::views::transform([&coarseMesh](const auto &vd) { return coarseMesh.point(vd); }) |
                       ranges::to_vector;

        // Create the J matrix. This matrix allows to convert the coordinates of the sampler point to an actual
        // position
        // and to compute derivatives
        Eigen::Matrix3d J;
        J << v[1][0] - v[0][0], v[2][0] - v[0][0], v[0][0],  //
            v[1][1] - v[0][1], v[2][1] - v[0][1], v[0][1],   //
            0, 0, 1;

        const Triangle coarseTriangle = {v[0], v[1], v[2]};
        struct MCRegion {
            std::vector<Point> points;
            std::optional<Mesh::Face_index> coarseFace;
        };

        auto bbox = coarseTriangle.bbox();
        auto xSize = bbox.xmax() - bbox.xmin();
        auto ySize = bbox.ymax() - bbox.ymin();

        std::vector<Point> mcPoints;
        const auto nbPoints = gridSize;

        mcPoints.reserve(nbPoints * nbPoints);
        const auto xStep = xSize / static_cast<Real>(nbPoints);
        const auto yStep = ySize / static_cast<Real>(nbPoints);
        for (int i = 0; i < nbPoints; ++i) {
            for (int j = 0; j < nbPoints; ++j) {
                mcPoints.emplace_back(bbox.xmin() + i * xStep, bbox.ymin() + j * yStep);
            }
        }

        // Remove points that are not on the triangle
        mcPoints = mcPoints | ranges::views::filter([&coarseTriangle](const auto &p) {
                       return !coarseTriangle.has_on_unbounded_side(p);
                   }) |
                   ranges::to_vector;
        const auto estimateOverlapedVolume =
            ranges::distance(mcPoints | ranges::views::filter([&fineIntersectionFunc](const Point &sampler) {
                                 return !fineIntersectionFunc.query(sampler).has_value();
                             }));

        if (estimateOverlapedVolume) {
            SamplerDefinition samplerDefinition;
            samplerDefinition.coordinates = {1.0 / 3.0, 1.0 / 3.0};
            samplerDefinition.baseFace = coarseFace;
            samplerDefinition.locationType = SamplerDefinition::LocationType::CoarseTriangle;
            samplerDefinition.kinematicType = SamplerDefinition::KinematicType::Coarse;
            samplerDefinition.weight = coarseTriangle.area() *  //
                                       static_cast<Real>(estimateOverlapedVolume) / static_cast<Real>(mcPoints.size());
            m_coarseSamplingPoints.push_back(
                generateSampler(samplerDefinition,  //
                                coarseMesh,
                                VNCS::Generators::Gen222::GenerateSampler::CoarseSampler{}));
        }
    }

    for (const auto [fineFaceId, fineFace] : fineMesh.faces() | ranges::views::enumerate) {
        const auto vIndices = fineMesh.vertices_around_face(fineMesh.halfedge(fineFace)) | ranges::to_vector;
        const auto v = vIndices | ranges::views::transform([&fineMesh](const auto &vd) { return fineMesh.point(vd); }) |
                       ranges::to_vector;

        const Triangle fineTriangle{v[0], v[1], v[2]};
        SamplerDefinition samplerDefinition;
        samplerDefinition.coordinates = {1.0 / 3.0, 1.0 / 3.0};
        samplerDefinition.baseFace = fineFace;
        samplerDefinition.locationType = SamplerDefinition::LocationType::FineTriangle;
        samplerDefinition.kinematicType = SamplerDefinition::KinematicType::Fine;
        samplerDefinition.weight = fineTriangle.area();
        m_fineSamplingPoints.push_back(generateSampler(samplerDefinition,  //
                                                       fineMesh,
                                                       VNCS::Generators::Gen222::GenerateSampler::FineSampler{}));
    }

    {
        json coarseSamplers = m_coarseSamplingPoints;
        std::ofstream oCoarse(m_coarseSamplersFilePath);
        oCoarse << std::setw(4) << coarseSamplers << std::endl;
    }

    {
        json fineSamplers = m_fineSamplingPoints;
        std::ofstream oFine(m_fineSamplersFilePath);
        oFine << std::setw(4) << fineSamplers << std::endl;
    }
}

std::filesystem::path VNCS::Generators::Gen222::SamplingPointsExporter::coarseSamplersFilePath() const
{
    return m_coarseSamplersFilePath;
}

void VNCS::Generators::Gen222::SamplingPointsExporter::setCoarseSamplersFilePath(
    const std::filesystem::path &coarseSamplersFilePath)
{
    m_coarseSamplersFilePath = coarseSamplersFilePath;
}

std::filesystem::path VNCS::Generators::Gen222::SamplingPointsExporter::fineSamplersFilePath() const
{
    return m_fineSamplersFilePath;
}

void VNCS::Generators::Gen222::SamplingPointsExporter::setFineSamplersFilePath(
    const std::filesystem::path &fineSamplersFilePath)
{
    m_fineSamplersFilePath = fineSamplersFilePath;
}

std::filesystem::path VNCS::Generators::Gen222::SamplingPointsExporter::fineNodeSamplersFilePath() const
{
    return m_fineNodeSamplersFilePath;
}

void VNCS::Generators::Gen222::SamplingPointsExporter::setFineNodeSamplersFilePath(
    const std::filesystem::path &fineNodeSamplersFilePath)
{
    m_fineNodeSamplersFilePath = fineNodeSamplersFilePath;
}
