#include "SamplingPointsExporter.h"

#include <VNCS/Generators/Gen202/GenerateSampler.h>

#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>

namespace VNCS
{
namespace Generators
{
namespace Gen202
{
void SamplingPointsExporter::operator()(const VNCS::Space2D::Mesh &coarseMesh)
{
    GenerateSampler generateSampler;
    for (const auto [coarseFaceId, coarseFace] : coarseMesh.faces() | ranges::views::enumerate) {
        // Get the indices and points of the face
        const auto v = coarseMesh.vertices_around_face(coarseMesh.halfedge(coarseFace)) |
                       ranges::views::transform([&coarseMesh](const auto &vd) { return coarseMesh.point(vd); }) |
                       ranges::to_vector;

        const VNCS::Space2D::Triangle coarseTriangle = {v[0], v[1], v[2]};

        CoarseSamplerDefinition samplerDefinition;
        samplerDefinition.coordinates = {1.0 / 3.0, 1.0 / 3.0};
        samplerDefinition.baseFace = coarseFace;
        samplerDefinition.weight = coarseTriangle.area();
        m_coarseSamplingPoints.push_back(generateSampler(samplerDefinition,  //
                                                         coarseMesh));
    }

    json coarseSamplers = m_coarseSamplingPoints;

    std::ofstream oCoarse(m_coarseSamplersFilePath);
    oCoarse << std::setw(4) << coarseSamplers << std::endl;
}

std::filesystem::path SamplingPointsExporter::coarseSamplersFilePath() const
{
    return m_coarseSamplersFilePath;
}

void SamplingPointsExporter::setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath)
{
    m_coarseSamplersFilePath = coarseSamplersFilePath;
}

}  // namespace Gen202
}  // namespace Generators
}  // namespace VNCS
