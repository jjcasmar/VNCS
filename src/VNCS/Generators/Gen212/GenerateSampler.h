#ifndef VNCS_GENERATORS_GEN212_GENERATE_SAMPLER_H
#define VNCS_GENERATORS_GEN212_GENERATE_SAMPLER_H

#include <VNCS/Spaces.h>
#include <VNCS/SamplingPoints.h>

#include <VNCS/Generators/Gen212/SamplingPointsExporter.h>

namespace VNCS
{
namespace Generators
{
class PMPTriangleIntersection;

namespace Gen212
{
struct GenerateSampler {
    GenerateSampler() = default;

    VNCS::Sampler<VNCS::Space2D> operator()(const SamplingPointsExporter::CoarseSamplerDefinition &samplerDefinition,
                                            const SamplingPointsExporter::Mesh &coarseMesh);

    VNCS::Sampler<VNCS::Space1D> operator()(const SamplingPointsExporter::FineSamplerDefinition &samplerDefinition);
};
}  // namespace Gen212
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_GENERATORS_GEN212_GENERATE_SAMPLER_H
