#ifndef VNCS_GENERATORS_GEN202_GENERATE_SAMPLER_H
#define VNCS_GENERATORS_GEN202_GENERATE_SAMPLER_H

#include <VNCS/Spaces.h>
#include <VNCS/SamplingPoints.h>

#include <VNCS/Generators/Gen202/SamplingPointsExporter.h>

namespace VNCS
{
namespace Generators
{
class PMPTriangleIntersection;

namespace Gen202
{
struct GenerateSampler {
    GenerateSampler() = default;

    VNCS::Sampler<VNCS::Space2D> operator()(const SamplingPointsExporter::CoarseSamplerDefinition &samplerDefinition,
                                            const SamplingPointsExporter::Mesh &coarseMesh);
};
}  // namespace Gen202
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_GENERATORS_GEN212_GENERATE_SAMPLER_H
