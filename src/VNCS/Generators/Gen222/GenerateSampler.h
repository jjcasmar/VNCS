#ifndef VNCS_GENERATORS_GEN222_GENERATESAMPLER_H
#define VNCS_GENERATORS_GEN222_GENERATESAMPLER_H

#include <VNCS/Spaces.h>
#include <VNCS/Generators/PMPTriangleIntersection.h>

#include <VNCS/Generators/Gen222/SamplingPointsExporter.h>

namespace VNCS
{
namespace Generators
{
namespace Gen222
{
struct GenerateSampler {
    const VNCS::Generators::PMPTriangleIntersection &m_coarseIntersectionFunc;
    const VNCS::Generators::PMPTriangleIntersection &m_fineIntersectionFunc;

    struct CoarseSampler {
    };
    struct FineSampler {
    };

    GenerateSampler(const VNCS::Generators::PMPTriangleIntersection &coarseIntersectionFunc,
                    const VNCS::Generators::PMPTriangleIntersection &fineIntersectionFunc);

    VNCS::Sampler<VNCS::Space2D> operator()(const SamplingPointsExporter::SamplerDefinition &samplerDefinition,
                                            const SamplingPointsExporter::Mesh &coarseMesh,
                                            CoarseSampler);

    VNCS::Sampler<VNCS::Space2D> operator()(const SamplingPointsExporter::SamplerDefinition &samplerDefinition,
                                            const SamplingPointsExporter::Mesh &mesh,
                                            FineSampler);
};
}  // namespace Gen222
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_GENERATORS_GEN222_GENERATESAMPLER_H
