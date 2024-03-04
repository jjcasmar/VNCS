#ifndef VNCS_GENERATORS_GEN323_GENERATESAMPLER_H
#define VNCS_GENERATORS_GEN323_GENERATESAMPLER_H

#include <VNCS/Spaces.h>
#include <VNCS/SamplingPoints.h>

#include <VNCS/Generators/Gen323/SamplingPointsExporter.h>

namespace VNCS
{
namespace Generators
{
namespace Gen323
{
class GenerateSampler
{
public:
    GenerateSampler() = default;

    VNCS::Sampler<VNCS::Space3D> operator()(const SamplingPointsExporter::CoarseSamplerDefinition &samplerDefinition,
                                            const VNCS::Space3D::TetraMesh &coarseMesh);

    VNCS::Sampler<VNCS::Space2D> operator()(const VNCS::Space3D::Mesh::Face_index &samplerDefinition,
                                            const VNCS::Space3D::Mesh &mesh);
};
}  // namespace Gen323
}  // namespace Generators
}  // namespace VNCS

#endif  //  VNCS_GENERATORS_GEN333_GENERATESAMPLER_H
