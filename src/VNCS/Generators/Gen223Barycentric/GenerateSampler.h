#ifndef VNCS_GENERATORS_GEN223BARYCENTRIC_GENERATESAMPLER_H
#define VNCS_GENERATORS_GEN223BARYCENTRIC_GENERATESAMPLER_H

#include <VNCS/Spaces.h>
#include <VNCS/SamplingPoints.h>

#include <VNCS/Generators/Gen223Barycentric/SamplingPointsExporter.h>

namespace VNCS
{
namespace Generators
{
namespace Gen223Barycentric
{
class GenerateSampler
{
public:
    GenerateSampler() = default;

    VNCS::Sampler<VNCS::Space2D> operator()(const VNCS::Space3D::Mesh::Face_index &samplerDefinition,
                                            const VNCS::Space3D::Mesh &mesh);
};
}  // namespace Gen223Barycentric
}  // namespace Generators
}  // namespace VNCS

#endif  //  VNCS_GENERATORS_GEN333_GENERATESAMPLER_H
