#ifndef VNCS_GENERATORS_GEN223_SAMPLINGPOINTSEXPORTER_H
#define VNCS_GENERATORS_GEN223_SAMPLINGPOINTSEXPORTER_H

#include <VNCS/Spaces.h>
#include <VNCS/SamplingPoints.h>

namespace VNCS
{
namespace Generators
{
namespace Gen223
{
class SamplingPointsExporter
{
public:
    using K = VNCS::Space3D::K;
    using Real = VNCS::Space3D::Real;
    using Mesh = VNCS::Space3D::Mesh;
    using Point = VNCS::Space3D::Point;

    SamplingPointsExporter() = default;

    void operator()(const Mesh &coarseMesh, const Mesh &fineMesh);

    std::filesystem::path coarseSamplersFilePath() const;
    void setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath);

    std::filesystem::path fineSamplersFilePath() const;
    void setFineSamplersFilePath(const std::filesystem::path &fineSamplersFilePath);

private:
    VNCS::SamplingPoints<VNCS::Space2D> m_coarseSamplingPoints;
    VNCS::SamplingPoints<VNCS::Space2D> m_fineSamplingPoints;

    std::filesystem::path m_coarseSamplersFilePath;
    std::filesystem::path m_fineSamplersFilePath;
};
}  // namespace Gen223
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_GENERATORS_GEN333_SAMPLINGPOINTSEXPORTER_H
