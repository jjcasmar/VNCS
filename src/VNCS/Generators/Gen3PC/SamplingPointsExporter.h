#ifndef VNCS_GENERATORS_GEN3PC_SAMPLINGPOINTSEXPORTER_H
#define VNCS_GENERATORS_GEN3PC_SAMPLINGPOINTSEXPORTER_H

#include <VNCS/Spaces.h>
#include <VNCS/SamplingPoints.h>

namespace VNCS
{
namespace Generators
{
namespace Gen333
{
class SamplingPointsExporter
{
public:
    using K = VNCS::Space3D::K;
    using Real = VNCS::Space3D::Real;
    using TetraMesh = VNCS::Space3D::TetraMesh;
    using Point = VNCS::Space3D::Point;

    SamplingPointsExporter() = default;

    void operator()(const TetraMesh &coarseMesh,
                    const TetraMesh &fineMesh);

    std::filesystem::path coarseSamplersFilePath() const { return m_coarseSamplersFilePath; }
    void setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath)
    {
        m_coarseSamplersFilePath = coarseSamplersFilePath;
    }

    std::filesystem::path fineSamplersFilePath() const { return m_fineSamplersFilePath; }
    void setFineSamplersFilePath(const std::filesystem::path &fineSamplersFilePath)
    {
        m_fineSamplersFilePath = fineSamplersFilePath;
    }

private:
    VNCS::SamplingPoints<VNCS::Space3D> m_coarseSamplingPoints;
    VNCS::SamplingPoints<VNCS::Space3D> m_fineSamplingPoints;

    std::filesystem::path m_coarseSamplersFilePath;
    std::filesystem::path m_fineSamplersFilePath;
};
}  // namespace Gen333
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_GENERATORS_GEN333_SAMPLINGPOINTSEXPORTER_H
