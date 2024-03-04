#ifndef VNCS_GENERATORS_GEN323_SAMPLINGPOINTSEXPORTER_H
#define VNCS_GENERATORS_GEN323_SAMPLINGPOINTSEXPORTER_H

#include <VNCS/Spaces.h>
#include <VNCS/SamplingPoints.h>

namespace VNCS
{
namespace Generators
{
namespace Gen323
{
class SamplingPointsExporter
{
public:
    using K = VNCS::Space3D::K;
    using Real = VNCS::Space3D::Real;
    using Mesh = VNCS::Space3D::Mesh;
    using TetraMesh = VNCS::Space3D::TetraMesh;
    using Point = VNCS::Space3D::Point;

    struct CoarseSamplerDefinition {
        Eigen::Vector3d coordinates;
        Real weight;
        int baseFace;
    };

    struct FineSamplerDefinition {
        Real weight;
        VNCS::Space2D::Mesh::Face_index faceHandle;
    };

    SamplingPointsExporter() = default;

    void operator()(const TetraMesh &coarseMesh, const Mesh &fineMesh);

    std::filesystem::path coarseSamplersFilePath() const;
    void setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath);

    std::filesystem::path fineSamplersFilePath() const;
    void setFineSamplersFilePath(const std::filesystem::path &fineSamplersFilePath);

private:
    VNCS::SamplingPoints<VNCS::Space3D> m_coarseSamplingPoints;
    VNCS::SamplingPoints<VNCS::Space2D> m_fineSamplingPoints;

    std::filesystem::path m_coarseSamplersFilePath;
    std::filesystem::path m_fineSamplersFilePath;
};
}  // namespace Gen323
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_GENERATORS_GEN333_SAMPLINGPOINTSEXPORTER_H
