#ifndef VNCS_GENERATORS_GEN202_SAMPLING_POINTS_EXPORTER_H
#define VNCS_GENERATORS_GEN202_SAMPLING_POINTS_EXPORTER_H

#include <VNCS/Spaces.h>
#include <VNCS/Generators/BlendingField.h>
#include <VNCS/SamplingPoints.h>
#include <filesystem>

namespace VNCS
{
namespace Generators
{
namespace Gen202
{
class SamplingPointsExporter
{
public:
    using K = VNCS::Space2D::K;
    using Real = VNCS::Space2D::Real;
    using Mesh = VNCS::Space2D::Mesh;
    using Point = VNCS::Space2D::Point;

    struct CoarseSamplerDefinition {
        Eigen::Vector2d coordinates;
        Real weight;
        Mesh::Face_index baseFace;
    };

    SamplingPointsExporter() = default;

    void operator()(const Mesh &coarseMesh);

    std::filesystem::path coarseSamplersFilePath() const;
    void setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath);

private:
    VNCS::SamplingPoints<VNCS::Space2D> m_coarseSamplingPoints;

    std::filesystem::path m_coarseSamplersFilePath;
};

}  // namespace Gen202
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_GENERATORS_GEN212_SAMPLING_POINTS_EXPORTER_H
