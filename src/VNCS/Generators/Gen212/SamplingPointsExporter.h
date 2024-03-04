#ifndef VNCS_GENERATORS_GEN212_SAMPLING_POINTS_EXPORTER_H
#define VNCS_GENERATORS_GEN212_SAMPLING_POINTS_EXPORTER_H

#include <VNCS/Spaces.h>
#include <VNCS/SamplingPoints.h>
#include <VNCS/EdgeMesh.h>
#include <filesystem>

namespace VNCS
{
namespace Generators
{
namespace Gen212
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

    struct FineSamplerDefinition {
        Real weight;
        EdgeMesh<VNCS::Space2D>::vertex_descriptor pointHandle;
    };

    SamplingPointsExporter() = default;

    void operator()(const Mesh &coarseMesh, const VNCS::EdgeMesh<VNCS::Space2D> &fineMesh);

    std::filesystem::path coarseSamplersFilePath() const;
    void setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath);

    std::filesystem::path fineSamplersFilePath() const;
    void setFineSamplersFilePath(const std::filesystem::path &fineSamplersFilePath);

private:
    VNCS::SamplingPoints<VNCS::Space2D> m_coarseSamplingPoints;
    VNCS::SamplingPoints<VNCS::Space1D> m_fineSamplingPoints;

    std::filesystem::path m_coarseSamplersFilePath;
    std::filesystem::path m_fineSamplersFilePath;
};

}  // namespace Gen212
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_GENERATORS_GEN212_SAMPLING_POINTS_EXPORTER_H
