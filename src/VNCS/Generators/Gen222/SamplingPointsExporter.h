#ifndef VNCS_GENERATORS_GEN222_SAMPLING_POINTS_EXPORTER_H
#define VNCS_GENERATORS_GEN222_SAMPLING_POINTS_EXPORTER_H

#include <VNCS/Spaces.h>
#include <VNCS/Generators/BlendingField.h>
#include <VNCS/SamplingPoints.h>
#include <filesystem>

namespace VNCS
{
namespace Generators
{
namespace Gen222
{
class SamplingPointsExporter
{
public:
    using K = VNCS::Space2D::K;
    using Real = VNCS::Space2D::Real;
    using Mesh = VNCS::Space2D::Mesh;
    using Point = VNCS::Space2D::Point;

    struct SamplerDefinition {
        enum LocationType {
            CoarseTriangle,  //
            FineTriangle
        };

        enum KinematicType {
            Coarse,
            Hybrid,
            Fine,
        };

        Eigen::Vector2d coordinates;
        Real weight;
        LocationType locationType;
        Mesh::Face_index baseFace;
        std::optional<Mesh::Face_index> supportFace;
        KinematicType kinematicType = KinematicType::Hybrid;
    };

    SamplingPointsExporter() = default;

    void operator()(const Mesh &coarseMesh, const Mesh &fineMesh, std::size_t gridSize);

    std::filesystem::path coarseSamplersFilePath() const;
    void setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath);

    std::filesystem::path fineSamplersFilePath() const;
    void setFineSamplersFilePath(const std::filesystem::path &fineSamplersFilePath);

    std::filesystem::path fineNodeSamplersFilePath() const;
    void setFineNodeSamplersFilePath(const std::filesystem::path &fineNodeSamplersFilePath);

private:
    VNCS::SamplingPoints<VNCS::Space2D> m_coarseSamplingPoints;
    VNCS::SamplingPoints<VNCS::Space2D> m_fineSamplingPoints;

    std::filesystem::path m_coarseSamplersFilePath;
    std::filesystem::path m_fineNodeSamplersFilePath;
    std::filesystem::path m_fineSamplersFilePath;
};

}  // namespace Gen222
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_GENERATORS_GEN222_SAMPLING_POINTS_EXPORTER_H
