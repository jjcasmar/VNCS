#ifndef VNCS_GENERATORS_GEN222_VISUALSAMPLERSEXPORTER_H
#define VNCS_GENERATORS_GEN222_VISUALSAMPLERSEXPORTER_H

#include <VNCS/Spaces.h>
#include <filesystem>

namespace VNCS
{
namespace Generators
{
namespace Gen222
{
class VisualSamplersExporter
{
    using K = VNCS::Space2D::K;
    using Real = VNCS::Space2D::Real;
    using Mesh = VNCS::Space2D::Mesh;

public:
    VisualSamplersExporter() = default;

    void operator()(const Mesh &visualMesh, const Mesh &coarseMesh, const Mesh &fineMesh);

    const std::filesystem::path &samplersFilePath() const;
    void setSamplersFilePath(const std::filesystem::path &samplersFilePath);

private:
    std::filesystem::path m_samplersFilePath;
};

}  // namespace Gen222
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_GENERATORS_GEN222_VISUALSAMPLERSEXPORTER_H
