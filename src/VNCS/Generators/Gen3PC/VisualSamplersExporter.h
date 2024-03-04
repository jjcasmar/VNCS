#ifndef VNCS_GENERATORS_GEN333_VISUALSAMPLERSEXPORTER_H
#define VNCS_GENERATORS_GEN333_VISUALSAMPLERSEXPORTER_H

#include <VNCS/Spaces.h>
#include <VNCS/Generators/BlendingField.h>
#include <filesystem>
#include <Eigen/Sparse>

namespace VNCS
{
namespace Generators
{
namespace Gen333
{
class VisualSamplersExporter
{
    using K = VNCS::Space3D::K;
    using Real = VNCS::Space3D::Real;
    using Mesh = VNCS::Space3D::Mesh;

public:
    VisualSamplersExporter() = default;

    std::pair<Eigen::SparseMatrix<VNCS::Real>, Eigen::SparseMatrix<VNCS::Real>> matrices() const;

    void setVisualMeshPath(const std::filesystem::path &visualPath);
    void setFineTetraMesh(const VNCS::Space3D::TetraMesh &mesh);
    void setCoarseTetraMesh(const VNCS::Space3D::TetraMesh &mesh);

    void setBlendingField(std::shared_ptr<BlendingField<VNCS::Space3D::Real, VNCS::Space3D::Point>> blendingField);

private:
    VNCS::Space3D::Mesh m_visualMesh;
    VNCS::Space3D::TetraMesh m_coarseMesh;
    VNCS::Space3D::TetraMesh m_fineMesh;

    std::shared_ptr<VNCS::Generators::BlendingField<Space3D::Real, Space3D::Point>> m_blendingField;
};
}  // namespace Gen333
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_GENERATORS_GEN333_VISUALSAMPLERSEXPORTER_H
