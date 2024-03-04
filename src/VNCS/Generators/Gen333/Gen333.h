#ifndef VNCS_GENERATORS_GEN333_GENERATOR_H
#define VNCS_GENERATORS_GEN333_GENERATOR_H

#include <filesystem>
#include <memory>

#include <VNCS/Generators/BlendingField.h>
#include <VNCS/Spaces.h>
#include <VNCS/Generators/TetraMeshGenerator.h>

namespace VNCS
{
namespace Generators
{
namespace Gen333
{
class Generator
{
public:
    Generator() = default;

    void operator()();

    void setBlendingField(std::shared_ptr<BlendingField<VNCS::Space3D::Real, VNCS::Space3D::Point>> blendingField);

    std::filesystem::path coarseSamplersFilePath() const;
    void setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath);

    std::filesystem::path fineSamplersFilePath() const;
    void setFineSamplersFilePath(const std::filesystem::path &fineSamplersFilePath);

    std::filesystem::path visualSamplersFilePath() const;
    void setVisualSamplersFilePath(const std::filesystem::path &visualSamplersFilePath);

    std::filesystem::path clusterMatrixFilePath() const;
    void setClusterMatrixFilePath(const std::filesystem::path &clusterMatrixFilePath);

    std::filesystem::path dofFilePath() const;
    void setDofFilePath(const std::filesystem::path &dofFilePath);

    std::filesystem::path coarseMeshPath() const;
    void setCoarseMeshPath(const std::filesystem::path &coarseMeshPath);

    std::filesystem::path fineMeshPath() const;
    void setFineMeshPath(const std::filesystem::path &fineMeshPath);

    std::filesystem::path visualMeshPath() const;
    void setVisualMeshPath(const std::filesystem::path &visualMeshPath);

    VNCS::Generators::C3t3Criteria coarseCriteria() const;
    void setCoarseCriteria(const VNCS::Generators::C3t3Criteria &coarseCriteria);

    VNCS::Generators::C3t3Criteria fineCriteria() const;
    void setFineCriteria(const VNCS::Generators::C3t3Criteria &fineCriteria);

    VNCS::Generators::RemeshCriteria visualCriteria() const;
    void setVisualCriteria(const VNCS::Generators::RemeshCriteria &visualCriteria);

    void setExportFunction(
        std::function<void(const Space3D::TetraMesh &,
                           const std::string &,
                           const std::unordered_map<std::string, std::vector<VNCS::Space3D::Real>> &)> exportFunction);

private:
    std::shared_ptr<VNCS::Generators::BlendingField<Space3D::Real, Space3D::Point>> m_blendingField;

    std::filesystem::path m_coarseMeshPath;
    std::filesystem::path m_fineMeshPath;
    std::filesystem::path m_visualMeshPath;

    std::filesystem::path m_coarseSamplersFilePath;
    std::filesystem::path m_fineSamplersFilePath;
    std::filesystem::path m_visualSamplersFilePath;
    std::filesystem::path m_clusterMatrixFilePath;
    std::filesystem::path m_dofFilePath;

    VNCS::Generators::C3t3Criteria m_coarseCriteria;
    VNCS::Generators::C3t3Criteria m_fineCriteria;
    VNCS::Generators::RemeshCriteria m_visualCriteria;

    std::function<void(const Space3D::TetraMesh &,
                       const std::string &,
                       const std::unordered_map<std::string, std::vector<VNCS::Space3D::Real>> &)>
        m_exportFunction;
};

}  // namespace Gen333
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_GENERATORS_GEN333_GENERATOR_H
