#ifndef VNCS_GENERATORS_GEN222_GEN222_H
#define VNCS_GENERATORS_GEN222_GEN222_H

#include <filesystem>
#include <memory>

#include <VNCS/Generators/BlendingField.h>
#include <VNCS/Spaces.h>
#include <VNCS/Generators/AdaptiveMeshCriteria.h>

namespace VNCS
{
namespace Generators
{
namespace Gen222
{
class Generator
{
public:
    Generator() = default;

    void operator()();

    void setBlendingField(std::shared_ptr<BlendingField<VNCS::Space2D::Real, VNCS::Space2D::Point>> blendingField);

    std::filesystem::path coarseSamplersFilePath() const;
    void setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath);

    std::filesystem::path fineSamplersFilePath() const;
    void setFineSamplersFilePath(const std::filesystem::path &fineSamplersFilePath);

    std::filesystem::path fineNodeSamplersFilePath() const;
    void setFineNodeSamplersFilePath(const std::filesystem::path &fineNodeSamplersFilePath);

    std::filesystem::path visualSamplersFilePath() const;
    void setVisualSamplersFilePath(const std::filesystem::path &visualSamplersFilePath);

    std::filesystem::path clusterMatrixFilePath() const;
    void setClusterMatrixFilePath(const std::filesystem::path &clusterMatrixFilePath);

    std::filesystem::path dofFilePath() const;
    void setDofFilePath(const std::filesystem::path &dofFilePath);

    void setGridSize(const std::size_t gridSize);

    std::filesystem::path coarseMeshPath() const;
    void setCoarseMeshPath(const std::filesystem::path &coarseMeshPath);

    std::filesystem::path fineMeshPath() const;
    void setFineMeshPath(const std::filesystem::path &fineMeshPath);

    std::filesystem::path visualMeshPath() const;
    void setVisualMeshPath(const std::filesystem::path &visualMeshPath);

    std::shared_ptr<AdaptiveMeshCriteriaFunctor> coarseCriteria() const;
    void setCoarseCriteria(const std::shared_ptr<AdaptiveMeshCriteriaFunctor> &coarseCriteria);

    std::shared_ptr<AdaptiveMeshCriteriaFunctor> fineCriteria() const;
    void setFineCriteria(const std::shared_ptr<AdaptiveMeshCriteriaFunctor> &fineCriteria);

    std::shared_ptr<AdaptiveMeshCriteriaFunctor> visualCriteria() const;
    void setVisualCriteria(const std::shared_ptr<AdaptiveMeshCriteriaFunctor> &visualCriteria);

    void setExportFunction(
        std::function<void(const Space2D::Mesh &, const std::vector<VNCS::Real> &, const std::string &)>
            exportFunction);

    bool allowBoundary() const;
    void setAllowBoundary(bool allowBoundary);

private:
    std::shared_ptr<VNCS::Generators::BlendingField<Space2D::Real, Space2D::Point>> m_blendingField;

    std::filesystem::path m_coarseMeshPath;
    std::filesystem::path m_fineMeshPath;
    std::filesystem::path m_visualMeshPath;

    std::filesystem::path m_coarseSamplersFilePath;
    std::filesystem::path m_fineSamplersFilePath;
    std::filesystem::path m_fineNodeSamplersFilePath;
    std::filesystem::path m_visualSamplersFilePath;
    std::filesystem::path m_clusterMatrixFilePath;
    std::filesystem::path m_dofFilePath;
    std::size_t m_gridSize;
    bool m_allowBoundary = false;

    std::shared_ptr<AdaptiveMeshCriteriaFunctor> m_coarseCriteria;
    std::shared_ptr<AdaptiveMeshCriteriaFunctor> m_fineCriteria;
    std::shared_ptr<AdaptiveMeshCriteriaFunctor> m_visualCriteria;

    std::function<void(const Space2D::Mesh &, const std::vector<VNCS::Real> &, const std::string &)> m_exportFunction;
};
}  // namespace Gen222
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_GENERATORS_GEN222_H
