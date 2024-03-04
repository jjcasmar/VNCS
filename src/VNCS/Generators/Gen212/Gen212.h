#ifndef VNCS_GENERATORS_GEN212_GEN212_H
#define VNCS_GENERATORS_GEN212_GEN212_H

#include <filesystem>
#include <memory>

#include <VNCS/Spaces.h>
#include <VNCS/Generators/AdaptiveMeshCriteria.h>

namespace VNCS
{
namespace Generators
{
namespace Gen212
{
class Generator
{
public:
    Generator() = default;

    void operator()();

    std::filesystem::path coarseMeshPath() const;
    void setCoarseMeshPath(const std::filesystem::path &coarseMeshPath);

    std::filesystem::path edgeMeshPath() const;
    void setEdgeMeshPath(const std::filesystem::path &edgeMeshPath);

    std::filesystem::path coarseSamplersFilePath() const;
    void setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath);

    std::filesystem::path fineSamplersFilePath() const;
    void setFineSamplersFilePath(const std::filesystem::path &fineSamplersFilePath);

    std::filesystem::path clusterMatrixFilePath() const;
    void setClusterMatrixFilePath(const std::filesystem::path &clusterMatrixFilePath);

    std::filesystem::path dofFilePath() const;
    void setDofFilePath(const std::filesystem::path &dofFilePath);

    std::shared_ptr<AdaptiveMeshCriteriaFunctor> coarseCriteria() const;
    void setCoarseCriteria(const std::shared_ptr<AdaptiveMeshCriteriaFunctor> &coarseCriteria);

    bool remeshCoarseMesh() const;
    void setRemeshCoarseMesh(bool newRemeshCoarseMesh);

private:
    std::filesystem::path m_coarseMeshPath;
    std::filesystem::path m_edgeMeshPath;
    bool m_remeshCoarseMesh = true;

    std::filesystem::path m_coarseSamplersFilePath;
    std::filesystem::path m_fineSamplersFilePath;
    std::filesystem::path m_clusterMatrixFilePath;
    std::filesystem::path m_dofFilePath;

    std::shared_ptr<AdaptiveMeshCriteriaFunctor> m_coarseCriteria;
};
}  // namespace Gen212
}  // namespace Generators
}  // namespace VNCS

#endif  //  VNCS_GENERATORS_GEN212_GEN212_H
