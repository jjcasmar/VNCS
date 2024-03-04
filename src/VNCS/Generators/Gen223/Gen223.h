#ifndef VNCS_GENERATORS_GEN223_GENERATOR_H
#define VNCS_GENERATORS_GEN223_GENERATOR_H

#include <filesystem>
#include <memory>

#include <VNCS/Spaces.h>
#include <VNCS/Generators/AdaptiveMeshCriteria.h>
#include <VNCS/Generators/TetraMeshGenerator.h>

namespace VNCS
{
namespace Generators
{
namespace Gen223
{
class Generator
{
public:
    Generator() = default;

    void operator()();

    std::filesystem::path coarseSamplersFilePath() const;
    void setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath);

    std::filesystem::path fineSamplersFilePath() const;
    void setFineSamplersFilePath(const std::filesystem::path &fineSamplersFilePath);

    std::filesystem::path clusterMatrixFilePath() const;
    void setClusterMatrixFilePath(const std::filesystem::path &clusterMatrixFilePath);

    std::filesystem::path dofFilePath() const;
    void setDofFilePath(const std::filesystem::path &dofFilePath);

    std::filesystem::path coarseMeshPath() const;
    void setCoarseMeshPath(const std::filesystem::path &coarseMeshPath);

    std::filesystem::path fineMeshPath() const;
    void setFineMeshPath(const std::filesystem::path &fineMeshPath);

    std::optional<VNCS::Generators::RemeshCriteria> coarseCriteria() const;
    void setCoarseCriteria(const VNCS::Generators::RemeshCriteria &coarseCriteria);

    std::optional<VNCS::Generators::RemeshCriteria> fineCriteria() const;
    void setFineCriteria(const VNCS::Generators::RemeshCriteria &fineCriteria);

private:
    std::filesystem::path m_coarseMeshPath;
    std::filesystem::path m_fineMeshPath;

    std::filesystem::path m_coarseSamplersFilePath;
    std::filesystem::path m_fineSamplersFilePath;
    std::filesystem::path m_clusterMatrixFilePath;
    std::filesystem::path m_dofFilePath;

    std::optional<VNCS::Generators::RemeshCriteria> m_coarseCriteria;
    std::optional<VNCS::Generators::RemeshCriteria> m_fineCriteria;
};

}  // namespace Gen223
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_SIMCREATOR_H
