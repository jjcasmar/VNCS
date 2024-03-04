#ifndef VNCS_GENERATORS_GEN323_GENERATOR_H
#define VNCS_GENERATORS_GEN323_GENERATOR_H

#include <filesystem>
#include <memory>

#include <VNCS/Spaces.h>
#include <VNCS/Generators/AdaptiveMeshCriteria.h>
#include <VNCS/Generators/TetraMeshGenerator.h>

namespace VNCS
{
namespace Generators
{
namespace Gen323
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

    VNCS::Generators::C3t3Criteria coarseCriteria() const;
    void setCoarseCriteria(const VNCS::Generators::C3t3Criteria &coarseCriteria);

    VNCS::Generators::RemeshCriteria fineCriteria() const;
    void setFineCriteria(const VNCS::Generators::RemeshCriteria &fineCriteria);

    void setExportFunction(
        std::function<void(const Space3D::TetraMesh &,
                           const std::string &,
                           const std::unordered_map<std::string, std::vector<VNCS::Space3D::Real>> &)> exportFunction);

    bool isBoundaryAllowed() const;
    void alloBoundary(bool allowBoundary);

private:
    std::filesystem::path m_coarseMeshPath;
    std::filesystem::path m_fineMeshPath;

    std::filesystem::path m_coarseSamplersFilePath;
    std::filesystem::path m_fineSamplersFilePath;
    std::filesystem::path m_clusterMatrixFilePath;
    std::filesystem::path m_dofFilePath;

    VNCS::Generators::C3t3Criteria m_coarseCriteria;
    VNCS::Generators::RemeshCriteria m_fineCriteria;

    bool m_allowBoundary;

    std::function<void(const Space3D::TetraMesh &,
                       const std::string &,
                       const std::unordered_map<std::string, std::vector<VNCS::Space3D::Real>> &)>
        m_exportFunction;
};

}  // namespace Gen323
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_SIMCREATOR_H
