#ifndef VNCS_GENERATORS_GEN3PC_GENERATOR_H
#define VNCS_GENERATORS_GEN3PC_GENERATOR_H

#include <filesystem>
#include <memory>

#include <VNCS/Generators/BlendingField.h>
#include <VNCS/Spaces.h>
#include <VNCS/Generators/TetraMeshGenerator.h>

namespace VNCS
{
namespace Generators
{
namespace Gen3PC
{
class Generator
{
public:
    Generator() = default;

    void operator()();

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

    VNCS::Generators::C3t3Criteria fineCriteria() const;
    void setFineCriteria(const VNCS::Generators::C3t3Criteria &fineCriteria);

    void setExportFunction(
        std::function<void(const Space3D::TetraMesh &,
                           const std::string &)>);

private:
    std::shared_ptr<VNCS::Generators::BlendingField<Space3D::Real, Space3D::Point>> m_blendingField;

    std::filesystem::path m_coarseMeshPath;
    std::filesystem::path m_fineMeshPath;

    std::filesystem::path m_fineSamplersFilePath;
    std::filesystem::path m_clusterMatrixFilePath;
    std::filesystem::path m_dofFilePath;

    VNCS::Generators::C3t3Criteria m_fineCriteria;

    std::function<void(const Space3D::TetraMesh &,
                       const std::string &)>
        m_exportFunction;
};

}  // namespace Gen3PC
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_GENERATORS_GEN333_GENERATOR_H
