#ifndef VNCS_GENERATORS_GEN202_GEN202_H
#define VNCS_GENERATORS_GEN202_GEN202_H

#include <filesystem>
#include <memory>

#include <VNCS/Generators/BlendingField.h>
#include <VNCS/Spaces.h>
#include <VNCS/Generators/AdaptiveMeshCriteria.h>

namespace VNCS
{
namespace Generators
{
namespace Gen202
{
class Generator
{
public:
    Generator() = default;

    void operator()();

    std::filesystem::path coarseMeshPath() const;
    void setCoarseMeshPath(const std::filesystem::path &coarseMeshPath);

    std::filesystem::path coarseSamplersFilePath() const;
    void setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath);

    std::filesystem::path dofFilePath() const;
    void setDofFilePath(const std::filesystem::path &dofFilePath);

    std::shared_ptr<AdaptiveMeshCriteriaFunctor> coarseCriteria() const;
    void setCoarseCriteria(const std::shared_ptr<AdaptiveMeshCriteriaFunctor> &coarseCriteria);

    bool remeshCoarseMesh() const;
    void setRemeshCoarseMesh(bool newRemeshCoarseMesh);

private:
    std::filesystem::path m_coarseMeshPath;
    bool m_remeshCoarseMesh = true;

    std::filesystem::path m_coarseSamplersFilePath;
    std::filesystem::path m_dofFilePath;

    std::shared_ptr<AdaptiveMeshCriteriaFunctor> m_coarseCriteria;
};
}  // namespace Gen202
}  // namespace Generators
}  // namespace VNCS

#endif  //  VNCS_GENERATORS_GEN212_GEN212_H
