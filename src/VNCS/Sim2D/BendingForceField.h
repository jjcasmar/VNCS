#ifndef VNCS_SIM2D_BENDINGFORCEFIELD_H
#define VNCS_SIM2D_BENDINGFORCEFIELD_H

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/ForceField.inl>

#include <VNCS/Spaces.h>
#include <VNCS/EdgeMesh.h>
#include <VNCS/SamplingPoints.h>
#include <VNCS/DeformationGradientTypes.h>

namespace VNCS
{
namespace Sim2D
{
class BendingForceField : public sofa::core::behavior::ForceField<VNCS::Space2D::VecType>
{
    using Inherit = sofa::core::behavior::ForceField<VNCS::Space2D::VecType>;
    using Real = VNCS::Space2D::Real;
    using DataVecCoord = Inherit::DataVecCoord;
    using DataVecDeriv = Inherit::DataVecDeriv;

public:
    SOFA_CLASS(BendingForceField, SOFA_TEMPLATE(sofa::core::behavior::ForceField, VNCS::Space2D::VecType));
    BendingForceField();

    void init() final;
    Real getPotentialEnergy(const sofa::core::MechanicalParams * /*mparams*/, const DataVecCoord &x) const final;
    void addForce(const sofa::core::MechanicalParams * /*mparams*/,
                  DataVecDeriv &f,
                  const DataVecCoord &x,
                  const DataVecDeriv &v) final;
    void addDForce(const sofa::core::MechanicalParams *mparams, DataVecDeriv &df, const DataVecDeriv &dx) final;

    Real bendingStiffness() const;
    void setBendingStiffness(Real newBendingStiffness);

    void setEdgeMeshPath(const std::filesystem::path &meshPath);

    const std::vector<Eigen::Matrix<VNCS::Real, 6, 6>> &hessians() const;

private:
    struct BendingElement {
        int indexLeft;
        int indexCenter;
        int indexRight;
        VNCS::Real restLength;
    };

    Real m_bendingStiffness;
    std::vector<BendingElement> m_elements;
    std::vector<Eigen::Matrix<VNCS::Real, 6, 6>> m_hessians;

    VNCS::EdgeMesh<VNCS::Space2D> m_edgeMesh;
};
}  // namespace Sim2D
}  // namespace VNCS

#endif  // VNCS_SPRINGFORCEFIELD_H
