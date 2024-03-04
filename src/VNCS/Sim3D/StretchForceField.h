#ifndef VNCS_SIM3D_STRETCHFORCEFIELD_H
#define VNCS_SIM3D_STRETCHFORCEFIELD_H

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/ForceField.inl>

#include <VNCS/Spaces.h>
#include <VNCS/EdgeMesh.h>
#include <VNCS/SamplingPoints.h>
#include <VNCS/DeformationGradientTypes.h>

namespace VNCS
{
namespace Sim3D
{
class StretchForceField : public sofa::core::behavior::ForceField<VNCS::Space3D::VecType>
{
    using Inherit = sofa::core::behavior::ForceField<VNCS::Space3D::VecType>;
    using Real = VNCS::Space3D::Real;
    using DataVecCoord = Inherit::DataVecCoord;
    using DataVecDeriv = Inherit::DataVecDeriv;

public:
    SOFA_CLASS(StretchForceField, SOFA_TEMPLATE(sofa::core::behavior::ForceField, VNCS::Space3D::VecType));
    StretchForceField();

    void init() final;
    Real getPotentialEnergy(const sofa::core::MechanicalParams * /*mparams*/, const DataVecCoord &x) const final;
    void addForce(const sofa::core::MechanicalParams * /*mparams*/,
                  DataVecDeriv &f,
                  const DataVecCoord &x,
                  const DataVecDeriv &v) final;
    void addDForce(const sofa::core::MechanicalParams *mparams, DataVecDeriv &df, const DataVecDeriv &dx) final;

    Real stretchStiffness() const;
    void setStretchStiffness(Real newStretchStiffness);

    void setEdgeMeshPath(const std::filesystem::path &meshPath);

    const std::vector<Eigen::Matrix<VNCS::Real, 6, 6>> &hessians() const;

private:
    struct StretchElement {
        int i0;
        int i1;
        Eigen::Vector3d restLength;
    };

    Real m_stretchStiffness;
    std::vector<StretchElement> m_elements;
    std::vector<Eigen::Matrix<VNCS::Real, 6, 6>> m_hessians;

    EdgeMesh<VNCS::Space3D> m_edgeMesh;
};
}  // namespace Sim3D
}  // namespace VNCS

#endif  // VNCS_SPRINGFORCEFIELD_H
