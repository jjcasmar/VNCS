#ifndef VNCS_SIM3D_STVKFORCEFIELD_H
#define VNCS_SIM3D_STVKFORCEFIELD_H

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/ForceField.inl>

#include <VNCS/Spaces.h>
#include <VNCS/SamplingPoints.h>
#include <VNCS/DeformationGradientTypes.h>

namespace VNCS
{
namespace Sim3D
{
class StVKForceField : public sofa::core::behavior::ForceField<VNCS::Space3D::F32>
{
    using Inherit = sofa::core::behavior::ForceField<VNCS::Space3D::F32>;
    using Real = VNCS::Space3D::Real;
    using DataVecCoord = Inherit::DataVecCoord;
    using DataVecDeriv = Inherit::DataVecDeriv;

public:
    SOFA_CLASS(StVKForceField, SOFA_TEMPLATE(sofa::core::behavior::ForceField, VNCS::Space3D::F32));
    StVKForceField();

    void init();
    Real getPotentialEnergy(const sofa::core::MechanicalParams * /*mparams*/, const DataVecCoord &x) const final;
    void addForce(const sofa::core::MechanicalParams * /*mparams*/,
                  DataVecDeriv &f,
                  const DataVecCoord &x,
                  const DataVecDeriv &v) final;
    void addDForce(const sofa::core::MechanicalParams *mparams, DataVecDeriv &df, const DataVecDeriv &dx) final;

    VNCS::Space3D::Real mu() const;
    void setMu(const VNCS::Space3D::Real &mu);

    VNCS::Space3D::Real lambda() const;
    void setLambda(const VNCS::Space3D::Real &lambda);

    std::shared_ptr<VNCS::SamplingPoints<VNCS::Space2D>> samplingPoints() const;
    void setSamplingPoints(const std::shared_ptr<VNCS::SamplingPoints<VNCS::Space2D>> &samplingPoints);

    Real beta() const;
    void setBeta(const Real &beta);

private:
    VNCS::Space3D::Real m_mu;
    VNCS::Space3D::Real m_lambda;
    VNCS::Space3D::Real m_ratio;

    std::optional<VNCS::Real> m_beta;

    using Hessian = Eigen::Matrix<VNCS::Real, 6, 6>;

    std::shared_ptr<VNCS::SamplingPoints<VNCS::Space2D>> m_samplingPoints;
    std::vector<Eigen::Matrix<VNCS::Space3D::Real, 6, 6>> m_hessians;

    std::vector<std::pair<Hessian, Hessian>> m_dissipationHessians;
};
}  // namespace Sim3D
}  // namespace VNCS

#endif  // VNCS_STVKFORCEFIELD_H
