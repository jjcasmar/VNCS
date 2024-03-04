#ifndef VNCS_SIM2D_STVKFORCEFIELD_H
#define VNCS_SIM2D_STVKFORCEFIELD_H

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/ForceField.inl>

#include <VNCS/Spaces.h>
#include <VNCS/SamplingPoints.h>
#include <VNCS/DeformationGradientTypes.h>

namespace VNCS
{
namespace Sim2D
{
class StVKForceField : public sofa::core::behavior::ForceField<VNCS::Space2D::F22>
{
    using Inherit = sofa::core::behavior::ForceField<VNCS::Space2D::F22>;
    using Real = VNCS::Space2D::Real;
    using DataVecCoord = Inherit::DataVecCoord;
    using DataVecDeriv = Inherit::DataVecDeriv;

public:
    SOFA_CLASS(StVKForceField, SOFA_TEMPLATE(sofa::core::behavior::ForceField, VNCS::Space2D::F22));
    StVKForceField();

    void init();
    Real getPotentialEnergy(const sofa::core::MechanicalParams * /*mparams*/, const DataVecCoord &x) const final;
    void addForce(const sofa::core::MechanicalParams * /*mparams*/,
                  DataVecDeriv &f,
                  const DataVecCoord &x,
                  const DataVecDeriv &v) final;
    void addDForce(const sofa::core::MechanicalParams *mparams, DataVecDeriv &df, const DataVecDeriv &dx) final;

    Real young() const;
    void setYoung(const Real &young);

    Real beta() const;
    void setBeta(const Real &beta);

    Real poisson() const;
    void setPoisson(const Real &poisson);

    std::shared_ptr<VNCS::SamplingPoints<VNCS::Space2D>> samplingPoints() const;
    void setSamplingPoints(const std::shared_ptr<VNCS::SamplingPoints<VNCS::Space2D>> &samplingPoints);

private:
    Real m_youngModulus;
    Real m_poissonRatio;
    Real muLame;
    Real lambdaLame;
    Real mu;
    Real lambda;
    Real alpha;
    std::optional<Real> m_beta;

    using Hessian = Eigen::Matrix<VNCS::Real, 4, 4>;
    std::vector<std::pair<Hessian, Hessian>> m_dissipationHessians;

    std::shared_ptr<VNCS::SamplingPoints<VNCS::Space2D>> m_samplingPoints;
};
}  // namespace Sim2D
}  // namespace VNCS

#endif  // VNCS_STVKFORCEFIELD_H
