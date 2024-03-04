#ifndef VNCS_SIM3D_STABLENEOHOOKEAN_H
#define VNCS_SIM3D_STABLENEOHOOKEAN_H

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/ForceField.inl>

#include <VNCS/Spaces.h>
#include <VNCS/SamplingPoints.h>
#include <VNCS/DeformationGradientTypes.h>

namespace VNCS
{
namespace Sim3D
{
class StableNeoHookean : public sofa::core::behavior::ForceField<VNCS::Space3D::F33>
{
    using Inherit = sofa::core::behavior::ForceField<VNCS::Space3D::F33>;
    using Real = VNCS::Space3D::Real;
    using DataVecCoord = Inherit::DataVecCoord;
    using DataVecDeriv = Inherit::DataVecDeriv;

public:
    SOFA_CLASS(StableNeoHookean, SOFA_TEMPLATE(sofa::core::behavior::ForceField, VNCS::Space3D::F33));
    StableNeoHookean() = default;

    void init() final;
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

    void setSamplingPoints(const std::shared_ptr<VNCS::SamplingPoints<VNCS::Space3D>> &samplingPoints);

private:
    VNCS::Space3D::Real psi(const VNCS::Space3D::F33::Coord &F) const;
    Eigen::Matrix3d pk1(const VNCS::Space3D::F33::Coord &F) const;
    Eigen::Matrix<VNCS::Space3D::Real, 9, 9> partialPpartialF(const VNCS::Space3D::F33::Coord &F) const;
    Eigen::Matrix<VNCS::Space3D::Real, 9, 9> clampedPartialPpartialF(const VNCS::Space3D::F33::Coord &F) const;

    VNCS::Space3D::Real m_mu;
    VNCS::Space3D::Real m_lambda;
    VNCS::Space3D::Real m_ratio;

    std::vector<Eigen::Matrix<VNCS::Space3D::Real, 9, 9>> m_hessians;

    std::shared_ptr<VNCS::SamplingPoints<VNCS::Space3D>> m_samplingPoints;
};
}  // namespace Sim3D
}  // namespace VNCS

#endif
