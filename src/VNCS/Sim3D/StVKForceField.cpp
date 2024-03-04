#include "StVKForceField.h"

#include <VNCS/DataExtensions.h>
#include <range/v3/view/zip.hpp>

#include <spdlog/fmt/ostr.h>

VNCS::Sim3D::StVKForceField::StVKForceField()
    : Inherit()
{
}

void VNCS::Sim3D::StVKForceField::init()
{
    assert(m_mu > 0.0);
    assert(m_lambda > 0.0);
    assert(m_ratio > 0.0);

    m_ratio = m_mu / m_lambda;

    Inherit::init();
}

VNCS::Space3D::Real VNCS::Sim3D::StVKForceField::getPotentialEnergy(const sofa::core::MechanicalParams *mparams,
                                                                    const DataVecCoord &xVec) const
{
    auto energy = Real(0.0);
    auto dissipation = Real(0.0);

    const auto lame1 = m_lambda;
    const auto lame2 = m_mu;

    const auto &vVec = getMState()->readVelocities();

    const auto &samplers = *m_samplingPoints;
    for (const auto &[F, dFdt, sampler] :
         ranges::views::zip(make_read_accessor(xVec), make_read_accessor(vVec), samplers)) {
        const auto w = sampler.w;

        std::array<VNCS::Real, 6> Fv;
        Eigen::Map<Eigen::Matrix<VNCS::Real, 3, 2, Eigen::RowMajor>> fvMap(&Fv[0]);
        fvMap = F;

#include "FEMMem_StVK_Energy.mcg"

        energy += w * t32;

        if (m_beta) {
            const auto beta = m_beta.value();
            const Eigen::Matrix2d strainRate = 1.0 / 2.0 * (F.transpose() * dFdt + dFdt.transpose() * F);
            const auto e00 = strainRate(0, 0);
            const auto e11 = strainRate(1, 1);
            const auto e10 = strainRate(1, 0);

            dissipation += 0.5 * w * beta * (e00 * e00 + e11 * e11 + e10 * e10);
        }
    }

    spdlog::info("StVK energy = {}", energy);

    return energy + mparams->dt() * dissipation;
}

void VNCS::Sim3D::StVKForceField::addForce(const sofa::core::MechanicalParams *mparams,
                                           DataVecDeriv &fVec,
                                           const DataVecCoord &xVec,
                                           const DataVecDeriv &vVec)
{
    const auto lame1 = m_lambda;
    const auto lame2 = m_mu;
    const auto &samplers = *m_samplingPoints;
    const auto &Fs = make_read_accessor(xVec);
    const auto &dFdts = make_read_accessor(vVec);
    auto &fs = make_write_accessor(fVec);
    m_hessians.resize(Fs.size());
    for (int i = 0; i < Fs.size(); ++i) {
        const auto &F = Fs[i];
        const auto &dFdt = dFdts[i];
        auto &f = fs[i];
        const auto &sampler = samplers[i];
        const auto w = sampler.w;

        std::array<VNCS::Real, 6> Fv;
        Eigen::Map<Eigen::Matrix<VNCS::Real, 3, 2, Eigen::RowMajor>> fvMap(&Fv[0]);
        fvMap = F;

        {
            std::array<VNCS::Real, 6> vgx;
#include "FEMMem_StVK_Gradient.mcg"

            Eigen::Map<Eigen::Matrix<VNCS::Real, 3, 2, Eigen::RowMajor>> vgxMap(vgx.data());

            f -= w * vgxMap;
        }

        {
            std::array<std::array<VNCS::Real, 6>, 6> mHx;
#include "FEMMem_StVK_Hessian.mcg"

            Eigen::Map<Eigen::Matrix<VNCS::Real, 6, 6>> hessianMap(mHx[0].data());
            m_hessians[i] = w * hessianMap;
        }

        if (m_beta) {
            const auto beta = m_beta.value();
            const auto F00 = F(0, 0);
            const auto F01 = F(0, 1);
            const auto F10 = F(1, 0);
            const auto F11 = F(1, 1);
            const auto F20 = F(2, 0);
            const auto F21 = F(2, 1);

            const auto dFdt00 = dFdt(0, 0);
            const auto dFdt01 = dFdt(0, 1);
            const auto dFdt10 = dFdt(1, 0);
            const auto dFdt11 = dFdt(1, 1);
            const auto dFdt20 = dFdt(2, 0);
            const auto dFdt21 = dFdt(2, 1);

            // Extracted from Matlab
            const auto rStress = Eigen::Matrix<VNCS::Real, 3, 2>(
                {{2 * F00 * (F00 * dFdt00 + F10 * dFdt10 + F20 * dFdt20) +
                      F01 * ((F00 * dFdt01) / 2 + (F01 * dFdt00) / 2 + (F10 * dFdt11) / 2 + (F11 * dFdt10) / 2 +
                             (F20 * dFdt21) / 2 + (F21 * dFdt20) / 2),
                  2 * F01 * (F01 * dFdt01 + F11 * dFdt11 + F21 * dFdt21) +
                      F00 * ((F00 * dFdt01) / 2 + (F01 * dFdt00) / 2 + (F10 * dFdt11) / 2 + (F11 * dFdt10) / 2 +
                             (F20 * dFdt21) / 2 + (F21 * dFdt20) / 2)},
                 {2 * F10 * (F00 * dFdt00 + F10 * dFdt10 + F20 * dFdt20) +
                      F11 * ((F00 * dFdt01) / 2 + (F01 * dFdt00) / 2 + (F10 * dFdt11) / 2 + (F11 * dFdt10) / 2 +
                             (F20 * dFdt21) / 2 + (F21 * dFdt20) / 2),
                  2 * F11 * (F01 * dFdt01 + F11 * dFdt11 + F21 * dFdt21) +
                      F10 * ((F00 * dFdt01) / 2 + (F01 * dFdt00) / 2 + (F10 * dFdt11) / 2 + (F11 * dFdt10) / 2 +
                             (F20 * dFdt21) / 2 + (F21 * dFdt20) / 2)},
                 {2 * F20 * (F00 * dFdt00 + F10 * dFdt10 + F20 * dFdt20) +
                      F21 * ((F00 * dFdt01) / 2 + (F01 * dFdt00) / 2 + (F10 * dFdt11) / 2 + (F11 * dFdt10) / 2 +
                             (F20 * dFdt21) / 2 + (F21 * dFdt20) / 2),
                  2 * F21 * (F01 * dFdt01 + F11 * dFdt11 + F21 * dFdt21) +
                      F20 * ((F00 * dFdt01) / 2 + (F01 * dFdt00) / 2 + (F10 * dFdt11) / 2 + (F11 * dFdt10) / 2 +
                             (F20 * dFdt21) / 2 + (F21 * dFdt20) / 2)}});

            f -= w * beta * rStress;

            // Hessians
            const Hessian dE2_X = w * beta *
                                  Hessian({
                                      {4 * F00 * dFdt00 + (F01 * dFdt01) / 2 + 2 * F10 * dFdt10 + 2 * F20 * dFdt20,
                                       (F00 * dFdt01) / 2 + F01 * dFdt00 + (F10 * dFdt11) / 2 + (F11 * dFdt10) / 2 +
                                           (F20 * dFdt21) / 2 + (F21 * dFdt20) / 2,
                                       2 * F00 * dFdt10 + (F01 * dFdt11) / 2,
                                       (F01 * dFdt10) / 2,
                                       2 * F00 * dFdt20 + (F01 * dFdt21) / 2,
                                       (F01 * dFdt20) / 2},
                                      {F00 * dFdt01 + (F01 * dFdt00) / 2 + (F10 * dFdt11) / 2 + (F11 * dFdt10) / 2 +
                                           (F20 * dFdt21) / 2 + (F21 * dFdt20) / 2,
                                       (F00 * dFdt00) / 2 + 4 * F01 * dFdt01 + 2 * F11 * dFdt11 + 2 * F21 * dFdt21,
                                       (F00 * dFdt11) / 2,
                                       (F00 * dFdt10) / 2 + 2 * F01 * dFdt11,
                                       (F00 * dFdt21) / 2,
                                       (F00 * dFdt20) / 2 + 2 * F01 * dFdt21},
                                      {2 * F10 * dFdt00 + (F11 * dFdt01) / 2,
                                       (F11 * dFdt00) / 2,
                                       2 * F00 * dFdt00 + 4 * F10 * dFdt10 + (F11 * dFdt11) / 2 + 2 * F20 * dFdt20,
                                       (F00 * dFdt01) / 2 + (F01 * dFdt00) / 2 + (F10 * dFdt11) / 2 + F11 * dFdt10 +
                                           (F20 * dFdt21) / 2 + (F21 * dFdt20) / 2,
                                       2 * F10 * dFdt20 + (F11 * dFdt21) / 2,
                                       (F11 * dFdt20) / 2},
                                      {(F10 * dFdt01) / 2,
                                       (F10 * dFdt00) / 2 + 2 * F11 * dFdt01,
                                       (F00 * dFdt01) / 2 + (F01 * dFdt00) / 2 + F10 * dFdt11 + (F11 * dFdt10) / 2 +
                                           (F20 * dFdt21) / 2 + (F21 * dFdt20) / 2,
                                       2 * F01 * dFdt01 + (F10 * dFdt10) / 2 + 4 * F11 * dFdt11 + 2 * F21 * dFdt21,
                                       (F10 * dFdt21) / 2,
                                       (F10 * dFdt20) / 2 + 2 * F11 * dFdt21},
                                      {2 * F20 * dFdt00 + (F21 * dFdt01) / 2,
                                       (F21 * dFdt00) / 2,
                                       2 * F20 * dFdt10 + (F21 * dFdt11) / 2,
                                       (F21 * dFdt10) / 2,
                                       2 * F00 * dFdt00 + 2 * F10 * dFdt10 + 4 * F20 * dFdt20 + (F21 * dFdt21) / 2,
                                       (F00 * dFdt01) / 2 + (F01 * dFdt00) / 2 + (F10 * dFdt11) / 2 +
                                           (F11 * dFdt10) / 2 + (F20 * dFdt21) / 2 + F21 * dFdt20},
                                      {(F20 * dFdt01) / 2,
                                       (F20 * dFdt00) / 2 + 2 * F21 * dFdt01,
                                       (F20 * dFdt11) / 2,
                                       (F20 * dFdt10) / 2 + 2 * F21 * dFdt11,
                                       (F00 * dFdt01) / 2 + (F01 * dFdt00) / 2 + (F10 * dFdt11) / 2 +
                                           (F11 * dFdt10) / 2 + F20 * dFdt21 + (F21 * dFdt20) / 2,
                                       2 * F01 * dFdt01 + 2 * F11 * dFdt11 + (F20 * dFdt20) / 2 + 4 * F21 * dFdt21},
                                  });

            const Hessian dE2_V = w * beta *
                                  Hessian({
                                      {2 * std::pow(F00, 2.0) + std::pow(F01, 2.0) / 2,
                                       (F00 * F01) / 2,
                                       2 * F00 * F10 + (F01 * F11) / 2,
                                       (F01 * F10) / 2,
                                       2 * F00 * F20 + (F01 * F21) / 2,
                                       (F01 * F20) / 2},
                                      {(F00 * F01) / 2,
                                       std::pow(F00, 2.0) / 2 + 2 * std::pow(F01, 2.0),
                                       (F00 * F11) / 2,
                                       (F00 * F10) / 2 + 2 * F01 * F11,
                                       (F00 * F21) / 2,
                                       (F00 * F20) / 2 + 2 * F01 * F21},
                                      {2 * F00 * F10 + (F01 * F11) / 2,
                                       (F00 * F11) / 2,
                                       2 * std::pow(F10, 2.0) + std::pow(F11, 2.0) / 2,
                                       (F10 * F11) / 2,
                                       2 * F10 * F20 + (F11 * F21) / 2,
                                       (F11 * F20) / 2},
                                      {(F01 * F10) / 2,
                                       (F00 * F10) / 2 + 2 * F01 * F11,
                                       (F10 * F11) / 2,
                                       std::pow(F10, 2.0) / 2 + 2 * std::pow(F11, 2.0),
                                       (F10 * F21) / 2,
                                       (F10 * F20) / 2 + 2 * F11 * F21},
                                      {2 * F00 * F20 + (F01 * F21) / 2,
                                       (F00 * F21) / 2,
                                       2 * F10 * F20 + (F11 * F21) / 2,
                                       (F10 * F21) / 2,
                                       2 * std::pow(F20, 2.0) + std::pow(F21, 2.0) / 2,
                                       (F20 * F21) / 2},
                                      {(F01 * F20) / 2,
                                       (F00 * F20) / 2 + 2 * F01 * F21,
                                       (F11 * F20) / 2,
                                       (F10 * F20) / 2 + 2 * F11 * F21,
                                       (F20 * F21) / 2,
                                       std::pow(F20, 2.0) / 2 + 2 * std::pow(F21, 2.0)},
                                  });

            m_dissipationHessians.push_back({dE2_X, dE2_V});
        }
    }
}

void VNCS::Sim3D::StVKForceField::addDForce(const sofa::core::MechanicalParams *mparams,
                                            DataVecDeriv &fDiffVec,
                                            const DataVecDeriv &xDiffVec)
{
    const auto &dxis = make_read_accessor(xDiffVec);
    auto &dfs = make_write_accessor(fDiffVec);
    for (int i = 0; i < dxis.size(); ++i) {
        const auto &hessian = m_dissipationHessians[i];
        const auto &dxi = dxis[i];
        auto &df = dfs[i];

        Eigen::Matrix<VNCS::Space3D::Real, 6, 1> dxi_;
        dxi_[0] = dxi(0, 0);
        dxi_[1] = dxi(0, 1);
        dxi_[2] = dxi(1, 0);
        dxi_[3] = dxi(1, 1);
        dxi_[4] = dxi(2, 0);
        dxi_[5] = dxi(2, 1);

        const Eigen::Matrix<VNCS::Space3D::Real, 6, 1> df_ = mparams->kFactor() * m_hessians[i] * dxi_;
        df(0, 0) -= df_[0];
        df(0, 1) -= df_[1];
        df(1, 0) -= df_[2];
        df(1, 1) -= df_[3];
        df(2, 0) -= df_[4];
        df(2, 1) -= df_[5];

        // Dissipation
        if (m_beta) {
            Eigen::Matrix<VNCS::Real, 6, 1> dStress =
                (mparams->kFactor() * hessian.first + mparams->bFactor() * hessian.second) * dxi_;
            df(0, 0) -= dStress[0];
            df(0, 1) -= dStress[1];
            df(1, 0) -= dStress[2];
            df(1, 1) -= dStress[3];
            df(2, 0) -= dStress[4];
            df(2, 1) -= dStress[5];
        }
    }
}

VNCS::Space3D::Real VNCS::Sim3D::StVKForceField::mu() const
{
    return m_mu;
}

void VNCS::Sim3D::StVKForceField::setMu(const Space3D::Real &mu)
{
    m_mu = mu;
}

VNCS::Space3D::Real VNCS::Sim3D::StVKForceField::lambda() const
{
    return m_lambda;
}

void VNCS::Sim3D::StVKForceField::setLambda(const Space3D::Real &lambda)
{
    m_lambda = lambda;
}

std::shared_ptr<VNCS::SamplingPoints<VNCS::Space2D>> VNCS::Sim3D::StVKForceField::samplingPoints() const
{
    return m_samplingPoints;
}

void VNCS::Sim3D::StVKForceField::setSamplingPoints(
    const std::shared_ptr<VNCS::SamplingPoints<VNCS::Space2D>> &samplingPoints)
{
    m_samplingPoints = samplingPoints;
}

void VNCS::Sim3D::StVKForceField::setBeta(const Real &beta)
{
    m_beta = beta;
}
VNCS::Sim3D::StVKForceField::Real VNCS::Sim3D::StVKForceField::beta() const
{
    return m_beta.value_or(0.0);
}