#include "StVKForceField.h"

#include <VNCS/DataExtensions.h>
#include <range/v3/view/zip.hpp>

VNCS::Sim2D::StVKForceField::StVKForceField()
    : Inherit()
{
}

void VNCS::Sim2D::StVKForceField::init()
{
    muLame = 0.5 * m_youngModulus / (1.0 + m_poissonRatio);
    lambdaLame = (m_youngModulus * m_poissonRatio) /
                 ((1.0 + m_poissonRatio) * (1.0 - 2.0 * m_poissonRatio));  // Plane strain condition (thick)
    // const Scalar lambdaLame = (young * poisson) / (1.0 - poisson*poisson); //
    // Plane stress condition (thin)

    mu = muLame;
    lambda = muLame + lambdaLame;

    alpha = 1.0 + (mu / lambda);

    Inherit::init();
}

VNCS::Space2D::Real VNCS::Sim2D::StVKForceField::getPotentialEnergy(const sofa::core::MechanicalParams *mparams,
                                                                    const DataVecCoord &xVec) const
{
    auto energy = Real(0.0);
    auto dissipation = Real(0.0);

    const auto dt = getContext()->getDt();
    const auto vVec = getMState()->readVelocities();

    const auto &samplers = *m_samplingPoints;
    for (const auto &[F, dFdt, sampler] :
         ranges::views::zip(make_read_accessor(xVec), make_read_accessor(vVec), samplers)) {
        const auto w = sampler.w;

        const auto Fn = F.squaredNorm();
        const auto dec = F.determinant() - alpha;

        energy += 0.5 * w * (mu * Fn + lambda * dec * dec - (mu * (1 + alpha)));

        if (m_beta) {
            const auto beta = m_beta.value();
            const Eigen::Matrix2d strainRate = 1.0 / 2.0 * (F.transpose() * dFdt + dFdt.transpose() * F);
            const auto e00 = strainRate(0, 0);
            const auto e11 = strainRate(1, 1);
            const auto e10 = strainRate(1, 0);

            dissipation += 0.5 * w * beta * (e00 * e00 + e11 * e11 + e10 * e10);
        }
    }

    return energy;  // + dt * dissipation;
}

void VNCS::Sim2D::StVKForceField::addForce(const sofa::core::MechanicalParams *mparams,
                                           DataVecDeriv &fVec,
                                           const DataVecCoord &xVec,
                                           const DataVecDeriv &vVec)
{
    const auto &samplers = *m_samplingPoints;

    m_dissipationHessians.clear();
    for (const auto &[F, dFdt, sampler, f] : ranges::views::zip(make_read_accessor(xVec),  //
                                                                make_read_accessor(vVec),
                                                                samplers,  //
                                                                make_write_accessor(fVec))) {
        const auto w = sampler.w;
        const Real dec = F.determinant() - alpha;

        Eigen::Matrix2d p_dec_p_F;
        p_dec_p_F << F(1, 1), -F(1, 0), -F(0, 1), F(0, 0);

        f += w * (-mu * F - lambda * dec * p_dec_p_F);

        // Dissipation
        if (m_beta) {
            const auto beta = m_beta.value();
            const auto F00 = F(0, 0);
            const auto F01 = F(0, 1);
            const auto F10 = F(1, 0);
            const auto F11 = F(1, 1);

            const auto dFdt00 = dFdt(0, 0);
            const auto dFdt01 = dFdt(0, 1);
            const auto dFdt10 = dFdt(1, 0);
            const auto dFdt11 = dFdt(1, 1);

            const auto diff_E_dFdt00 =
                2.0 * F00 * (F00 * dFdt00 + F10 * dFdt10) +
                F01 * ((F00 * dFdt01) / 2.0 + (F01 * dFdt00) / 2.0 + (F10 * dFdt11) / 2.0 + (F11 * dFdt10) / 2.0);
            const auto diff_E_dFdt11 =
                2.0 * F11 * (F01 * dFdt01 + F11 * dFdt11) +
                F10 * ((F00 * dFdt01) / 2.0 + (F01 * dFdt00) / 2.0 + (F10 * dFdt11) / 2.0 + (F11 * dFdt10) / 2.0);

            const auto diff_E_dFdt10 =
                2.0 * F10 * (F00 * dFdt00 + F10 * dFdt10) +
                F11 * ((F00 * dFdt01) / 2.0 + (F01 * dFdt00) / 2.0 + (F10 * dFdt11) / 2.0 + (F11 * dFdt10) / 2.0);

            const auto diff_E_dFdt01 =
                2.0 * F01 * (F01 * dFdt01 + F11 * dFdt11) +
                F00 * ((F00 * dFdt01) / 2.0 + (F01 * dFdt00) / 2.0 + (F10 * dFdt11) / 2.0 + (F11 * dFdt10) / 2.0);

            const auto dE00_dFdt = Eigen::Matrix2d({{F(0, 0), 0}, {F(1, 0), 0}});
            const auto dE11_dFdt = Eigen::Matrix2d({{0, F(0, 1)}, {0, F(1, 1)}});
            const auto dE10_dFdt = (1.0 / 2.0 * Eigen::Matrix2d({{F(0, 1), F(0, 0)}, {F(1, 1), F(1, 0)}})).eval();

            const Eigen::Matrix2d strainRate = 1.0 / 2.0 * (F.transpose() * dFdt + dFdt.transpose() * F);
            const auto e00 = strainRate(0, 0);
            const auto e11 = strainRate(1, 1);
            const auto e10 = strainRate(1, 0);

            const auto rStress = Eigen::Matrix2d({{diff_E_dFdt00, diff_E_dFdt01}, {diff_E_dFdt10, diff_E_dFdt11}});

            f -= w * beta * rStress;

            // Hessians
            // clang-format off
                Hessian dE2_wrt_X = w * beta * Hessian({
                    {4*F00*dFdt00 + (F01*dFdt01)/2 + 2*F10*dFdt10, (F00*dFdt01)/2 + F01*dFdt00 + (F10*dFdt11)/2 + (F11*dFdt10)/2, 2*F00*dFdt10 + (F01*dFdt11)/2, (F01*dFdt10)/2}, //
                    {F00*dFdt01 + (F01*dFdt00)/2 + (F10*dFdt11)/2 + (F11*dFdt10)/2, (F00*dFdt00)/2 + 4*F01*dFdt01 + 2*F11*dFdt11, (F00*dFdt11)/2, (F00*dFdt10)/2 + 2*F01*dFdt11}, //
                    {2*F10*dFdt00 + (F11*dFdt01)/2, (F11*dFdt00)/2, 2*F00*dFdt00 + 4*F10*dFdt10 + (F11*dFdt11)/2, (F00*dFdt01)/2 + (F01*dFdt00)/2 + (F10*dFdt11)/2 + F11*dFdt10}, //
                    {(F10*dFdt01)/2, (F10*dFdt00)/2 + 2*F11*dFdt01, (F00*dFdt01)/2 + (F01*dFdt00)/2 + F10*dFdt11 + (F11*dFdt10)/2, 2*F01*dFdt01 + (F10*dFdt10)/2 + 4*F11*dFdt11}  //
                });

                Hessian dE2_wrt_V = w * beta *  Hessian({
                    {3*std::pow(F00,2.0) + std::pow(F01,2.0)/2, (F00*F01)/2, 2*F00*F10 + (F01*F11)/2, (F01*F10)/2},
                    {(F00*F01)/2, std::pow(F00,2.0)/2 + 2*std::pow(F01,2.0), (F00*F11)/2, (F00*F10)/2 + 2*F01*F11},
                    {2*F00*F10 + (F01*F11)/2, (F00*F11)/2, 2*std::pow(F10,2.0) + std::pow(F11,2.0)/2, (F10*F11)/2},
                    {(F01*F10)/2, (F00*F10)/2 + 2*F01*F11, (F10*F11)/2, std::pow(F10,2.0)/2 + 2*std::pow(F11, 2.0)}
                });
            // clang-format on

            m_dissipationHessians.push_back({dE2_wrt_X, dE2_wrt_V});
        }
    }
}

void VNCS::Sim2D::StVKForceField::addDForce(const sofa::core::MechanicalParams *mparams,
                                            DataVecDeriv &fDiffVec,
                                            const DataVecDeriv &xDiffVec)
{
    auto mstate = this->getMState();

    const auto &samplers = *m_samplingPoints;
    const auto &Fs = mstate->readPositions();

    // If beta is not defined, the dissipation hessians is empty!
    if (m_beta) {
        for (const auto &[F, dxi, sampler, df, hessian] : ranges::views::zip(Fs,                            //
                                                                             make_read_accessor(xDiffVec),  //
                                                                             samplers,                      //
                                                                             make_write_accessor(fDiffVec),
                                                                             m_dissipationHessians)) {
            const auto w = sampler.w;

            const Real dec = F.determinant() - alpha;

            Eigen::Matrix2d p_dec_p_F;
            p_dec_p_F << F(1, 1), -F(1, 0), -F(0, 1), F(0, 0);

            Eigen::Matrix2d p_dec_p_FDiff;
            p_dec_p_FDiff << dxi(1, 1), -dxi(1, 0), -dxi(0, 1), dxi(0, 0);

            const Real decDiff = p_dec_p_F.cwiseProduct(dxi).sum();

            df += mparams->kFactor() * w * (-mu * dxi - lambda * (decDiff * p_dec_p_F + dec * p_dec_p_FDiff));

            // Dissipation
            Eigen::Matrix<VNCS::Real, 4, 1> dF = {dxi(0, 0), dxi(0, 1), dxi(1, 0), dxi(1, 1)};

            Eigen::Matrix<VNCS::Real, 4, 1> dStress =
                (mparams->kFactor() * hessian.first + mparams->bFactor() * hessian.second) * dF;
            df(0, 0) -= dStress[0];
            df(0, 1) -= dStress[1];
            df(1, 0) -= dStress[2];
            df(1, 1) -= dStress[3];
        }
    } else {
        for (const auto &[F, dxi, sampler, df] : ranges::views::zip(Fs,                            //
                                                                    make_read_accessor(xDiffVec),  //
                                                                    samplers,                      //
                                                                    make_write_accessor(fDiffVec))) {
            const auto w = sampler.w;

            const Real dec = F.determinant() - alpha;

            Eigen::Matrix2d p_dec_p_F;
            p_dec_p_F << F(1, 1), -F(1, 0), -F(0, 1), F(0, 0);

            Eigen::Matrix2d p_dec_p_FDiff;
            p_dec_p_FDiff << dxi(1, 1), -dxi(1, 0), -dxi(0, 1), dxi(0, 0);

            const Real decDiff = p_dec_p_F.cwiseProduct(dxi).sum();

            df += mparams->kFactor() * w * (-mu * dxi - lambda * (decDiff * p_dec_p_F + dec * p_dec_p_FDiff));
        }
    }
}

VNCS::Space2D::Real VNCS::Sim2D::StVKForceField::poisson() const
{
    return m_poissonRatio;
}

void VNCS::Sim2D::StVKForceField::setPoisson(const VNCS::Space2D::Real &poisson)
{
    m_poissonRatio = poisson;
}

VNCS::Space2D::Real VNCS::Sim2D::StVKForceField::young() const
{
    return m_youngModulus;
}

void VNCS::Sim2D::StVKForceField::setYoung(const VNCS::Space2D::Real &young)
{
    m_youngModulus = young;
}

std::shared_ptr<VNCS::SamplingPoints<VNCS::Space2D>> VNCS::Sim2D::StVKForceField::samplingPoints() const
{
    return m_samplingPoints;
}

void VNCS::Sim2D::StVKForceField::setSamplingPoints(
    const std::shared_ptr<VNCS::SamplingPoints<VNCS::Space2D>> &samplingPoints)
{
    m_samplingPoints = samplingPoints;
}
void VNCS::Sim2D::StVKForceField::setBeta(const Real &beta)
{
    m_beta = beta;
}
VNCS::Sim2D::StVKForceField::Real VNCS::Sim2D::StVKForceField::beta() const
{
    return m_beta.value_or(0.0);
}
