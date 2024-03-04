#include "StableNeoHookean.h"

#include <VNCS/DataExtensions.h>
#include <range/v3/view/zip.hpp>

#include <spdlog/fmt/ostr.h>

namespace VNCS
{
namespace Sim3D
{
void StableNeoHookean::init()
{
    assert(m_mu > 0.0);
    assert(m_lambda > 0.0);
    assert(m_ratio > 0.0);

    m_ratio = m_mu / m_lambda;

    Inherit::init();

    spdlog::get("VNCS")->info("lambda = {}", m_lambda);
    spdlog::get("VNCS")->info("mu = {}", m_mu);
    spdlog::get("VNCS")->info("ratio = {}", m_ratio);
}

VNCS::Space3D::Real StableNeoHookean::psi(const VNCS::Space3D::F33::Coord &F) const
{
    const VNCS::Space3D::Real Ic = F.squaredNorm();
    const VNCS::Space3D::Real Jminus1 = F.determinant() - 1.0 - m_ratio;
    return 0.5 * (m_mu * (Ic - 3.0) + m_lambda * Jminus1 * Jminus1);
}

static Eigen::Matrix3d PartialJpartialF(const Eigen::Matrix3d &F)
{
    Eigen::Matrix3d pJpF;

    pJpF.col(0) = F.col(1).cross(F.col(2));
    pJpF.col(1) = F.col(2).cross(F.col(0));
    pJpF.col(2) = F.col(0).cross(F.col(1));

    return pJpF;
}

Eigen::Matrix3d StableNeoHookean::pk1(const VNCS::Space3D::F33::Coord &F) const
{
    const Eigen::Matrix3d pJpF = PartialJpartialF(F);
    const VNCS::Space3D::Real Jminus1 = F.determinant() - 1.0 - m_ratio;
    return m_mu * F + m_lambda * Jminus1 * pJpF;
}

// static void BuildTwistAndFlipEigenvectors(const Matrix3 &U, const Matrix3 &V, Matrix9 &Q)
//{
//    static const VNCS::Space3D::Real scale = 1.0 / std::sqrt(2.0);
//    const VNCS::S sV = scale * V;

//    using M3 = Eigen::Matrix<VNCS::Space3D::Real, 3, 3, Eigen::ColMajor>;

//    M3 A;
//    A << sV(0, 2) * U(0, 1), sV(1, 2) * U(0, 1), sV(2, 2) * U(0, 1), sV(0, 2) * U(1, 1), sV(1, 2) * U(1, 1),
//        sV(2, 2) * U(1, 1), sV(0, 2) * U(2, 1), sV(1, 2) * U(2, 1), sV(2, 2) * U(2, 1);

//    M3 B;
//    B << sV(0, 1) * U(0, 2), sV(1, 1) * U(0, 2), sV(2, 1) * U(0, 2), sV(0, 1) * U(1, 2), sV(1, 1) * U(1, 2),
//        sV(2, 1) * U(1, 2), sV(0, 1) * U(2, 2), sV(1, 1) * U(2, 2), sV(2, 1) * U(2, 2);

//    M3 C;
//    C << sV(0, 2) * U(0, 0), sV(1, 2) * U(0, 0), sV(2, 2) * U(0, 0), sV(0, 2) * U(1, 0), sV(1, 2) * U(1, 0),
//        sV(2, 2) * U(1, 0), sV(0, 2) * U(2, 0), sV(1, 2) * U(2, 0), sV(2, 2) * U(2, 0);

//    M3 D;
//    D << sV(0, 0) * U(0, 2), sV(1, 0) * U(0, 2), sV(2, 0) * U(0, 2), sV(0, 0) * U(1, 2), sV(1, 0) * U(1, 2),
//        sV(2, 0) * U(1, 2), sV(0, 0) * U(2, 2), sV(1, 0) * U(2, 2), sV(2, 0) * U(2, 2);

//    M3 E;
//    E << sV(0, 1) * U(0, 0), sV(1, 1) * U(0, 0), sV(2, 1) * U(0, 0), sV(0, 1) * U(1, 0), sV(1, 1) * U(1, 0),
//        sV(2, 1) * U(1, 0), sV(0, 1) * U(2, 0), sV(1, 1) * U(2, 0), sV(2, 1) * U(2, 0);

//    M3 F;
//    F << sV(0, 0) * U(0, 1), sV(1, 0) * U(0, 1), sV(2, 0) * U(0, 1), sV(0, 0) * U(1, 1), sV(1, 0) * U(1, 1),
//        sV(2, 0) * U(1, 1), sV(0, 0) * U(2, 1), sV(1, 0) * U(2, 1), sV(2, 0) * U(2, 1);

//    // Twist eigenvectors
//    Eigen::Map<M3>(Q.data()) = B - A;
//    Eigen::Map<M3>(Q.data() + 9) = D - C;
//    Eigen::Map<M3>(Q.data() + 18) = F - E;

//    // Flip eigenvectors
//    Eigen::Map<M3>(Q.data() + 27) = A + B;
//    Eigen::Map<M3>(Q.data() + 36) = C + D;
//    Eigen::Map<M3>(Q.data() + 45) = E + F;
//}

// static Matrix9 ProjectHessianWithAnalyticalFormulasNew(const Scalar &mu,
//                                                       const Scalar &lambda,
//                                                       const Matrix3 &F,
//                                                       const Matrix3 &U,
//                                                       const Matrix3 &V,
//                                                       const Vector3 &S)
//{
//    Vector9 eigenvalues;
//    Matrix9 eigenvectors;

//    const Scalar J = F.determinant();

//    // Compute the twist and flip eigenvalues
//    {
//        // Twist eigenvalues
//        eigenvalues.segment<3>(0) = S;
//        // Flip eigenvalues
//        eigenvalues.segment<3>(3) = -S;
//        const Scalar evScale = lambda * (J - 1.0) - mu;
//        eigenvalues.segment<6>(0) *= evScale;
//        eigenvalues.segment<6>(0).array() += mu;
//    }

//    // Compute the twist and flip eigenvectors
//    BuildTwistAndFlipEigenvectors(U, V, eigenvectors);

//    // Compute the remaining three eigenvalues and eigenvectors
//    {
//        Matrix3 A;
//        const Scalar s0s0 = S(0) * S(0);
//        const Scalar s1s1 = S(1) * S(1);
//        const Scalar s2s2 = S(2) * S(2);
//        A(0, 0) = mu + lambda * s1s1 * s2s2;
//        A(1, 1) = mu + lambda * s0s0 * s2s2;
//        A(2, 2) = mu + lambda * s0s0 * s1s1;
//        const Scalar evScale = lambda * (2.0 * J - 1.0) - mu;
//        A(0, 1) = evScale * S(2);
//        A(1, 0) = A(0, 1);
//        A(0, 2) = evScale * S(1);
//        A(2, 0) = A(0, 2);
//        A(1, 2) = evScale * S(0);
//        A(2, 1) = A(1, 2);

//        const Eigen::SelfAdjointEigenSolver<Matrix3> Aeigs(A);
//        eigenvalues.segment<3>(6) = Aeigs.eigenvalues();

//        Eigen::Map<Matrix3>(eigenvectors.data() + 54) = U * Aeigs.eigenvectors().col(0).asDiagonal() * V.transpose();
//        Eigen::Map<Matrix3>(eigenvectors.data() + 63) = U * Aeigs.eigenvectors().col(1).asDiagonal() * V.transpose();
//        Eigen::Map<Matrix3>(eigenvectors.data() + 72) = U * Aeigs.eigenvectors().col(2).asDiagonal() * V.transpose();
//    }

//    // Clamp the eigenvalues
//    for (int i = 0; i < 9; i++) {
//        if (eigenvalues(i) < 0.0) {
//            eigenvalues(i) = 0.0;
//        }
//    }

//    return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose();
//}

static Eigen::Matrix<VNCS::Space3D::Real, 9, 1> partialJpartialFVec(const Eigen::Matrix3d &F)
{
    Eigen::Matrix<VNCS::Space3D::Real, 9, 1> pJpF;
    pJpF.segment<3>(0) = F.col(1).cross(F.col(2));
    pJpF.segment<3>(3) = F.col(2).cross(F.col(0));
    pJpF.segment<3>(6) = F.col(0).cross(F.col(1));
    return pJpF;
}

static Eigen::Matrix3d crossProductMatrix(const Eigen::Matrix3d &F, const int col, const VNCS::Space3D::Real &scale)
{
    Eigen::Matrix3d cpm;
    cpm << 0, -scale * F(2, col), scale * F(1, col), scale * F(2, col), 0, -scale * F(0, col), -scale * F(1, col),
        scale * F(0, col), 0;
    return cpm;
}

static Eigen::Matrix<VNCS::Space3D::Real, 9, 9> computeFJSecondDerivContribs(const VNCS::Space3D::Real &lambda,
                                                                             const VNCS::Space3D::Real &ratio,
                                                                             const Eigen::Matrix3d &F)
{
    const VNCS::Space3D::Real scale = lambda * (F.determinant() - 1.0 - ratio);

    const Eigen::Matrix3d ahat = crossProductMatrix(F, 0, scale);
    const Eigen::Matrix3d bhat = crossProductMatrix(F, 1, scale);
    const Eigen::Matrix3d chat = crossProductMatrix(F, 2, scale);

    Eigen::Matrix<VNCS::Space3D::Real, 9, 9> FJ;
    FJ.block<3, 3>(0, 0).setZero();
    FJ.block<3, 3>(0, 3) = -chat;
    FJ.block<3, 3>(0, 6) = bhat;

    FJ.block<3, 3>(3, 0) = chat;
    FJ.block<3, 3>(3, 3).setZero();
    FJ.block<3, 3>(3, 6) = -ahat;

    FJ.block<3, 3>(6, 0) = -bhat;
    FJ.block<3, 3>(6, 3) = ahat;
    FJ.block<3, 3>(6, 6).setZero();

    return FJ;
}

Eigen::Matrix<VNCS::Space3D::Real, 9, 9> StableNeoHookean::partialPpartialF(const VNCS::Space3D::F33::Coord &F) const
{
    const Eigen::Matrix<VNCS::Space3D::Real, 9, 1> pjpf = partialJpartialFVec(F);
    return m_mu * Eigen::Matrix<VNCS::Space3D::Real, 9, 9>::Identity() + m_lambda * pjpf * pjpf.transpose() +
           computeFJSecondDerivContribs(m_lambda, m_ratio, F);
}

void StableNeoHookean::setSamplingPoints(const std::shared_ptr<VNCS::SamplingPoints<VNCS::Space3D>> &samplingPoints)
{
    m_samplingPoints = samplingPoints;
}

// Eigen::Matrix<VNCS::Space3D::Real, 9, 9> StableNeoHookean::clampedPartialPpartialF(
//    const VNCS::Space3D::F33::Coord &F) const
//{
//    return ProjectHessianWithAnalyticalFormulasNew(m_mu, m_lambda, F);
//}

VNCS::Space3D::Real StableNeoHookean::lambda() const
{
    return m_lambda;
}

void StableNeoHookean::setLambda(const VNCS::Space3D::Real &lambda)
{
    m_lambda = lambda;
}

VNCS::Space3D::Real StableNeoHookean::mu() const
{
    return m_mu;
}

void StableNeoHookean::setMu(const VNCS::Space3D::Real &mu)
{
    m_mu = mu;
}

VNCS::Space3D::Real StableNeoHookean::getPotentialEnergy(const sofa::core::MechanicalParams *,
                                                         const DataVecCoord &xVec) const
{
    auto energy = Real(0.0);

    const auto &samplers = *m_samplingPoints;
    for (const auto &[F, sampler] : ranges::views::zip(make_read_accessor(xVec), samplers)) {
        const auto w = sampler.w;
        energy += w * this->psi(F);
    }

    return energy;
}

void StableNeoHookean::addForce(const sofa::core::MechanicalParams *mparams,
                                DataVecDeriv &fVec,
                                const DataVecCoord &xVec,
                                const DataVecDeriv &vVec)
{
    const auto &samplers = *m_samplingPoints;
    const auto &Fs = make_read_accessor(xVec);
    auto &fs = make_write_accessor(fVec);
    m_hessians.resize(Fs.size());
    const auto lame1 = m_lambda;
    const auto lame2 = m_mu;
    for (int i = 0; i < Fs.size(); ++i) {
        const auto &F = Fs[i];
        auto &f = fs[i];
        const auto &sampler = samplers[i];
        const auto w = sampler.w;

        f -= w * this->pk1(F);
        m_hessians[i] = w * this->partialPpartialF(F);
        //const auto &F = Fs[i];
        //auto &f = fs[i];
        //const auto &sampler = samplers[i];
        //const auto w = sampler.w;

        //std::array<VNCS::Real, 9> Fv;
        //Fv[0] = F(0, 0);
        //Fv[1] = F(1, 0);
        //Fv[2] = F(2, 0);
        //Fv[3] = F(0, 1);
        //Fv[4] = F(1, 1);
        //Fv[5] = F(2, 1);
        //Fv[6] = F(0, 2);
        //Fv[7] = F(1, 2);
        //Fv[8] = F(2, 2);

        //{
        //    std::array<VNCS::Real, 9> vgx;
//#include "FEMVol_StVK_Gradient.mcg"

        //    Eigen::Map<Eigen::Matrix<VNCS::Real, 3, 3>> vgxMap(vgx.data());

        //    f -= w * vgxMap;
        //}

        //{
        //    std::array<std::array<VNCS::Real, 9>, 9> mHx;
//#include "FEMVol_StVK_Hessian.mcg"

        //    Eigen::Map<Eigen::Matrix<VNCS::Real, 9, 9>> hessianMap(mHx[0].data());
        //    m_hessians[i] = w * hessianMap;
        //}
    }
}

void StableNeoHookean::addDForce(const sofa::core::MechanicalParams *mparams,
                                 DataVecDeriv &fDiffVec,
                                 const DataVecDeriv &xDiffVec)
{
    const auto &samplers = *m_samplingPoints;
    const auto &dxis = make_read_accessor(xDiffVec);
    auto &dfs = make_write_accessor(fDiffVec);
    for (int i = 0; i < dxis.size(); ++i) {
        const auto &dxi = dxis[i];
        auto &df = dfs[i];

        Eigen::Map<Eigen::Matrix<VNCS::Space3D::Real, 9, 1>> df_(df.data());
        Eigen::Map<const Eigen::Matrix<VNCS::Space3D::Real, 9, 1>> dxi_(dxi.data());

        df_ -= mparams->kFactor() * m_hessians[i] * dxi_;
    }
}

}  // namespace Sim3D
}  // namespace VNCS
