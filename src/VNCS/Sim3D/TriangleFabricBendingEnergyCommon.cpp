#include "TriangleFabricBendingEnergyCommon.h"

#include <VNCS/Spaces.h>
#include <VNCS/DataExtensions.h>
#include <VNCS/Logger.h>

#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <sofa/core/VecId.h>
#include <spdlog/spdlog.h>

#include <range/v3/view/enumerate.hpp>

using namespace VNCS::Sim3D;
namespace
{
const auto make_eigen = [](const VNCS::Space3D::VecType::Coord &v) -> Eigen::Vector3d { return {v[0], v[1], v[2]}; };
const auto make_coord = [](const Eigen::Vector3d &v) -> VNCS::Space3D::VecType::Coord { return {v[0], v[1], v[2]}; };
const auto make_deriv = [](const Eigen::Vector3d &v) -> VNCS::Space3D::VecType::Deriv { return {v[0], v[1], v[2]}; };
}  // namespace

void TriangleFabricBendingEnergy::init()
{
    Inherit::init();

    std::vector<VNCS::Space3D::Point> coarsePoints;
    std::vector<std::vector<std::size_t>> coarseFaces;
    std::ifstream coarseMeshIn(m_meshFilepath);
    CGAL::read_OBJ(coarseMeshIn, coarsePoints, coarseFaces);

    namespace PMP = CGAL::Polygon_mesh_processing;
    VNCS::Space3D::Mesh mesh;
    PMP::polygon_soup_to_polygon_mesh(coarsePoints, coarseFaces, mesh);

    computeHinges(mesh);
    computeRestPhis();
}

TriangleFabricBendingEnergy::Real TriangleFabricBendingEnergy::getPotentialEnergy(const sofa::core::MechanicalParams *,
                                                                                  const DataVecCoord &x) const
{
    const auto &positions = make_read_accessor(x);
    const auto &v = *getMState()->read(sofa::core::ConstVecDerivId::velocity());
    const auto &velocities = make_read_accessor(v);

    auto energy = 0.0;
    auto dissipation = 0.0;
    for (const auto &stencil : m_stencils) {
        energy += computeEnergy(stencil, positions);
        dissipation += computeDissipationEnergy(stencil, positions, velocities);
    }

    spdlog::info("Bending energy = {}", energy);

    return energy + getContext()->getDt() * dissipation;
}

void TriangleFabricBendingEnergy::addForce(const sofa::core::MechanicalParams *,
                                           DataVecDeriv &f,
                                           const DataVecCoord &x,
                                           const DataVecDeriv &v)
{
    const auto &positions = make_read_accessor(x);
    const auto &velocities = make_read_accessor(v);
    auto &forces = make_write_accessor(f);

    m_hessians.clear();
    for (const auto &stencil : m_stencils) {
        const auto f = accumulate_dEdxStencil(stencil, positions);

        const auto fDamping = accumulate_dEdvStencil(stencil, positions, velocities);

        forces[stencil.nodeIndices()[0]] += make_deriv(f.segment<3>(0) + fDamping.segment<3>(0));
        forces[stencil.nodeIndices()[1]] += make_deriv(f.segment<3>(3) + fDamping.segment<3>(3));
        forces[stencil.nodeIndices()[2]] += make_deriv(f.segment<3>(6) + fDamping.segment<3>(6));
        forces[stencil.nodeIndices()[3]] += make_deriv(f.segment<3>(9) + fDamping.segment<3>(9));

        m_hessians.push_back(std::make_pair(accumulate_d2Edx2Stencil(stencil, positions, velocities),
                                            accumulate_d2Edv2Stencil(stencil, positions, velocities)));
    }
}

void TriangleFabricBendingEnergy::addDForce(const sofa::core::MechanicalParams *mparams,
                                            DataVecDeriv &fDiffVec,
                                            const DataVecDeriv &xDiffVec)
{
    const auto &dxis = make_read_accessor(xDiffVec);
    auto &dfs = make_write_accessor(fDiffVec);
    for (const auto &[stencilIdx, stencil] : ranges::views::enumerate(m_stencils)) {
        Eigen::Matrix<VNCS::Real, 12, 1> dxi;

        dxi.segment<3>(0) = make_eigen(dxis[stencil.nodeIndices()[0]]);
        dxi.segment<3>(3) = make_eigen(dxis[stencil.nodeIndices()[1]]);
        dxi.segment<3>(6) = make_eigen(dxis[stencil.nodeIndices()[2]]);
        dxi.segment<3>(9) = make_eigen(dxis[stencil.nodeIndices()[3]]);

        const Eigen::Matrix<VNCS::Real, 12, 1> df =
            (mparams->kFactor() * m_hessians[stencilIdx].first + mparams->bFactor() * m_hessians[stencilIdx].second) *
            dxi;

        dfs[stencil.nodeIndices()[0]] += make_deriv(df.segment<3>(0));
        dfs[stencil.nodeIndices()[1]] += make_deriv(df.segment<3>(3));
        dfs[stencil.nodeIndices()[2]] += make_deriv(df.segment<3>(6));
        dfs[stencil.nodeIndices()[3]] += make_deriv(df.segment<3>(9));
    }
}

VNCS::Real TriangleFabricBendingEnergy::calculateAngleUnitVectors(const Eigen::Vector3d &v1,
                                                                  const Eigen::Vector3d &v2,
                                                                  const Eigen::Vector3d &axis)
{
    // axis does not need to be normalized as we only need it for the sign
    const VNCS::Real sign = ((v1.cross(v2)).dot(axis)) >= 0 ? 1.0 : -1.0;
    return 2. * sign * (std::atan2((v1 - v2).norm(), (v1 + v2).norm()));
}

VNCS::Real TriangleFabricBendingEnergy::computeEffectiveStiffness(const VNCS::Real kBending)
{
    return kBending;
}

void TriangleFabricBendingEnergy::computeRestPhis()
{
    /* -------------------------------------------------------------------------------
     * Note: rest angles different from 0 aren't supported by the damping implementation
     * ---------------------------------------------------------------------------- */

    const auto nbStencils = m_stencils.size();
    const auto &restPositions = getMState()->readRestPositions();

    for (int stencilIdx = 0; stencilIdx < nbStencils; stencilIdx++) {
        TriangleFabricBendingStencil &stencil = m_stencils[stencilIdx];

        const Eigen::Vector3d x0 = make_eigen(restPositions[stencil.nodeIndices()[0]]);
        const Eigen::Vector3d x1 = make_eigen(restPositions[stencil.nodeIndices()[1]]);
        const Eigen::Vector3d x2 = make_eigen(restPositions[stencil.nodeIndices()[2]]);
        const Eigen::Vector3d x3 = make_eigen(restPositions[stencil.nodeIndices()[3]]);

        const Eigen::Vector3d e0 = x2 - x1;
        const Eigen::Vector3d e1 = x0 - x2;
        const Eigen::Vector3d e1_tilde = x3 - x2;

        Eigen::Vector3d n = e0.cross(e1);
        Eigen::Vector3d n_tilde = -e0.cross(e1_tilde);
        const VNCS::Real n_length = n.norm();
        const VNCS::Real n_tilde_length = n_tilde.norm();
        n /= n_length;
        n_tilde /= n_tilde_length;

#if defined(HINGE_BENDING_USE_TANTHETA)
        const VNCS::Real tan_half_theta = (n - n_tilde).norm() / (n + n_tilde).norm();
        VNCS::Real sign_angle = n.cross(n_tilde).dot(e0);
        if (std::isnan(sign_angle))
            sign_angle = 1.;
        else
            sign_angle = sign_angle > 0 ? 1. : -1.;

        stencil.restPhi = 2. * sign_angle * tan_half_theta;
#elif defined(HINGE_BENDING_USE_SINTHETA)
        const VNCS::Real theta = std::atan2((n.cross(n_tilde)).dot(e0.normalized()), n.dot(n_tilde));

        stencil.restPhi = std::sin(theta / 2.);
#else
        const VNCS::Real theta = std::atan2((n.cross(n_tilde)).dot(e0.normalized()), n.dot(n_tilde));

        stencil.restPhi = theta;
#endif
    }
}

const std::filesystem::path &TriangleFabricBendingEnergy::meshFilepath() const
{
    return m_meshFilepath;
}

void TriangleFabricBendingEnergy::setMeshFilepath(const std::filesystem::path &newMeshFilepath)
{
    m_meshFilepath = newMeshFilepath;
}

const std::vector<std::pair<TriangleFabricBendingEnergy::Hessian, TriangleFabricBendingEnergy::Hessian>>
    &TriangleFabricBendingEnergy::hessians() const
{
    return m_hessians;
}

VNCS::Real TriangleFabricBendingEnergy::bendingStiffness() const
{
    return m_bendingStiffness;
}

void TriangleFabricBendingEnergy::setBendingStiffness(VNCS::Real newBendingStiffness)
{
    m_bendingStiffness = newBendingStiffness;
}

void TriangleFabricBendingEnergy::compute_dPsi(const Eigen::Vector3d &n1,
                                               const Eigen::Vector3d &n2,
                                               const Eigen::Vector3d &e0,
                                               const VNCS::Real rest_phi,
                                               const VNCS::Real ka,
                                               VNCS::Real &dPsi_dTheta,
                                               VNCS::Real &dPsi_dTheta_dTheta) const
{
#if defined(HINGE_BENDING_USE_TANTHETA)
    const VNCS::Real psi = (n1 - n2).norm() / (n1 + n2).norm();
    const VNCS::Real rest_psi = rest_phi * 0.5;
    VNCS::Real sign_angle = n1.cross(n2).dot(e0);
    if (std::isnan(sign_angle))
        sign_angle = 1.;
    else
        sign_angle = sign_angle > 0 ? 1. : -1.;
    const VNCS::Real phi_theta = 2. * sign_angle * psi;
    const VNCS::Real dPhi_dTheta = 1 + psi * psi;

    dPsi_dTheta = ka * (phi_theta - rest_phi) * dPhi_dTheta;

    dPsi_dTheta_dTheta = ka * dPhi_dTheta *              //
                         (                               //
                             2 * (psi - rest_psi) * psi  //
                             + dPhi_dTheta               //
                         );                              //
#elif defined(HINGE_BENDING_USE_SINTHETA)
    const VNCS::Real theta = std::atan2((n1.cross(n2)).dot(e0.normalized()), n1.dot(n2));
    const VNCS::Real phi_theta = 2 * std::sin(theta / 2.);
    dPsi_dTheta = 2. * ka * (phi_theta - rest_phi) * (0.5 * std::cos(theta / 2.));
    dPsi_dTheta_dTheta = 2. * ka *
                         ((0.5 * std::cos(theta / 2.)) * (0.5 * std::cos(theta / 2.)) +
                          (phi_theta - rest_phi) * (0.5 * 0.5 * -std::sin(theta / 2.)));
#else
    const VNCS::Real theta = std::atan2((n1.cross(n2)).dot(e0.normalized()), n1.dot(n2));
    dPsi_dTheta = 2. * ka * (theta - rest_phi);
    dPsi_dTheta_dTheta = 2. * ka;
#endif
}

VNCS::Real TriangleFabricBendingEnergy::calculateKa(const TriangleFabricBendingStencil &stencil) const
{
    const VNCS::Real &e0RestLength = stencil.e0RestNorm;
    const VNCS::Real &effectiveStiffness = stencil.effStiffness;
    const VNCS::Real &restArea = stencil.restArea;
    const VNCS::Real e0_rest_sqnorm = e0RestLength * e0RestLength;
    const VNCS::Real ka = effectiveStiffness * 3. * e0_rest_sqnorm / restArea;
    return ka;
}

std::array<Eigen::Vector3d, 5> TriangleFabricBendingEnergy::calculateEs(const TriangleFabricBendingStencil &stencil,
                                                                        const VecCoord &vecCoord) const
{
    const auto &positions = make_read_accessor(vecCoord);

    const auto &restPositions = getMState()->readRestPositions();

    const Eigen::Vector3d &x2 =
        make_eigen(restPositions[stencil.nodeIndices()[0]] + positions[stencil.nodeIndices()[0]]);
    const Eigen::Vector3d &x0 =
        make_eigen(restPositions[stencil.nodeIndices()[1]] + positions[stencil.nodeIndices()[1]]);
    const Eigen::Vector3d &x1 =
        make_eigen(restPositions[stencil.nodeIndices()[2]] + positions[stencil.nodeIndices()[2]]);
    const Eigen::Vector3d &x3 =
        make_eigen(restPositions[stencil.nodeIndices()[3]] + positions[stencil.nodeIndices()[3]]);

    return {
        x1 - x0,  // e0
        x2 - x0,  // e1
        x3 - x0,  // e2
        x2 - x1,  // e3
        x3 - x1   // e4
    };
}

void TriangleFabricBendingEnergy::calculateAngleGradients(const Eigen::Vector3d &n1,
                                                          const VNCS::Real n1_length,
                                                          const Eigen::Vector3d &n2,
                                                          const VNCS::Real n2_length,
                                                          const std::array<Eigen::Vector3d, 5> &e,
                                                          Eigen::Matrix<VNCS::Real, 12, 1> &thetaGradients) const
{
    const Eigen::Vector3d &e0 = e[0];
    const Eigen::Vector3d &e1 = e[1];
    const Eigen::Vector3d &e2 = e[2];
    const Eigen::Vector3d &e3 = e[3];
    const Eigen::Vector3d &e4 = e[4];

    const VNCS::Real e0_length = e[0].norm();

    thetaGradients.segment(0, 3) = -e0_length / n1_length * n1.transpose();
    thetaGradients.segment(3, 3) =
        (-e0 / e0_length).dot(e3) / n1_length * n1.transpose() + (-e0 / e0_length).dot(e4) / n2_length * n2.transpose();

    thetaGradients.segment(6, 3) =
        (e0 / e0_length).dot(e1) / n1_length * n1.transpose() + (e0 / e0_length).dot(e2) / n2_length * n2.transpose();

    thetaGradients.segment(9, 3) = -e0_length / n2_length * n2.transpose();
}

VNCS::Real TriangleFabricBendingEnergy::computeEnergy(const TriangleFabricBendingStencil &stencil,
                                                      const VecCoord &positions) const
{
    const VNCS::Real ka = calculateKa(stencil);
    const std::array<Eigen::Vector3d, 5> e = calculateEs(stencil, positions);

    // Analytic forces

    Eigen::Vector3d n1 = e[0].cross(e[3]);
    const VNCS::Real n1_length = n1.norm();
    Eigen::Vector3d n2 = -e[0].cross(e[4]);
    const VNCS::Real n2_length = n2.norm();
    n1 /= n1_length;
    n2 /= n2_length;

    const auto e0 = e[0];
    const VNCS::Real theta = std::atan2((n1.cross(n2)).dot(e0.normalized()), n1.dot(n2));
    auto psi = ka * std::pow(theta - stencil.restPhi, 2.0);

    return psi;
}

VNCS::Real TriangleFabricBendingEnergy::computeDissipationEnergy(const TriangleFabricBendingStencil &stencil,
                                                                 const VecCoord &positions,
                                                                 const VecDeriv &velocities) const
{
    const Eigen::Vector3d &v2 = make_eigen(velocities[stencil.nodeIndices()[0]]);
    const Eigen::Vector3d &v0 = make_eigen(velocities[stencil.nodeIndices()[1]]);
    const Eigen::Vector3d &v1 = make_eigen(velocities[stencil.nodeIndices()[2]]);
    const Eigen::Vector3d &v3 = make_eigen(velocities[stencil.nodeIndices()[3]]);

    const VNCS::Real ka = calculateKa(stencil);
    const std::array<Eigen::Vector3d, 5> e = calculateEs(stencil, positions);

    Eigen::Vector3d n1 = e[0].cross(e[3]);
    const VNCS::Real n1_length = n1.norm();
    Eigen::Vector3d n2 = -e[0].cross(e[4]);
    const VNCS::Real n2_length = n2.norm();
    n1 /= n1_length;
    n2 /= n2_length;

    Eigen::Matrix<VNCS::Real, 12, 1> thetaGradients;
    calculateAngleGradients(n1, n1_length, n2, n2_length, e, thetaGradients);

    const VNCS::Real thetaRate = thetaGradients.segment(0, 3).dot(v2) + thetaGradients.segment(3, 3).dot(v0) +
                                 thetaGradients.segment(6, 3).dot(v1) + thetaGradients.segment(9, 3).dot(v3);
    return -2.0 * m_dampingStiffness * ka * thetaRate * thetaRate;
}

// Forces

Eigen::Matrix<VNCS::Real, 12, 1> TriangleFabricBendingEnergy::calculateStencilForce(
    const TriangleFabricBendingStencil &stencil,
    const VecCoord &positions) const
{
    Eigen::Matrix<VNCS::Real, 12, 1> stencilForce;
    const VNCS::Real ka = calculateKa(stencil);
    const std::array<Eigen::Vector3d, 5> e = calculateEs(stencil, positions);

    // Analytic forces

    Eigen::Vector3d n1 = e[0].cross(e[3]);
    const VNCS::Real n1_length = n1.norm();
    Eigen::Vector3d n2 = -e[0].cross(e[4]);
    const VNCS::Real n2_length = n2.norm();
    n1 /= n1_length;
    n2 /= n2_length;

    VNCS::Real dPsi_dTheta;
    VNCS::Real dPsi_dTheta_dTheta;

    compute_dPsi(n1, n2, e[0], stencil.restPhi, ka, dPsi_dTheta, dPsi_dTheta_dTheta);

    const VNCS::Real common = -dPsi_dTheta;

    calculateAngleGradients(n1, n1_length, n2, n2_length, e, stencilForce);
    stencilForce *= common;
    return stencilForce;
}

Eigen::Matrix<VNCS::Real, 12, 1> TriangleFabricBendingEnergy::calculateStencilDampingForce(
    const TriangleFabricBendingStencil &stencil,
    const VecCoord &positions,
    const VecCoord &velocities) const
{
    const Eigen::Vector3d &v2 = make_eigen(velocities[stencil.nodeIndices()[0]]);
    const Eigen::Vector3d &v0 = make_eigen(velocities[stencil.nodeIndices()[1]]);
    const Eigen::Vector3d &v1 = make_eigen(velocities[stencil.nodeIndices()[2]]);
    const Eigen::Vector3d &v3 = make_eigen(velocities[stencil.nodeIndices()[3]]);

    const VNCS::Real ka = calculateKa(stencil);
    const std::array<Eigen::Vector3d, 5> e = calculateEs(stencil, positions);

    Eigen::Vector3d n1 = e[0].cross(e[3]);
    const VNCS::Real n1_length = n1.norm();
    Eigen::Vector3d n2 = -e[0].cross(e[4]);
    const VNCS::Real n2_length = n2.norm();
    n1 /= n1_length;
    n2 /= n2_length;

    Eigen::Matrix<VNCS::Real, 12, 1> thetaGradients;
    calculateAngleGradients(n1, n1_length, n2, n2_length, e, thetaGradients);

    const VNCS::Real thetaRate = thetaGradients.segment(0, 3).dot(v2) + thetaGradients.segment(3, 3).dot(v0) +
                                 thetaGradients.segment(6, 3).dot(v1) + thetaGradients.segment(9, 3).dot(v3);
    return -2.0 * m_dampingStiffness * ka * thetaRate * thetaGradients;
}

Eigen::Matrix<VNCS::Real, 12, 1> TriangleFabricBendingEnergy::accumulate_dEdxStencil(
    const TriangleFabricBendingStencil &stencil,
    const VecCoord &positions) const
{
    const auto stencilForce = calculateStencilForce(stencil, positions);
    return stencilForce;
}

Eigen::Matrix<VNCS::Real, 12, 1> TriangleFabricBendingEnergy::accumulate_dEdvStencil(
    const TriangleFabricBendingStencil &stencil,
    const VecCoord &positions,
    const VecDeriv &velocities) const
{
    const auto stencilForce = calculateStencilDampingForce(stencil, positions, velocities);
    return stencilForce;
}

void TriangleFabricBendingEnergy::computeHinges(const VNCS::Space3D::Mesh &mesh)
{
    m_stencils.reserve(mesh.number_of_edges());

    const auto &restPositions = getMState()->readRestPositions();

    for (const auto &edge : mesh.edges()) {
        if (!mesh.is_border(edge)) {
            const auto he0 = mesh.halfedge(edge, 0);
            const auto he1 = mesh.halfedge(edge, 1);

            const auto v0Id = mesh.source(he0);
            const auto v1Id = mesh.source(he1);

            const auto he2 = mesh.prev(he0);
            const auto he3 = mesh.prev(he1);

            const auto v2Id = mesh.source(he2);
            const auto v3Id = mesh.source(he3);

            const Eigen::Vector3d x0 = make_eigen(restPositions[v0Id.idx()]);
            const Eigen::Vector3d x1 = make_eigen(restPositions[v1Id.idx()]);
            const Eigen::Vector3d x2 = make_eigen(restPositions[v2Id.idx()]);
            const Eigen::Vector3d x3 = make_eigen(restPositions[v3Id.idx()]);

            const Eigen::Vector3d e0 = x1 - x0;
            const Eigen::Vector3d e1 = x2 - x1;
            const Eigen::Vector3d e1_tilde = x3 - x1;

            const int v0Idx = v0Id.idx();
            const int v1Idx = v1Id.idx();
            const int v2Idx = v2Id.idx();
            const int v3Idx = v3Id.idx();

            const std::array<int, 4> stNodes = {v2Idx, v0Idx, v1Idx, v3Idx};

            m_stencils.push_back(TriangleFabricBendingStencil{
                stNodes,                                                 // indices
                e0.norm(),                                               // e0RestNorm
                m_bendingStiffness,                                      // effStiffness
                0.5 * (e0.cross(e1).norm() + e0.cross(e1_tilde).norm())  // restArea
            });
        }
    }
}

// end Forces

// Hessians

void TriangleFabricBendingEnergy::calculateAngleHessian(const std::array<Eigen::Vector3d, 2> &ns,
                                                        const std::array<Eigen::Vector3d, 5> &e,
                                                        Eigen::Matrix<VNCS::Real, 12, 12> &theta_hessian) const
{
    const VNCS::Real e0_length = e[0].norm();
    const VNCS::Real e1_length = e[1].norm();
    const VNCS::Real e2_length = e[2].norm();
    const VNCS::Real e3_length = e[3].norm();
    const VNCS::Real e4_length = e[4].norm();

    const Eigen::Vector3d e0_normalized = e[0] / e0_length;
    const Eigen::Vector3d e1_normalized = e[1] / e1_length;
    const Eigen::Vector3d e2_normalized = e[2] / e2_length;
    const Eigen::Vector3d e3_normalized = e[3] / e3_length;
    const Eigen::Vector3d e4_normalized = e[4] / e4_length;

    // This is storable per vertex --> map with n cos per vertex where n is the number of incident triangles
    const VNCS::Real cos_alphas[4] = {e1_normalized.dot(e0_normalized),
                                      e2_normalized.dot(e0_normalized),
                                      e3_normalized.dot(-e0_normalized),
                                      e4_normalized.dot(-e0_normalized)

    };

    // this is storable per face
    const VNCS::Real A1 = 0.5 * (e[1]).cross(e[0]).norm();
    const VNCS::Real A2 = 0.5 * (e[2]).cross(e[0]).norm();

    // This is storable per half-edge (or vertex with an abstraction) --> two values per edge
    const VNCS::Real inv_hs[6] = {e1_length / (2. * A1),
                                  e2_length / (2. * A2),
                                  e3_length / (2. * A1),
                                  e4_length / (2. * A2),
                                  e0_length / (2. * A1),
                                  e0_length / (2. * A2)};

    const Eigen::Vector3d ms[6] = {-e1_normalized.cross(ns[0]).normalized(),   // 1
                                   e2_normalized.cross(ns[1]).normalized(),    // 2
                                   e3_normalized.cross(ns[0]).normalized(),    // 3
                                   -e4_normalized.cross(ns[1]).normalized(),   // 4
                                   e0_normalized.cross(ns[0]).normalized(),    // 01
                                   -e0_normalized.cross(ns[1]).normalized()};  // 02

    const Eigen::Matrix3d B1 = ns[0] * ms[4].transpose() / (e0_length * e0_length);
    const Eigen::Matrix3d B2 = ns[1] * ms[5].transpose() / (e0_length * e0_length);

    const auto l_S_operator = [](const Eigen::Matrix3d &m) { return m + m.transpose(); };

    const auto l_M_operator = [&cos_alphas, &inv_hs, &ms, &ns](
                                  const unsigned int i, const unsigned int j, const unsigned int k) {
        return cos_alphas[i] * inv_hs[i] * inv_hs[j] *  // scalars
               ms[j] * ns[k].transpose();               // outer product
    };

    const auto l_N_operator = [&inv_hs, &ms, &ns](const unsigned int i, const unsigned int j) {
        // 4+i because h01 & h02 are the last elements
        return inv_hs[4 + i] * inv_hs[j] *  // scalars
               ns[i] * ms[j].transpose();   // outer product
    };

    // Indices of matrix operators are shifted 1 less from Rasmus' notation to comply to computing conventions
    // row 0
    theta_hessian.block(3, 3, 3, 3) =               //
        l_S_operator(l_M_operator(2, 2, 0)) - B1 +  //
        l_S_operator(l_M_operator(3, 3, 1)) - B2;   //

    theta_hessian.block(3, 6, 3, 3) =                                     //
        l_M_operator(2, 0, 0) + l_M_operator(0, 2, 0).transpose() + B1 +  //
        l_M_operator(3, 1, 1) + l_M_operator(1, 3, 1).transpose() + B2;   //

    theta_hessian.block(3, 0, 3, 3) = l_M_operator(2, 4, 0) - l_N_operator(0, 2);
    theta_hessian.block(3, 9, 3, 3) = l_M_operator(3, 5, 1) - l_N_operator(1, 3);

    // row 1
    theta_hessian.block(6, 6, 3, 3) =
        l_S_operator(l_M_operator(0, 0, 0)) - B1 + l_S_operator(l_M_operator(1, 1, 1)) - B2;
    theta_hessian.block(6, 0, 3, 3) = l_M_operator(0, 4, 0) - l_N_operator(0, 0);
    theta_hessian.block(6, 9, 3, 3) = l_M_operator(1, 5, 1) - l_N_operator(1, 1);

    // row 2
    theta_hessian.block(0, 0, 3, 3) = -l_S_operator(l_N_operator(0, 4));
    theta_hessian.block(0, 9, 3, 3) = Eigen::Matrix3d::Zero();

    // row 3
    theta_hessian.block(9, 9, 3, 3) = -l_S_operator(l_N_operator(1, 5));

    // Symmetri
    theta_hessian.block(6, 3, 3, 3) = theta_hessian.block(3, 6, 3, 3).transpose();

    theta_hessian.block(0, 3, 3, 3) = theta_hessian.block(3, 0, 3, 3).transpose();
    theta_hessian.block(9, 3, 3, 3) = theta_hessian.block(3, 9, 3, 3).transpose();

    // row 1
    theta_hessian.block(0, 6, 3, 3) = theta_hessian.block(6, 0, 3, 3).transpose();
    theta_hessian.block(9, 6, 3, 3) = theta_hessian.block(6, 9, 3, 3).transpose();

    // row 2
    theta_hessian.block(9, 0, 3, 3) = theta_hessian.block(0, 9, 3, 3).transpose();
}

Eigen::Matrix<VNCS::Real, 12, 12> TriangleFabricBendingEnergy::calculateEdgeStencilHessian(
    const TriangleFabricBendingStencil &stencil,
    const VecCoord &positions) const
{
    const VNCS::Real ka = calculateKa(stencil);

    const std::array<Eigen::Vector3d, 5> e = calculateEs(stencil, positions);

    const Eigen::Vector3d n1_full = e[0].cross(e[3]);
    const Eigen::Vector3d n2_full = -e[0].cross(e[4]);
    const VNCS::Real n1_length = n1_full.norm();
    const VNCS::Real n2_length = n2_full.norm();

    const std::array<Eigen::Vector3d, 2> ns = {n1_full / n1_length, n2_full / n2_length};

    Eigen::Matrix<VNCS::Real, 12, 12> theta_hessian;
    calculateAngleHessian(ns, e, theta_hessian);

    VNCS::Real dPsi_dTheta;
    VNCS::Real dPsi_dTheta_dTheta;
    compute_dPsi(ns[0], ns[1], e[0], stencil.restPhi, ka, dPsi_dTheta, dPsi_dTheta_dTheta);

    Eigen::Matrix<VNCS::Real, 12, 1> theta_der;
    calculateAngleGradients(ns[0], n1_length, ns[1], n2_length, e, theta_der);

    Eigen::Matrix<VNCS::Real, 12, 12> stencilHessian;
    stencilHessian = -(dPsi_dTheta * theta_hessian + dPsi_dTheta_dTheta * theta_der * theta_der.transpose());
    return stencilHessian;
}

Eigen::Matrix<VNCS::Real, 12, 12> TriangleFabricBendingEnergy::calculateStencilDampingPositionHessian(
    const TriangleFabricBendingStencil &stencil,
    const VecCoord &positions,
    const VecDeriv &velocities) const
{
    const Eigen::Vector3d &v2 = make_eigen(velocities[stencil.nodeIndices()[0]]);
    const Eigen::Vector3d &v0 = make_eigen(velocities[stencil.nodeIndices()[1]]);
    const Eigen::Vector3d &v1 = make_eigen(velocities[stencil.nodeIndices()[2]]);
    const Eigen::Vector3d &v3 = make_eigen(velocities[stencil.nodeIndices()[3]]);
    Eigen::Matrix<VNCS::Real, 12, 1> velocityVector;
    velocityVector << v0, v1, v2, v3;

    const VNCS::Real ka = calculateKa(stencil);

    const std::array<Eigen::Vector3d, 5> e = calculateEs(stencil, positions);

    const Eigen::Vector3d n1_full = e[0].cross(e[3]);
    const Eigen::Vector3d n2_full = -e[0].cross(e[4]);
    const VNCS::Real n1_length = n1_full.norm();
    const VNCS::Real n2_length = n2_full.norm();

    const std::array<Eigen::Vector3d, 2> ns = {n1_full / n1_length, n2_full / n2_length};

    Eigen::Matrix<VNCS::Real, 12, 1> thetaGradients;
    calculateAngleGradients(ns[0], n1_length, ns[1], n2_length, e, thetaGradients);

    Eigen::Matrix<VNCS::Real, 12, 12> thetaHessian;
    calculateAngleHessian(ns, e, thetaHessian);
    const VNCS::Real thetaRate = thetaGradients.segment(0, 3).dot(v2) + thetaGradients.segment(3, 3).dot(v0) +
                                 thetaGradients.segment(6, 3).dot(v1) + thetaGradients.segment(9, 3).dot(v3);

    return -2. * ka * m_dampingStiffness * thetaRate * thetaHessian;
}

Eigen::Matrix<VNCS::Real, 12, 12> TriangleFabricBendingEnergy::calculateStencilDampingVelocityHessian(
    const TriangleFabricBendingStencil &stencil,
    const VecCoord &positions,
    const VecDeriv &velocities) const
{
    const VNCS::Real ka = calculateKa(stencil);

    const std::array<Eigen::Vector3d, 5> e = calculateEs(stencil, positions);

    const Eigen::Vector3d n1_full = e[0].cross(e[3]);
    const Eigen::Vector3d n2_full = -e[0].cross(e[4]);
    const VNCS::Real n1_length = n1_full.norm();
    const VNCS::Real n2_length = n2_full.norm();

    const Eigen::Vector3d ns[2] = {n1_full / n1_length, n2_full / n2_length};

    Eigen::Matrix<VNCS::Real, 12, 1> thetaGradients;
    calculateAngleGradients(ns[0], n1_length, ns[1], n2_length, e, thetaGradients);

    return -2. * ka * m_dampingStiffness * thetaGradients * thetaGradients.transpose();
}

Eigen::Matrix<VNCS::Real, 12, 12> TriangleFabricBendingEnergy::accumulate_d2Edx2Stencil(
    const TriangleFabricBendingStencil &stencil,
    const VecCoord &positions,
    const VecDeriv &velocities) const
{
    Eigen::Matrix<VNCS::Real, 12, 12> hess = calculateEdgeStencilHessian(stencil, positions);
    Eigen::Matrix<VNCS::Real, 12, 12> dampingHess =
        calculateStencilDampingPositionHessian(stencil, positions, velocities);
    return hess + dampingHess;
}

Eigen::Matrix<VNCS::Real, 12, 12> TriangleFabricBendingEnergy::accumulate_d2Edv2Stencil(
    const TriangleFabricBendingStencil &stencil,
    const VecCoord &positions,
    const VecDeriv &velocities) const
{
    return calculateStencilDampingVelocityHessian(stencil, positions, velocities);
}

// end Hessians

TriangleFabricBendingEnergy::TriangleFabricBendingStencil::TriangleFabricBendingStencil(const std::array<int, 4> &ids,
                                                                                        const Real e0RNorm,
                                                                                        const Real effStiff,
                                                                                        const Real rArea)
    : ids(ids)
    , e0RestNorm(e0RNorm)
    , effStiffness(effStiff)
    , restArea(rArea)
{
}

const std::array<int, 4> &TriangleFabricBendingEnergy::TriangleFabricBendingStencil::nodeIndices() const
{
    return ids;
}
void VNCS::Sim3D::TriangleFabricBendingEnergy::setBeta(VNCS::Real beta)
{
    m_dampingStiffness = beta;
}
VNCS::Real VNCS::Sim3D::TriangleFabricBendingEnergy::beta() const
{
    return m_dampingStiffness;
}
