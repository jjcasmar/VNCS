#ifndef VNCS_SIM3D_TRIANGLEFABRICBENDINGENERGYCOMMON_H
#define VNCS_SIM3D_TRIANGLEFABRICBENDINGENERGYCOMMON_H

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/ForceField.inl>

#include <filesystem>

#include <VNCS/Spaces.h>
#include <Eigen/Core>

namespace VNCS
{
namespace Sim3D
{
class TriangleFabricBendingEnergy : public sofa::core::behavior::ForceField<VNCS::Space3D::VecType>
{
    using Inherit = sofa::core::behavior::ForceField<VNCS::Space3D::VecType>;
    using Real = VNCS::Space3D::Real;
    using DataVecCoord = Inherit::DataVecCoord;
    using DataVecDeriv = Inherit::DataVecDeriv;

public:
    struct TriangleFabricBendingStencil {
        explicit TriangleFabricBendingStencil(const std::array<int, 4> &ids,
                                              const VNCS::Real e0RNorm,
                                              const VNCS::Real effStiff,
                                              const VNCS::Real rArea);

        const std::array<int, 4> &nodeIndices() const;

        const std::array<int, 4> ids;
        const VNCS::Real e0RestNorm;
        const VNCS::Real effStiffness;
        const VNCS::Real restArea;
        VNCS::Real restPhi = 0;
    };
    SOFA_CLASS(TriangleFabricBendingEnergy, SOFA_TEMPLATE(sofa::core::behavior::ForceField, VNCS::Space3D::VecType));

    void init();
    Real getPotentialEnergy(const sofa::core::MechanicalParams * /*mparams*/, const DataVecCoord &x) const final;
    void addForce(const sofa::core::MechanicalParams * /*mparams*/,
                  DataVecDeriv &f,
                  const DataVecCoord &x,
                  const DataVecDeriv &v) final;
    void addDForce(const sofa::core::MechanicalParams *mparams, DataVecDeriv &df, const DataVecDeriv &dx) final;

    VNCS::Real bendingStiffness() const;
    void setBendingStiffness(VNCS::Real newBendingStiffness);

    VNCS::Real beta() const;
    void setBeta(VNCS::Real beta);

    const std::filesystem::path &meshFilepath() const;
    void setMeshFilepath(const std::filesystem::path &newMeshFilepath);

    const std::vector<std::pair<Eigen::Matrix<VNCS::Real, 12, 12>, Eigen::Matrix<VNCS::Real, 12, 12>>> &hessians()
        const;

private:
    VNCS::Real computeEnergy(const TriangleFabricBendingStencil &stencil, const VecCoord &positions) const;
    VNCS::Real computeDissipationEnergy(const TriangleFabricBendingStencil &stencil,
                                        const VecCoord &positions,
                                        const VecDeriv &velocities) const;

    static VNCS::Real calculateAngleUnitVectors(const Eigen::Vector3d &v1,
                                                const Eigen::Vector3d &v2,
                                                const Eigen::Vector3d &axis);

    static VNCS::Real computeEffectiveStiffness(const VNCS::Real kBending);

    Eigen::Matrix<VNCS::Real, 12, 1> accumulate_dEdxStencil(const TriangleFabricBendingStencil &stencil,
                                                            const VecCoord &positions) const;

    Eigen::Matrix<VNCS::Real, 12, 1> accumulate_dEdvStencil(const TriangleFabricBendingStencil &stencil,
                                                            const VecCoord &positions,
                                                            const VecDeriv &velocities) const;

    Eigen::Matrix<VNCS::Real, 12, 12> accumulate_d2Edx2Stencil(const TriangleFabricBendingStencil &stencil,
                                                               const VecCoord &positions,
                                                               const VecDeriv &velocities) const;

    Eigen::Matrix<VNCS::Real, 12, 12> accumulate_d2Edv2Stencil(const TriangleFabricBendingStencil &stencil,
                                                               const VecCoord &positions,
                                                               const VecDeriv &velocities) const;

    void compute_dPsi(const Eigen::Vector3d &n1,
                      const Eigen::Vector3d &n2,
                      const Eigen::Vector3d &e0,
                      const VNCS::Real rest_phi,
                      const VNCS::Real ka,
                      VNCS::Real &dPsi_dTheta,
                      VNCS::Real &dPsi_dTheta_dTheta) const;

    void calculateAngleHessian(const std::array<Eigen::Vector3d, 2> &ns,
                               const std::array<Eigen::Vector3d, 5> &e,
                               Eigen::Matrix<VNCS::Real, 12, 12> &theta_hessian) const;

    void calculateAngleGradients(const Eigen::Vector3d &n1,
                                 const VNCS::Real n1_length,
                                 const Eigen::Vector3d &n2,
                                 const VNCS::Real n2_length,
                                 const std::array<Eigen::Vector3d, 5> &e,
                                 Eigen::Matrix<VNCS::Real, 12, 1> &thetaGradients) const;

    Eigen::Matrix<VNCS::Real, 12, 1> calculateStencilForce(const TriangleFabricBendingStencil &stencil,
                                                           const VecCoord &positions) const;
    Eigen::Matrix<VNCS::Real, 12, 1> calculateStencilDampingForce(const TriangleFabricBendingStencil &stencil,
                                                                  const VecCoord &positions,
                                                                  const VecDeriv &velocities) const;

    Eigen::Matrix<VNCS::Real, 12, 12> calculateEdgeStencilHessian(const TriangleFabricBendingStencil &stencil,
                                                                  const VecCoord &positions) const;

    Eigen::Matrix<VNCS::Real, 12, 12> calculateStencilDampingPositionHessian(
        const TriangleFabricBendingStencil &stencil,
        const VecCoord &positions,
        const VecDeriv &velocities) const;

    Eigen::Matrix<VNCS::Real, 12, 12> calculateStencilDampingVelocityHessian(
        const TriangleFabricBendingStencil &stencil,
        const VecCoord &positions,
        const VecDeriv &velocities) const;

    VNCS::Real calculateKa(const TriangleFabricBendingStencil &stencil) const;

    std::array<Eigen::Vector3d, 5> calculateEs(const TriangleFabricBendingStencil &stencil,
                                               const VecCoord &positions) const;

    void computeHinges(const VNCS::Space3D::Mesh &mesh);
    void computeRestPhis();

    std::vector<TriangleFabricBendingStencil> m_stencils;

    using Hessian = Eigen::Matrix<VNCS::Real, 12, 12>;

    VNCS::Real m_bendingStiffness;
    VNCS::Real m_dampingStiffness;
    std::filesystem::path m_meshFilepath;
    std::vector<std::pair<Hessian, Hessian>> m_hessians;
};

}  // namespace Sim3D
}  // namespace VNCS

#endif  //  VNCS_SIM3D_TRIANGLEFABRICBENDINGENERGYCOMMON_H
