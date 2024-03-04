#ifndef VNCS_PROJECTION_H
#define VNCS_PROJECTION_H

#include <sofa/core/behavior/ProjectiveConstraintSet.h>
#include <chrono>
#include <ratio>
#include <sofa/core/behavior/ProjectiveConstraintSet.inl>

#include <VNCS/DataExtensions.h>
#include <VNCS/Spaces.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <spdlog/spdlog.h>

namespace VNCS
{
template <typename WorldSpace>
class Projection : public sofa::core::behavior::ProjectiveConstraintSet<typename WorldSpace::VecType>
{
    using Inherit = typename sofa::core::behavior::ProjectiveConstraintSet<typename WorldSpace::VecType>;

public:
    enum class ProjectionMode { Barycentric, Cluster };

    using Real = typename Inherit::Real;
    SOFA_CLASS(SOFA_TEMPLATE(VNCS::Projection, WorldSpace),  //
               SOFA_TEMPLATE(sofa::core::behavior::ProjectiveConstraintSet, typename WorldSpace::VecType));

    Projection()
        : Inherit()
    {
    }

    void init() { Inherit::init(); }

    void projectPosition(const sofa::core::MechanicalParams *mparams, typename Inherit::DataVecCoord &x) final
    {
        apply(x);
    }

    void projectVelocity(const sofa::core::MechanicalParams *mparams, typename Inherit::DataVecDeriv &v) final
    {
        apply(v);
    }

    void projectResponse(const sofa::core::MechanicalParams *mparams, typename Inherit::DataVecDeriv &dx) final
    {
        apply(dx);
    }

    void projectJacobianMatrix(const sofa::core::MechanicalParams * /*mparams*/,
                               typename Inherit::DataMatrixDeriv &cData) final
    {
    }

    const Eigen::SparseMatrix<Real> &barycentricMatrix() const { return m_barycentricMatrix; }

    void setBarycentricMatrix(Eigen::SparseMatrix<Real> barycentricMatrix)
    {
        m_barycentricMatrix = std::move(barycentricMatrix);
        Eigen::SparseMatrix<Real> BtB = (m_barycentricMatrix).transpose() * (m_barycentricMatrix);
        m_BtBLDLT.compute(BtB);

        m_projectionMode = ProjectionMode::Barycentric;
    }

    const Eigen::SparseMatrix<Real> &clusterMatrix() const { return m_clusterMatrix; }

    void setClusterMatrix(Eigen::SparseMatrix<Real> clusterMatrix)
    {
        m_clusterMatrix = std::move(clusterMatrix);
        Eigen::SparseMatrix<Real> CtC = (m_clusterMatrix).transpose() * (m_clusterMatrix);
        m_CtCInverse.resize(CtC.rows());
        for (int i = 0; i < CtC.rows(); ++i) {
            m_CtCInverse.diagonal()[i] = 1.0 / CtC.diagonal()[i];
        }

        m_projectionMode = ProjectionMode::Cluster;
    }

    const std::vector<int> &durations() const { return m_durations; }
    void clearDurations() { m_durations.clear(); }

private:
    template <class D>
    void apply(D &v)
    {
        auto start = std::chrono::steady_clock::now();
        // auto &d = make_write_accessor(v);
        // Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>> dxEig(&(d[0][0]), WorldSpace::dim * d.size());

        // const auto &C = m_clusterMatrix;
        // dxEig -= C * (m_CtCInverse * (C.transpose() * dxEig));
        switch (m_projectionMode) {
            case ProjectionMode::Cluster: {
                auto &d = make_write_accessor(v);
                Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>> dxEig(&(d[0][0]), WorldSpace::dim * d.size());

                const auto &C = m_clusterMatrix;
                dxEig -= C * (m_CtCInverse * (C.transpose() * dxEig));
                break;
            }

            case ProjectionMode::Barycentric: {
                auto &d = make_write_accessor(v);
                Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>> dxEig(&(d[0][0]), WorldSpace::dim * d.size());

                const auto &C = m_barycentricMatrix;
                dxEig -= C * (m_BtBLDLT.solve(C.transpose() * dxEig));
                break;
            }
        }

        auto end = std::chrono::steady_clock::now();
        m_durations.emplace_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    }

    ProjectionMode m_projectionMode;

    Eigen::SparseMatrix<Real> m_clusterMatrix;
    Eigen::DiagonalMatrix<Real, Eigen::Dynamic> m_CtCInverse;

    Eigen::SparseMatrix<Real> m_barycentricMatrix;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>> m_BtBLDLT;

    std::vector<int> m_durations;
};

namespace Sim1D
{
using Projection = VNCS::Projection<VNCS::Space1D>;
}  // namespace Sim1D

namespace Sim2D
{
using Projection = VNCS::Projection<VNCS::Space2D>;
}  // namespace Sim2D

namespace Sim3D
{
using Projection = VNCS::Projection<VNCS::Space3D>;
}  // namespace Sim3D

}  // namespace VNCS

#endif  // VNCS_PROJECTION_H
