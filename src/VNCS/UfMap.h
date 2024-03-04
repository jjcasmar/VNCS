#ifndef VNCS_UFMAP_H
#define VNCS_UFMAP_H

#include <VNCS/DataExtensions.h>
#include <sofa/core/MultiMapping.h>
#include <sofa/core/MultiMapping.inl>
#include <spdlog/spdlog.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>

#include <VNCS/Spaces.h>
#include <unsupported/Eigen/src/SparseExtra/MarketIO.h>

#include <numeric>

namespace VNCS
{
// uf = C * uc + P * duf

template <typename WorldSpace>
class UfMap : public sofa::core::MultiMapping<typename WorldSpace::VecType, typename WorldSpace::VecType>
{
    using Inherit = sofa::core::MultiMapping<typename WorldSpace::VecType, typename WorldSpace::VecType>;

public:
    using Real = typename WorldSpace::VecType::Real;
    using Vector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
    SOFA_CLASS(SOFA_TEMPLATE(VNCS::UfMap, WorldSpace),  //
               SOFA_TEMPLATE2(sofa::core::MultiMapping, typename WorldSpace::VecType, typename WorldSpace::VecType));

    UfMap()
        : Inherit()
    {
    }

    void init() override
    {
        const auto &C = m_clusterMatrix;
        this->toModels[0]->resize(C.rows() / WorldSpace::dim);
        Inherit::init();
    }

    void apply(const sofa::core::MechanicalParams *mparams,
               const sofa::helper::vector<typename Inherit::OutDataVecCoord *> &dataVecOutPos,
               const sofa::helper::vector<const typename Inherit::InDataVecCoord *> &dataVecInPos) final
    {
        const auto &C = m_clusterMatrix;
        if (C.rows()) {
            const auto &base = make_read_accessor(*dataVecInPos[0]);
            const auto &detail = make_read_accessor(*dataVecInPos[1]);
            auto &final = make_write_accessor(*dataVecOutPos[0]);

            Eigen::Map<const Vector> baseEig(&(base[0][0]), WorldSpace::dim * base.size());
            Eigen::Map<const Vector> detailEig(&(detail[0][0]), WorldSpace::dim * detail.size());
            Eigen::Map<Vector> finalEig(&(final[0][0]), WorldSpace::dim * final.size());

            finalEig = C * baseEig + detailEig;
        }
    }

    void applyJ(const sofa::core::MechanicalParams *mparams,
                const sofa::helper::vector<typename Inherit::OutDataVecDeriv *> &dataVecOutVel,
                const sofa::helper::vector<const typename Inherit::InDataVecDeriv *> &dataVecInVel) final
    {
        const auto &C = m_clusterMatrix;
        if (C.rows()) {
            const auto &base = make_read_accessor(*dataVecInVel[0]);
            const auto &detail = make_read_accessor(*dataVecInVel[1]);
            auto &final = make_write_accessor(*dataVecOutVel[0]);

            Eigen::Map<const Vector> baseEig(&(base[0][0]), WorldSpace::dim * base.size());
            Eigen::Map<const Vector> detailEig(&(detail[0][0]), WorldSpace::dim * detail.size());
            Eigen::Map<Vector> finalEig(&(final[0][0]), WorldSpace::dim * final.size());

            finalEig = C * baseEig + detailEig;
        }
    }

    void applyJT(const sofa::core::MechanicalParams *mparams,
                 const sofa::helper::vector<typename Inherit::OutDataVecDeriv *> &dataVecOutForce,
                 const sofa::helper::vector<const typename Inherit::InDataVecDeriv *> &dataVecInForce) final
    {
        const auto &C = m_clusterMatrix;
        if (C.rows()) {
            auto &base = make_write_accessor(*dataVecOutForce[0]);
            auto &detail = make_write_accessor(*dataVecOutForce[1]);
            const auto &final = make_read_accessor(*dataVecInForce[0]);

            Eigen::Map<Vector> baseEig(&(base[0][0]), WorldSpace::dim * base.size());
            Eigen::Map<Vector> detailEig(&(detail[0][0]), WorldSpace::dim * detail.size());
            Eigen::Map<const Vector> finalEig(&(final[0][0]), WorldSpace::dim * final.size());

            baseEig += C.transpose() * finalEig;
            detailEig += finalEig;
        }
    }

    void applyDJT(const sofa::core::MechanicalParams * /*mparams*/,
                  sofa::core::MultiVecDerivId /*inForce*/,
                  sofa::core::ConstMultiVecDerivId /*outForce*/) final
    {
    }

    void applyJT(const sofa::core::ConstraintParams * /* cparams */,
                 const sofa::helper::vector<typename Inherit::InDataMatrixDeriv *> & /* dataMatOutConst */,
                 const sofa::helper::vector<const typename Inherit::OutDataMatrixDeriv *> & /* dataMatInConst */) final
    {
    }

    std::shared_ptr<const Eigen::SparseMatrix<Real>> clusterMatrix() const { return m_clusterMatrix; }

    void setClusterMatrix(const Eigen::SparseMatrix<Real> &clusterMatrix) { m_clusterMatrix = clusterMatrix; }

private:
    Eigen::SparseMatrix<Real> m_clusterMatrix;
};

namespace Sim1D
{
using UfMap = VNCS::UfMap<VNCS::Space1D>;
}

namespace Sim2D
{
using UfMap = VNCS::UfMap<VNCS::Space2D>;
}  // namespace Sim2D

namespace Sim3D
{
using UfMap = VNCS::UfMap<VNCS::Space3D>;
}  // namespace Sim3D
}  // namespace VNCS

#endif  //  VNCS_UFMAP_H
