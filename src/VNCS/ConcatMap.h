#ifndef VNCS_CONCATMAP_H
#define VNCS_CONCATMAP_H

#include <sofa/core/MultiMapping.h>
#include <sofa/core/MultiMapping.inl>
#include <VNCS/DataExtensions.h>
#include <spdlog/spdlog.h>

#include <VNCS/Spaces.h>

#include <Eigen/Dense>
#include <numeric>

namespace VNCS
{
template <typename WorldSpace>
class ConcatMap : public sofa::core::MultiMapping<typename WorldSpace::VecType, typename WorldSpace::VecType>
{
    using Inherit = sofa::core::MultiMapping<typename WorldSpace::VecType, typename WorldSpace::VecType>;
    using Real = typename WorldSpace::VecType::Real;
    using Vector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

public:
    SOFA_CLASS(SOFA_TEMPLATE(VNCS::ConcatMap, WorldSpace),  //
               SOFA_TEMPLATE2(sofa::core::MultiMapping, typename WorldSpace::VecType, typename WorldSpace::VecType));

    ConcatMap()
        : Inherit()
    {
    }

    void init() override
    {
        auto totalSize = 0;
        for (int i = 0; i < this->fromModels.size(); ++i) {
            totalSize += this->fromModels[i]->getSize();
        }
        this->toModels[0]->resize(totalSize);
        Inherit::init();
    }

    void apply(const sofa::core::MechanicalParams *mparams,
               const sofa::helper::vector<typename Inherit::OutDataVecCoord *> &dataVecOutPos,
               const sofa::helper::vector<const typename Inherit::InDataVecCoord *> &dataVecInPos) final
    {
        const auto &base = make_read_accessor(*dataVecInPos[0]);
        const auto &detail = make_read_accessor(*dataVecInPos[1]);
        auto &final = make_write_accessor(*dataVecOutPos[0]);

        auto nextChunk = std::copy(std::begin(base), std::end(base), std::begin(final));
        std::copy(std::begin(detail), std::end(detail), nextChunk);
    }

    void applyJ(const sofa::core::MechanicalParams *mparams,
                const sofa::helper::vector<typename Inherit::OutDataVecDeriv *> &dataVecOutVel,
                const sofa::helper::vector<const typename Inherit::InDataVecDeriv *> &dataVecInVel) final
    {
        const auto &base = make_read_accessor(*dataVecInVel[0]);
        const auto &detail = make_read_accessor(*dataVecInVel[1]);
        auto &final = make_write_accessor(*dataVecOutVel[0]);

        auto nextChunk = std::copy(std::begin(base), std::end(base), std::begin(final));
        std::copy(std::begin(detail), std::end(detail), nextChunk);
    }

    void applyJT(const sofa::core::MechanicalParams *mparams,
                 const sofa::helper::vector<typename Inherit::InDataVecDeriv *> &dataVecOutForce,
                 const sofa::helper::vector<const typename Inherit::InDataVecDeriv *> &dataVecInForce) final
    {
        auto &base = make_write_accessor(*dataVecOutForce[0]);
        auto &detail = make_write_accessor(*dataVecOutForce[1]);
        const auto &final = make_read_accessor(*dataVecInForce[0]);

        std::transform(std::begin(final),  //
                       std::begin(final) + std::distance(std::begin(base), std::end(base)),
                       std::begin(base),
                       std::begin(base),
                       std::plus<>());

        std::transform(std::begin(final) + std::distance(std::begin(base), std::end(base)),  //
                       std::end(final),
                       std::begin(detail),
                       std::begin(detail),
                       std::plus<>());
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
};

namespace Sim1D
{
using ConcatMap = VNCS::ConcatMap<VNCS::Space1D>;
}
namespace Sim2D
{
using ConcatMap = VNCS::ConcatMap<VNCS::Space2D>;
}
namespace Sim3D
{
using ConcatMap = VNCS::ConcatMap<VNCS::Space3D>;
}
}  // namespace VNCS

#endif  // VNCS_CONCATMAP_H
