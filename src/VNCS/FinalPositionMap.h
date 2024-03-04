#ifndef VNCS_FINALPOSITIONMAP_H
#define VNCS_FINALPOSITIONMAP_H

#include <sofa/core/Mapping.h>
#include <sofa/core/Mapping.inl>

#include <VNCS/DataExtensions.h>
#include <spdlog/spdlog.h>

namespace VNCS
{
template <typename WorldSpace>
class FinalPositionMap : public sofa::core::Mapping<typename WorldSpace::VecType, typename WorldSpace::VecType>
{
public:
    using Inherit = sofa::core::Mapping<typename WorldSpace::VecType, typename WorldSpace::VecType>;
    using Real = typename WorldSpace::Real;
    using Vector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

public:
    SOFA_CLASS(SOFA_TEMPLATE(VNCS::FinalPositionMap, WorldSpace),  //
               SOFA_TEMPLATE2(sofa::core::Mapping, typename WorldSpace::VecType, typename WorldSpace::VecType));

    FinalPositionMap() = default;

    void init()
    {
        if (this->toModel && this->fromModel)
            this->toModel->resize(this->fromModel->getSize());

        Inherit::init();
    }

    void apply(const sofa::core::MechanicalParams *mparams,
               typename Inherit::OutDataVecCoord &dataOutPos,
               const typename Inherit::InDataVecCoord &dataInPos) final
    {
        if (this->toModel && this->fromModel) {
            // Get rest position
            const auto &restPositions = this->fromModel->readRestPositions();
            auto &outPos = make_write_accessor(dataOutPos);
            const auto &inPos = make_read_accessor(dataInPos);
            Eigen::Map<const Vector> rest((&restPositions[0][0]), WorldSpace::dim * restPositions.size());
            Eigen::Map<Vector> out(&(outPos[0][0]), WorldSpace::dim * outPos.size());
            Eigen::Map<const Vector> in(&(inPos[0][0]), WorldSpace::dim * inPos.size());

            out = rest + in;
        }
    }

    void applyJ(const sofa::core::MechanicalParams *mparams,
                typename Inherit::OutDataVecDeriv &dataOutVel,
                const typename Inherit::InDataVecDeriv &dataInVel) final
    {
        if (this->toModel && this->fromModel) {
            auto &outVel = make_write_accessor(dataOutVel);
            const auto &inVel = make_read_accessor(dataInVel);

            std::copy(std::begin(inVel), std::end(inVel), std::begin(outVel));
        }
    }

    void applyJT(const sofa::core::MechanicalParams *mparams,
                 typename Inherit::InDataVecDeriv &dataOutForce,
                 const typename Inherit::OutDataVecDeriv &dataInForce) final
    {
        if (this->toModel && this->fromModel) {
            auto &outForce = make_write_accessor(dataOutForce);
            const auto &inForce = make_read_accessor(dataInForce);

            std::transform(std::begin(inForce),  //
                           std::end(inForce),
                           std::begin(outForce),
                           std::begin(outForce),
                           std::plus<>{});
        }
    }
};
}  // namespace VNCS

#endif  // VNCS_FINALPOSITIONMAP_H
