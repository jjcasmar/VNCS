#ifndef VNCS_VEC2D23D_H
#define VNCS_VEC2D23D_H

#include <sofa/core/Mapping.h>
#include <sofa/core/Mapping.inl>

#include <VNCS/DataExtensions.h>
#include <spdlog/spdlog.h>

namespace VNCS
{
namespace Sim2D
{
class Vec2D23D : public sofa::core::Mapping<sofa::defaulttype::Vec2Types, sofa::defaulttype::Vec3Types>
{
public:
    using Inherit = sofa::core::Mapping<sofa::defaulttype::Vec2dTypes, sofa::defaulttype::Vec3Types>;
    using Real = sofa::defaulttype::Vec2Types::Real;

public:
    SOFA_CLASS(VNCS::Sim2D::Vec2D23D,  //
               SOFA_TEMPLATE2(sofa::core::Mapping, sofa::defaulttype::Vec2Types, sofa::defaulttype::Vec3Types));

    Vec2D23D();

    void init();

    void apply(const sofa::core::MechanicalParams *mparams,
               Inherit::OutDataVecCoord &dataOutPos,
               const Inherit::InDataVecCoord &dataInPos) final;

    void applyJ(const sofa::core::MechanicalParams *mparams,
                Inherit::OutDataVecDeriv &dataOutPos,
                const Inherit::InDataVecDeriv &dataInPos) final;

    void applyJT(const sofa::core::MechanicalParams *mparams,
                 Inherit::InDataVecDeriv &dataOutPos,
                 const Inherit::OutDataVecDeriv &dataInPos) final;

private:
};
}  // namespace Sim2D
}  // namespace VNCS

#endif  // VNCS_DEFORMATIONGRADIENTMAP_H
