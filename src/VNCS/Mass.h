#ifndef VNCS_MASS_H
#define VNCS_MASS_H

#include <sofa/core/behavior/Mass.h>
#include <sofa/core/behavior/Mass.inl>

#include <VNCS/SamplingPoints.h>
#include <VNCS/DataExtensions.h>
#include <VNCS/Spaces.h>

#include <memory>

namespace VNCS
{
template <typename WorldSpace, typename MaterialSpace>
class Mass : public sofa::core::behavior::Mass<typename WorldSpace::VecType>
{
    using Inherit = sofa::core::behavior::Mass<typename WorldSpace::VecType>;
    using ForceField = sofa::core::behavior::ForceField<typename WorldSpace::VecType>;

public:
    SOFA_CLASS(SOFA_TEMPLATE2(Mass, WorldSpace, MaterialSpace),  //
               SOFA_TEMPLATE(sofa::core::behavior::Mass, typename WorldSpace::VecType));
    Mass()
        : Inherit()
    {
    }

    void addMDx(const sofa::core::MechanicalParams *mparams,
                typename Inherit::DataVecDeriv &f,
                const typename Inherit::DataVecDeriv &dx,
                typename Inherit::Real factor)
    {
        auto &f_ = make_write_accessor(f);
        const auto &dx_ = make_read_accessor(dx);

        Eigen::Map<Eigen::VectorXd> fEig(&(f_[0][0]), WorldSpace::dim * f_.size());
        const auto &samplers = *m_samplingPoints;
        for (int i = 0; i < samplers.size(); ++i)
            f_[i] += factor * samplers[i].w * m_density * dx_[i];
    }

    void addForce(const sofa::core::MechanicalParams * /*mparams*/,
                  typename Inherit::DataVecDeriv &f,
                  const typename Inherit::DataVecCoord &x,
                  const typename Inherit::DataVecDeriv &v) final
    {
        auto &f_ = make_write_accessor(f);
        const auto gravity = typename WorldSpace::VecType::Coord(this->getContext()->getGravity());

        const auto &samplers = *m_samplingPoints;
        for (int i = 0; i < samplers.size(); ++i)
            f_[i] += samplers[i].w * m_density * gravity;
    }

    typename WorldSpace::VecType::Real density() const { return m_density; }
    void setDensity(const typename WorldSpace::VecType::Real &density) { m_density = density; }

    VNCS::Real getKineticEnergy(const sofa::core::MechanicalParams *,
                                const typename Inherit::DataVecDeriv &v) const final
    {
        const auto &v_ = make_read_accessor(v);

        const auto &samplers = *m_samplingPoints;
        VNCS::Real energy = 0;
        for (int i = 0; i < samplers.size(); ++i)

            energy += 0.5 * v_[i] * (samplers[i].w * m_density * v_[i]);
        return energy;
    }

    VNCS::Real getPotentialEnergy(const sofa::core::MechanicalParams *mparams,
                                  const typename Inherit::DataVecCoord &xVec) const final
    {
        const auto &msState = static_cast<const ForceField *>(this)->getMState();

        const auto &u = make_read_accessor(xVec);
        const auto &x = msState->readRestPositions();
        const auto gravity = typename WorldSpace::VecType::Coord(this->getContext()->getGravity());

        const auto &samplers = *m_samplingPoints;
        VNCS::Real energy = 0;
        for (int i = 0; i < samplers.size(); ++i)
            energy += samplers[i].w * m_density * (-1.0 * gravity) * (x[i] + u[i]);
        return energy;
    }

    std::shared_ptr<VNCS::SamplingPoints<MaterialSpace>> samplingPoints() const { return m_samplingPoints; }
    void setSamplingPoints(const std::shared_ptr<VNCS::SamplingPoints<MaterialSpace>> &samplingPoints)
    {
        m_samplingPoints = samplingPoints;
    }

private:
    std::shared_ptr<VNCS::SamplingPoints<MaterialSpace>> m_samplingPoints;
    typename WorldSpace::VecType::Real m_density;
};

namespace Sim1D
{
using Mass11 = VNCS::Mass<VNCS::Space1D, VNCS::Space1D>;
}

namespace Sim2D
{
using Mass21 = VNCS::Mass<VNCS::Space2D, VNCS::Space1D>;
using Mass22 = VNCS::Mass<VNCS::Space2D, VNCS::Space2D>;
}  // namespace Sim2D

namespace Sim3D
{
using Mass31 = VNCS::Mass<VNCS::Space3D, VNCS::Space1D>;
using Mass32 = VNCS::Mass<VNCS::Space3D, VNCS::Space2D>;
using Mass33 = VNCS::Mass<VNCS::Space3D, VNCS::Space3D>;
}  // namespace Sim3D

}  // namespace VNCS

#endif  // VNCS_MASS_H
