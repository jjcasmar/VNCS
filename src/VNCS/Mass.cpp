#include "Mass.h"

VNCS::Mass::Mass()
    : Inherit()
    , m_samplingPointsLink(initLink("samplingPoints", "Sampling points"))
    , m_density(initData(&m_density, "density", "Density for this mass"))
{
    m_logger = spdlog::get("VNCS");
}

void VNCS::Mass::init()
{
    m_logger->trace("VNCS::Mass::init");

    m_samplingPoints = m_samplingPointsLink.get();
}

void VNCS::Mass::addMDx(const sofa::core::MechanicalParams *mparams,
                        Inherit::DataVecDeriv &f,
                        const Inherit::Mass::DataVecDeriv &dx,
                        SReal factor)
{
}

void VNCS::Mass::addForce(const sofa::core::MechanicalParams *,
                          Inherit::DataVecDeriv &f,
                          const Inherit::DataVecCoord &x,
                          const Inherit::DataVecDeriv &v)
{
}


namespace VNCS {
Real Mass::density() const
{
    return m_density;
}

void Mass::setDensity(const Real &density)
{
    m_density = density;
}

std::shared_ptr<VNCS::SamplingPoints<T> > Mass::samplingPoints() const
{
    return m_samplingPoints;
}

void Mass::setSamplingPoints(const std::shared_ptr<VNCS::SamplingPoints<T> > &samplingPoints)
{
    m_samplingPoints = samplingPoints;
}

}
