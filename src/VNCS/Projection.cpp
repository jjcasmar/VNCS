#include "Projection.h"

VNCS::Projection::Projection(sofa::core::behavior::MechanicalState<VNCS::Projection::DataType> *ms)
    : Inherit(ms)
    , m_C(initData(&m_C, Eigen::SparseMatrix<Real>{}, "C", "Cluster matrix"))
{
}

void VNCS::Projection::init()
{
    m_logger->trace("VNCS::Projection::init");
    const auto &C = make_read_accessor(m_C);
    Eigen::SparseMatrix<Real> CtC = (C.transpose() * C);
    m_CtCInverse.resize(CtC.rows());
    for (int i = 0; i < CtC.rows(); ++i) {
        if (CtC.coeff(i, i) == Real{0.0})
            m_logger->error("Coarse node {} doesnt have associated d nodes", i);
        m_CtCInverse.diagonal()[i] = 1.0 / CtC.coeff(i, i);
    }

    Inherit::init();
}

void VNCS::Projection::projectResponse(const sofa::core::MechanicalParams *mparams, VNCS::Projection::DataVecDeriv &dx)
{
    m_logger->trace("VNCS::Projection::projectResponse");
    apply(dx);
}

void VNCS::Projection::projectVelocity(const sofa::core::MechanicalParams *mparams, VNCS::Projection::DataVecDeriv &v)
{
    m_logger->trace("VNCS::Projection::projectVelocity");
    apply(v);
}

void VNCS::Projection::projectPosition(const sofa::core::MechanicalParams *mparams, VNCS::Projection::DataVecCoord &x)
{
    m_logger->trace("VNCS::Projection::projectPosition");
    apply(x);
}

void VNCS::Projection::projectJacobianMatrix(const sofa::core::MechanicalParams *,
                                             VNCS::Projection::DataMatrixDeriv &cData)
{
    m_logger->warn(
        "VNCS::Projection::projectJacobianMatrix not implement\nThis projection doesn't support LM constraints");
}
