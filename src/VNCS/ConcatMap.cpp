#include "ConcatMap.h"

void VNCS::ConcatMap::init()
{
}

void VNCS::ConcatMap::apply(const sofa::core::MechanicalParams *mparams,
                            const sofa::helper::vector<VNCS::ConcatMap::DataVecCoord *> &dataVecOutPos,
                            const sofa::helper::vector<const VNCS::ConcatMap::DataVecCoord *> &dataVecInPos)
{
}

void VNCS::ConcatMap::applyJ(const sofa::core::MechanicalParams *mparams,
                             const sofa::helper::vector<VNCS::ConcatMap::DataVecDeriv *> &dataVecOutVel,
                             const sofa::helper::vector<const VNCS::ConcatMap::DataVecDeriv *> &dataVecInVel)
{
    m_logger->trace("VNCS::ConcatMap::applyJ");
    const VecCoord &base = make_read_accessor(*dataVecInVel[0]);
    const VecCoord &detail = make_read_accessor(*dataVecInVel[1]);
    VecCoord &final = make_write_accessor(*dataVecOutVel[0]);

    auto nextChunk = std::copy(std::begin(base), std::end(base), std::begin(final));
    std::copy(std::begin(detail), std::end(detail), nextChunk);
}

void VNCS::ConcatMap::applyJT(const sofa::core::MechanicalParams *mparams,
                              const sofa::helper::vector<VNCS::ConcatMap::DataVecDeriv *> &dataVecOutForce,
                              const sofa::helper::vector<const VNCS::ConcatMap::DataVecDeriv *> &dataVecInForce)
{
}

void VNCS::ConcatMap::applyDJT(const sofa::core::MechanicalParams *,
                               sofa::core::MultiVecDerivId,
                               sofa::core::ConstMultiVecDerivId)
{
}

void VNCS::ConcatMap::applyJT(const sofa::core::ConstraintParams *cparams,
                              const sofa::helper::vector<VNCS::ConcatMap::DataMatrixDeriv *> &dataMatOutConst,
                              const sofa::helper::vector<const VNCS::ConcatMap::DataMatrixDeriv *> &dataMatInConst)
{
}
