#include "Vec2D23D.h"

VNCS::Sim2D::Vec2D23D::Vec2D23D()
{
}

void VNCS::Sim2D::Vec2D23D::init()
{
    if (this->toModel && this->fromModel)
        this->toModel->resize(this->fromModel->getSize());

    Inherit::init();
}

void VNCS::Sim2D::Vec2D23D::apply(const sofa::core::MechanicalParams *mparams,
                                  Inherit::OutDataVecCoord &dataOutPos,
                                  const Inherit::InDataVecCoord &dataInPos)
{
    if (this->toModel && this->fromModel) {
        auto &outPos = make_write_accessor(dataOutPos);
        const auto &inPos = make_read_accessor(dataInPos);

        for (int i = 0; i < inPos.size(); ++i) {
            outPos[i][0] = inPos[i][0];
            outPos[i][1] = inPos[i][1];
            outPos[i][2] = 0;
        }
    }
}

void VNCS::Sim2D::Vec2D23D::applyJ(const sofa::core::MechanicalParams *mparams,
                                   Inherit::OutDataVecDeriv &dataOutPos,
                                   const Inherit::InDataVecDeriv &dataInPos)
{
    if (this->toModel && this->fromModel) {
        auto &outPos = make_write_accessor(dataOutPos);
        const auto &inPos = make_read_accessor(dataInPos);

        for (int i = 0; i < inPos.size(); ++i) {
            outPos[i][0] = inPos[i][0];
            outPos[i][1] = inPos[i][1];
            outPos[i][2] = 0;
        }
    }
}

void VNCS::Sim2D::Vec2D23D::applyJT(const sofa::core::MechanicalParams *mparams,
                                    Inherit::InDataVecDeriv &dataOutPos,
                                    const Inherit::OutDataVecDeriv &dataInPos)
{
    if (this->toModel && this->fromModel) {
        auto &outPos = make_write_accessor(dataOutPos);
        const auto &inPos = make_read_accessor(dataInPos);

        for (int i = 0; i < inPos.size(); ++i) {
            outPos[i][0] += inPos[i][0];
            outPos[i][1] += inPos[i][1];
        }
    }
}
