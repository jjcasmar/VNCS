#ifndef CPF_TETRAMESHGENERATOR_H
#define CPF_TETRAMESHGENERATOR_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <SofaBaseTopology/TriangleSetTopologyContainer.h>

#include <CPFSofaPlugin/types.h>

namespace CPF
{
class TetraMeshGenerator : public sofa::core::objectmodel::BaseObject
{
    using TriangleSetTCLink =
        typename sofa::core::objectmodel::SingleLink<TetraMeshGenerator,
                                                     sofa::component::topology::TriangleSetTopologyContainer,
                                                     sofa::core::objectmodel::BaseLink::FLAG_STOREPATH |
                                                         sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>;

public:
    SOFA_CLASS(TetraMeshGenerator, sofa::core::objectmodel::BaseObject);
    TetraMeshGenerator();

    void init() final;

private:
    TriangleSetTCLink m_triangleTopologyLink;
    sofa::Data<sofa::helper::vector<CPF::SofaTypes::Point>> m_outputPoints;
    sofa::Data<sofa::helper::vector<CPF::SofaTypes::Tetra>> m_outputTetras;
};
}  // namespace CPF

#endif  // CPF_TETRAMESHGENERATOR_H
