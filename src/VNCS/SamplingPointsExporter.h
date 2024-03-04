#ifndef SAMPLINGPOINTSEXPORTER_H
#define SAMPLINGPOINTSEXPORTER_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/helper/vector.h>
#include <sofa/core/objectmodel/Link.h>
#include <VNCS/Types.h>
#include <VNCS/BlendingField.h>
#include <filesystem>

namespace VNCS
{
class SamplingPointsExporter : public sofa::core::objectmodel::BaseObject
{
    typedef sofa::defaulttype::Vec3Types DType;

    typedef DType::VecCoord VecCoord;
    typedef DType::Coord Coord;
    typedef DType::Real Real;
    typedef sofa::Data<VecCoord> DataVecCoord;

    using BlendingFieldLink =
        typename sofa::core::objectmodel::SingleLink<VNCS::SamplingPointsExporter,
                                                     VNCS::BlendingField,
                                                     sofa::core::objectmodel::BaseLink::FLAG_STOREPATH |
                                                         sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>;

public:
    SOFA_CLASS(SamplingPointsExporter, sofa::core::objectmodel::BaseObject);

    SamplingPointsExporter();

    void init();

    void setCoarseFilePath(const std::filesystem::path &coarseFilePath);
    void setFineFilePath(const std::filesystem::path &fineFilePath);

    void setBlendingField(std::shared_ptr<BlendingField> blendingField);
    std::shared_ptr<BlendingField> blendingField() const;

private:
    DataVecCoord m_coarsePoints;
    DataVecCoord m_finePoints;

    sofa::Data<sofa::helper::vector<VNCS::SofaTypes::Tetra>> m_lowResTetrahedra;
    sofa::Data<sofa::helper::vector<VNCS::SofaTypes::Tetra>> m_highResTetrahedra;

    std::filesystem::path m_coarseFilePath;
    std::filesystem::path m_fineFilePath;

    std::shared_ptr<BlendingField> m_blendingField;
};
}  // namespace VNCS

#endif  // SAMPLINGPOINTSEXPORTER_H
