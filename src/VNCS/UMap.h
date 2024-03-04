#ifndef VNCS_UMAP_H
#define VNCS_UMAP_H

#include <sofa/core/MultiMapping.h>
#include <sofa/core/MultiMapping.inl>

#include <VNCS/SamplingPoints.h>
#include <VNCS/DataExtensions.h>
#include <VNCS/Spaces.h>
#include <unsupported/Eigen/src/SparseExtra/MarketIO.h>

#include <Eigen/Sparse>
#include <type_traits>
#include <sofa/core/VecId.h>

namespace VNCS
{
template <typename WorldSpace, typename MaterialSpace>
class UMap : public sofa::core::MultiMapping<typename WorldSpace::VecType, typename WorldSpace::VecType>
{
    using Inherit = sofa::core::MultiMapping<typename WorldSpace::VecType, typename WorldSpace::VecType>;
    using Real = typename WorldSpace::VecType::Real;
    using Vector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

public:
    SOFA_CLASS(SOFA_TEMPLATE2(VNCS::UMap, WorldSpace, MaterialSpace),  //
               SOFA_TEMPLATE2(sofa::core::MultiMapping, typename WorldSpace::VecType, typename WorldSpace::VecType));

    UMap()
        : Inherit()
    {
    }

    void init()
    {
        if (m_coarseInput) {
            this->addInputModel(m_coarseInput);
            m_coarseIndex = 0;
        }

        if (m_fineInput) {
            this->addInputModel(m_fineInput);
            m_fineIndex = m_coarseIndex.has_value() ? 1 : 0;
        }

        if (m_output)
            this->addOutputModel(m_output);

        const auto &samplers = *m_samplingPoints;
        this->toModels[0]->resize(samplers.size());

        std::vector<Eigen::Triplet<Real>> phiCTriplets;
        std::vector<Eigen::Triplet<Real>> phiFTriplets;

        for (int i = 0; i < samplers.size(); ++i) {
            const auto &sampler = samplers[i];
            const auto a = sampler.a;
            const auto &coarseShapeFunctions = sampler.coarseShapeFunctions;
            for (int j = 0; j < coarseShapeFunctions.size(); ++j) {
                const auto &shapeFunction = coarseShapeFunctions[j];

                for (int k = 0; k < WorldSpace::dim; ++k)
                    phiCTriplets.emplace_back(
                        WorldSpace::dim * i + k, WorldSpace::dim * shapeFunction.nodeIndex + k, a * shapeFunction.v);
            }

            const auto &fineShapeFunctions = sampler.fineShapeFunctions;
            for (int j = 0; j < fineShapeFunctions.size(); ++j) {
                const auto &shapeFunction = fineShapeFunctions[j];

                for (int k = 0; k < WorldSpace::dim; ++k)
                    phiFTriplets.emplace_back(WorldSpace::dim * i + k,
                                              WorldSpace::dim * shapeFunction.nodeIndex + k,
                                              (1 - a) * shapeFunction.v);
            }
        }

        if (m_coarseInput) {
            m_phiC.resize(WorldSpace::dim * samplers.size(),
                          WorldSpace::dim * this->fromModels[m_coarseIndex.value()]->getSize());
            m_phiC.setFromTriplets(std::begin(phiCTriplets), std::end(phiCTriplets));
        }

        if (m_fineInput) {
            m_phiF.resize(WorldSpace::dim * samplers.size(),
                          WorldSpace::dim * this->fromModels[m_fineIndex.value()]->getSize());
            m_phiF.setFromTriplets(std::begin(phiFTriplets), std::end(phiFTriplets));
        }

        Inherit::init();
    }

    void apply(const sofa::core::MechanicalParams *mparams,
               const sofa::helper::vector<typename Inherit::OutDataVecCoord *> &dataVecOutPos,
               const sofa::helper::vector<const typename Inherit::InDataVecCoord *> &dataVecInPos) final
    {
        auto &deformationAtSampler = make_write_accessor(*dataVecOutPos[0]);

        Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>> u(&(deformationAtSampler[0][0]),  //
                                                             WorldSpace::dim * deformationAtSampler.size());
        u.setZero();

        if (m_coarseInput) {
            const auto &coarseDeformationAtNodes = make_read_accessor(*dataVecInPos[m_coarseIndex.value()]);
            if (coarseDeformationAtNodes.size()) {
                Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> uc(
                    &(coarseDeformationAtNodes[0][0]),  //
                    WorldSpace::dim * coarseDeformationAtNodes.size());
                u += m_phiC * uc;
            }
        }

        if (m_fineInput) {
            const auto &fineDeformationAtNodes = make_read_accessor(*dataVecInPos[m_fineIndex.value()]);
            if (fineDeformationAtNodes.size()) {
                Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> uf(
                    &(fineDeformationAtNodes[0][0]),  //
                    WorldSpace::dim * fineDeformationAtNodes.size());
                u += m_phiF * uf;
            }
        }
    }

    void applyJ(const sofa::core::MechanicalParams *mparams,
                const sofa::helper::vector<typename Inherit::OutDataVecDeriv *> &dataVecOutVel,
                const sofa::helper::vector<const typename Inherit::InDataVecDeriv *> &dataVecInVel) final
    {
        auto &deformationAtSampler = make_write_accessor(*dataVecOutVel[0]);

        Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>> u(&(deformationAtSampler[0][0]),  //
                                                             WorldSpace::dim * deformationAtSampler.size());
        u.setZero();

        if (m_coarseInput) {
            const auto &coarseDeformationAtNodes = make_read_accessor(*dataVecInVel[m_coarseIndex.value()]);
            if (coarseDeformationAtNodes.size()) {
                Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> uc(
                    &(coarseDeformationAtNodes[0][0]),  //
                    WorldSpace::dim * coarseDeformationAtNodes.size());
                u += m_phiC * uc;
            }
        }

        if (m_fineInput) {
            const auto &fineDeformationAtNodes = make_read_accessor(*dataVecInVel[m_fineIndex.value()]);
            if (fineDeformationAtNodes.size()) {
                Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> uf(
                    &(fineDeformationAtNodes[0][0]),  //
                    WorldSpace::dim * fineDeformationAtNodes.size());
                u += m_phiF * uf;
            }
        }
    }

    void applyJT(const sofa::core::MechanicalParams *mparams,
                 const sofa::helper::vector<typename Inherit::OutDataVecDeriv *> &dataVecOutForce,
                 const sofa::helper::vector<const typename Inherit::InDataVecDeriv *> &dataVecInForce) final
    {
        const auto &samplerForce = make_read_accessor(*dataVecInForce[0]);

        Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> fSampler(&(samplerForce[0][0]),  //
                                                                          WorldSpace::dim * samplerForce.size());

        if (m_coarseInput) {
            auto &coarseForce = make_write_accessor(*dataVecOutForce[m_coarseIndex.value()]);

            if (coarseForce.size()) {
                Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>> fCoarse(&(coarseForce[0][0]),  //
                                                                           WorldSpace::dim * coarseForce.size());
                fCoarse += m_phiC.transpose() * fSampler;
            }
        }

        if (m_fineInput) {
            auto &fineForce = make_write_accessor(*dataVecOutForce[m_fineIndex.value()]);
            if (fineForce.size()) {
                Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>> fFine(&(fineForce[0][0]),  //
                                                                         WorldSpace::dim * fineForce.size());

                fFine += m_phiF.transpose() * fSampler;
            }
        }
    }

    void applyDJT(const sofa::core::MechanicalParams * /*mparams*/,
                  sofa::core::MultiVecDerivId /*inForce*/,
                  sofa::core::ConstMultiVecDerivId /*outForce*/) final
    {
    }

    // This is for mapping constraints. We dont need it
    void applyJT(const sofa::core::ConstraintParams *cparams,
                 const sofa::helper::vector<typename Inherit::InDataMatrixDeriv *> &dataMatOutConst,
                 const sofa::helper::vector<const typename Inherit::OutDataMatrixDeriv *> &dataMatInConst) final
    {
    }

    void setSamplingPoints(const std::shared_ptr<SamplingPoints<MaterialSpace>> &samplingPoints)
    {
        m_samplingPoints = samplingPoints;
    }

    const Eigen::SparseMatrix<Real> &phiC() const { return m_phiC; }
    const Eigen::SparseMatrix<Real> &phiF() const { return m_phiC; }

    sofa::core::objectmodel::BaseObject *coarseInput() const { return m_coarseInput; }
    void setCoarseInput(sofa::core::objectmodel::BaseObject *newCoarseInput)
    {
        m_coarseInput = dynamic_cast<sofa::core::BaseState *>(newCoarseInput);
    }

    sofa::core::objectmodel::BaseObject *fineInput() const { return m_fineInput; }
    void setFineInput(sofa::core::objectmodel::BaseObject *newFineInput)
    {
        m_fineInput = dynamic_cast<sofa::core::BaseState *>(newFineInput);
    }

    sofa::core::objectmodel::BaseObject *output() const { return m_output; }
    void setOutput(sofa::core::objectmodel::BaseObject *newOutput)
    {
        m_output = dynamic_cast<sofa::core::BaseState *>(newOutput);
    }

private:
    Eigen::SparseMatrix<Real> m_phiC;
    Eigen::SparseMatrix<Real> m_phiF;

    sofa::core::BaseState *m_output = nullptr;
    sofa::core::BaseState *m_coarseInput = nullptr;
    sofa::core::BaseState *m_fineInput = nullptr;

    std::optional<int> m_coarseIndex;
    std::optional<int> m_fineIndex;

    std::shared_ptr<SamplingPoints<MaterialSpace>> m_samplingPoints;
};

namespace Sim1D
{
using UMap11 = VNCS::UMap<VNCS::Space1D, VNCS::Space1D>;
}  // namespace Sim1D

namespace Sim2D
{
using UMap21 = VNCS::UMap<VNCS::Space2D, VNCS::Space1D>;
using UMap22 = VNCS::UMap<VNCS::Space2D, VNCS::Space2D>;
}  // namespace Sim2D

namespace Sim3D
{
using UMap31 = VNCS::UMap<VNCS::Space3D, VNCS::Space1D>;
using UMap32 = VNCS::UMap<VNCS::Space3D, VNCS::Space2D>;
using UMap33 = VNCS::UMap<VNCS::Space3D, VNCS::Space3D>;
}  // namespace Sim3D

// namespace Sim3D

// namespace Sim3D
}  // namespace VNCS

#endif  // VNCS_DEFORMATIONGRADIENTMAP_H
