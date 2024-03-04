#ifndef VNCS_FMAP_H
#define VNCS_FMAP_H

#include <sofa/core/MultiMapping.h>
#include <sofa/core/MultiMapping.inl>

#include <VNCS/Spaces.h>
#include <VNCS/DataExtensions.h>
#include <VNCS/DeformationGradientTypes.h>

#include <VNCS/SamplingPoints.h>
#include <sofa/core/VecId.h>

#include <spdlog/fmt/ostr.h>

namespace VNCS
{
template <typename F>
class FMap : public sofa::core::MultiMapping<typename F::WorldSpace::VecType, F>
{
    using Inherit = sofa::core::MultiMapping<typename F::WorldSpace::VecType, F>;
    using Real = typename F::WorldSpace::VecType::Real;
    using Vector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

public:
    SOFA_CLASS(SOFA_TEMPLATE(VNCS::FMap, F),  //
               SOFA_TEMPLATE2(sofa::core::MultiMapping, typename F::WorldSpace::VecType, F));

    FMap()
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

                // A deformation gradient is a WorldSpace::dim x MaterialSpace::dim matrix,
                //
                //     W0/M0   W1/M0 ...
                // F = W0/M1   W1/M1 ...
                //     ...     ...   ...
                //
                // I get a vector view of it which will be
                //
                //      W0/M0
                //  F = W1/M0
                //      ...
                //      W0/M1
                //      W1/M1
                //      ...

                for (auto dM = 0; dM < F::MaterialSpace::dim; ++dM) {
                    for (auto dW = 0; dW < F::WorldSpace::dim; ++dW) {
                        phiCTriplets.emplace_back(
                            F::WorldSpace::dim * F::MaterialSpace::dim * i + dM * F::WorldSpace::dim + dW,
                            F::WorldSpace::dim * shapeFunction.nodeIndex + dW,
                            a * shapeFunction.dv.value()[dM]);
                    }
                }
                //                const auto dX = 0;
                //                const auto dY = 1;

                //                const auto ux = 2 * shapeFunction.nodeIndex + 0;
                //                const auto uy = 2 * shapeFunction.nodeIndex + 1;

                //                // dux / dX
                //                phiCTriplets.emplace_back(4 * i + 0, ux, a * shapeFunction.dv[dX]);
                //                phiCTriplets.emplace_back(4 * i + 0, ux, da[dX] * shapeFunction.v);

                //                // duy / dX
                //                phiCTriplets.emplace_back(4 * i + 1, uy, a * shapeFunction.dv[dX]);
                //                phiCTriplets.emplace_back(4 * i + 1, uy, da[dX] * shapeFunction.v);

                //                // dux / dY
                //                phiCTriplets.emplace_back(4 * i + 2, ux, a * shapeFunction.dv[dY]);
                //                phiCTriplets.emplace_back(4 * i + 2, ux, da[dY] * shapeFunction.v);

                //                // duy / dY
                //                phiCTriplets.emplace_back(4 * i + 3, uy, a * shapeFunction.dv[dY]);
                //                phiCTriplets.emplace_back(4 * i + 3, uy, da[dY] * shapeFunction.v);
            }

            const auto &fineShapeFunctions = sampler.fineShapeFunctions;
            for (int j = 0; j < fineShapeFunctions.size(); ++j) {
                const auto &shapeFunction = fineShapeFunctions[j];

                for (auto dM = 0; dM < F::MaterialSpace::dim; ++dM) {
                    for (auto dW = 0; dW < F::WorldSpace::dim; ++dW) {
                        phiFTriplets.emplace_back(
                            F::WorldSpace::dim * F::MaterialSpace::dim * i + dM * F::WorldSpace::dim + dW,
                            F::WorldSpace::dim * shapeFunction.nodeIndex + dW,
                            (1 - a) * shapeFunction.dv.value()[dM]);
                    }
                }

                //                const auto dX = 0;
                //                const auto dY = 1;

                //                const auto ux = 2 * shapeFunction.nodeIndex + 0;
                //                const auto uy = 2 * shapeFunction.nodeIndex + 1;

                //                // dux / dX
                //                phiFTriplets.emplace_back(4 * i + 0, ux, (1 - a) * shapeFunction.dv[dX]);
                //                phiFTriplets.emplace_back(4 * i + 0, ux, -da[dX] * shapeFunction.v);

                //                // duy / dX
                //                phiFTriplets.emplace_back(4 * i + 1, uy, (1 - a) * shapeFunction.dv[dX]);
                //                phiFTriplets.emplace_back(4 * i + 1, uy, -da[dX] * shapeFunction.v);

                //                // dux / dY
                //                phiFTriplets.emplace_back(4 * i + 2, ux, (1 - a) * shapeFunction.dv[dY]);
                //                phiFTriplets.emplace_back(4 * i + 2, ux, -da[dY] * shapeFunction.v);

                //                // duy / dY
                //                phiFTriplets.emplace_back(4 * i + 3, uy, (1 - a) * shapeFunction.dv[dY]);
                //                phiFTriplets.emplace_back(4 * i + 3, uy, -da[dY] * shapeFunction.v);
            }
        }

        if (m_coarseInput) {
            m_phiC.resize(F::MaterialSpace::dim * F::WorldSpace::dim * samplers.size(),
                          F::WorldSpace::dim * this->fromModels[m_coarseIndex.value()]->getSize());

            m_phiC.setFromTriplets(std::begin(phiCTriplets), std::end(phiCTriplets));
        }

        if (m_fineInput) {
            m_phiF.resize(F::MaterialSpace::dim * F::WorldSpace::dim * samplers.size(),
                          F::WorldSpace::dim * this->fromModels[m_fineIndex.value()]->getSize());
            m_phiF.setFromTriplets(std::begin(phiFTriplets), std::end(phiFTriplets));
        }

        Inherit::init();
    }

    void apply(const sofa::core::MechanicalParams *mparams,
               const sofa::helper::vector<typename Inherit::OutDataVecCoord *> &dataVecOutPos,
               const sofa::helper::vector<const typename Inherit::InDataVecCoord *> &dataVecInPos) final
    {
        auto &deformationGradientAtSampler = make_write_accessor(*dataVecOutPos[0]);

        Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>> u(
            (deformationGradientAtSampler[0].data()),  //
            F::coord_total_size * F::spatial_dimensions * deformationGradientAtSampler.size());
        u.setZero();

        if (m_coarseInput) {
            const auto &coarseDeformationAtNodes = make_read_accessor(*dataVecInPos[m_coarseIndex.value()]);
            if (coarseDeformationAtNodes.size()) {
                const auto &x0 =
                    *this->fromModels[m_coarseIndex.value()]->read(sofa::core::ConstVecCoordId::restPosition());
                const auto &restPositions = make_read_accessor(x0);

                Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> uc(
                    &(coarseDeformationAtNodes[0][0]),  //
                    F::WorldSpace::dim * coarseDeformationAtNodes.size());

                Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> x0Eig(
                    &(restPositions[0][0]),  //
                    F::WorldSpace::dim * coarseDeformationAtNodes.size());
                u += m_phiC * (x0Eig + uc);
            }
        }

        if (m_fineInput) {
            const auto &fineDeformationAtNodes = make_read_accessor(*dataVecInPos[m_fineIndex.value()]);
            if (fineDeformationAtNodes.size()) {
                const auto &x0 =
                    *this->fromModels[m_fineIndex.value()]->read(sofa::core::ConstVecCoordId::restPosition());
                const auto &restPositions = make_read_accessor(x0);

                Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> uf(
                    &(fineDeformationAtNodes[0][0]),  //
                    F::WorldSpace::dim * fineDeformationAtNodes.size());

                Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> x0Eig(
                    &(restPositions[0][0]),  //
                    F::WorldSpace::dim * fineDeformationAtNodes.size());
                u += m_phiF * (x0Eig + uf);
            }
        }
    }

    void applyJ(const sofa::core::MechanicalParams *mparams,
                const sofa::helper::vector<typename Inherit::OutDataVecDeriv *> &dataVecOutVel,
                const sofa::helper::vector<const typename Inherit::InDataVecDeriv *> &dataVecInVel) final
    {
        auto &deformationGradientAtSampler = make_write_accessor(*dataVecOutVel[0]);

        Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>> u(
            (deformationGradientAtSampler[0].data()),  //
            F::WorldSpace::dim * F::MaterialSpace::dim * deformationGradientAtSampler.size());
        u.setZero();

        if (m_coarseInput) {
            const auto &coarseDeformationAtNodes = make_read_accessor(*dataVecInVel[m_coarseIndex.value()]);
            if (coarseDeformationAtNodes.size()) {
                Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> uc(
                    &(coarseDeformationAtNodes[0][0]),  //
                    F::WorldSpace::dim * coarseDeformationAtNodes.size());
                u += m_phiC * uc;
            }
        }

        if (m_fineInput) {
            const auto &fineDeformationAtNodes = make_read_accessor(*dataVecInVel[m_fineIndex.value()]);
            if (fineDeformationAtNodes.size()) {
                Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> uf(
                    &(fineDeformationAtNodes[0][0]),  //
                    F::WorldSpace::dim * fineDeformationAtNodes.size());
                u += m_phiF * uf;
            }
        }
    }

    void applyJT(const sofa::core::MechanicalParams *mparams,
                 const sofa::helper::vector<typename Inherit::InDataVecDeriv *> &dataVecOutForce,
                 const sofa::helper::vector<const typename Inherit::OutDataVecDeriv *> &dataVecInForce) final
    {
        const auto &deformationGradientAtSampler = make_read_accessor(*dataVecInForce[0]);

        Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> u(
            (deformationGradientAtSampler[0].data()),  //
            F::WorldSpace::dim * F::MaterialSpace::dim * deformationGradientAtSampler.size());

        if (m_coarseInput) {
            auto &coarseDeformationAtNodes = make_write_accessor(*dataVecOutForce[m_coarseIndex.value()]);
            if (coarseDeformationAtNodes.size()) {
                Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>> uc(
                    &(coarseDeformationAtNodes[0][0]),  //
                    F::WorldSpace::dim * coarseDeformationAtNodes.size());
                uc += m_phiC.transpose() * u;
            }
        }

        if (m_fineInput) {
            auto &fineDeformationAtNodes = make_write_accessor(*dataVecOutForce[m_fineIndex.value()]);
            if (fineDeformationAtNodes.size()) {
                Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>> uf(
                    &(fineDeformationAtNodes[0][0]),  //
                    F::WorldSpace::dim * fineDeformationAtNodes.size());
                uf += m_phiF.transpose() * u;
            }
        }
    }

    void applyDJT(const sofa::core::MechanicalParams * /*mparams*/,
                  sofa::core::MultiVecDerivId /*inForce*/,
                  sofa::core::ConstMultiVecDerivId /*outForce*/) final
    {
    }

    void applyJT(const sofa::core::ConstraintParams * /* cparams */,
                 const sofa::helper::vector<typename Inherit::InDataMatrixDeriv *> & /* dataMatOutConst */,
                 const sofa::helper::vector<const typename Inherit::OutDataMatrixDeriv *> & /* dataMatInConst */) final
    {
    }

    const Eigen::SparseMatrix<Real> &phiC() const { return m_phiC; }
    const Eigen::SparseMatrix<Real> &phiF() const { return m_phiF; }

    std::shared_ptr<SamplingPoints<typename F::MaterialSpace>> samplingPoints() const { return m_samplingPoints; }
    void setSamplingPoints(const std::shared_ptr<SamplingPoints<typename F::MaterialSpace>> &samplingPoints)
    {
        m_samplingPoints = samplingPoints;
    }

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

    std::shared_ptr<SamplingPoints<typename F::MaterialSpace>> m_samplingPoints;
};
}  // namespace VNCS

#endif  // VNCS_DEFORMATIONGRADIENTMAP_H
