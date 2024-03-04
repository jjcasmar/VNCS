#ifndef VNCS_CGSOLVER_H
#define VNCS_CGSOLVER_H

#include <chrono>
#include <optional>

#include <VNCS/Types.h>

#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/core/behavior/MultiVec.h>
#include <sofa/core/MechanicalParams.h>
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>
#include <sofa/helper/OptionsGroup.h>

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>

namespace VNCS
{
class ConjugateGradientSolver : public sofa::core::behavior::LinearSolver
{
public:
    SOFA_CLASS(ConjugateGradientSolver, LinearSolver);
    ConjugateGradientSolver();

    void resetSystem() final;

    void setSystemMBKMatrix(const sofa::core::MechanicalParams *mparams) final;

    void setSystemRHVector(sofa::core::MultiVecDerivId b_id) final;

    void setSystemLHVector(sofa::core::MultiVecDerivId x_id) final;

    void solveSystem() final;

    int iterations() const;
    void setIterations(int newIterations);

    VNCS::Real toleranceThreshold() const;
    void setToleranceThreshold(VNCS::Real newToleranceThreshold);

    const std::optional<Eigen::SparseMatrix<VNCS::Real>> &preconditioner() const;
    void setPreconditioner(std::optional<Eigen::SparseMatrix<VNCS::Real>> newPreconditioner);

    const std::vector<VNCS::Real> &squaredResiduals() const;

    bool isVerbose() const;
    void setIsVerbose(bool newIsVerbose);

    const std::vector<int> &durations() const;
    void clearDurations();

private:
    void solve(sofa::core::behavior::MultiVecDeriv &b, sofa::core::behavior::MultiVecDeriv &x);

    /// INPUTS
    int m_iterations;
    VNCS::Real m_toleranceThreshold;

    /// Private members
    ///< The mechanical parameters containing the m, b and k coefficients.
    sofa::core::MechanicalParams m_mparams;

    ///< The identifier of the b vector
    sofa::core::MultiVecDerivId m_bId;

    ///< The identifier of the x vector
    sofa::core::MultiVecDerivId m_xId;

    std::vector<VNCS::Real> m_squaredResiduals;

    std::optional<Eigen::SparseMatrix<VNCS::Real>> m_preconditioner;

    std::vector<int> m_durations;
};

}  // namespace VNCS

#endif  // VNCS_CGSOLVER_H
