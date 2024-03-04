#ifndef VNCS_PRECONDITIONER_H
#define VNCS_PRECONDITIONER_H

#include <SofaBaseLinearSolver/config.h>

#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/simulation/VectorOperations.h>

#include <spdlog/spdlog.h>

namespace VNCS
{
class Preconditioner : public sofa::core::behavior::LinearSolver
{
    using Inherit = sofa::core::behavior::LinearSolver;

public:
    SOFA_CLASS(Preconditioner, sofa::core::behavior::LinearSolver);

    Preconditioner() = default;
    void init() override {}

    void setSystemRHVector(sofa::core::MultiVecDerivId b) { m_b = b; }
    void setSystemLHVector(sofa::core::MultiVecDerivId x) { m_x = x; }

    void solveSystem()
    {
        spdlog::get("VNCS")->info("Preconditioner");
        msg_info("b") << m_b;
        msg_info("x") << m_x;
    }

    void resetSystem() {}
    void setSystemMBKMatrix(const sofa::core::MechanicalParams *mparams) { m_mparams = *mparams; }

private:
    sofa::core::MechanicalParams m_mparams;
    sofa::core::MultiVecDerivId m_b;
    sofa::core::MultiVecDerivId m_x;
};
}  // namespace VNCS

#endif  // VNCS_PRECONDITIONER_H
