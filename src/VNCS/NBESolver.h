#ifndef VNCS_NBESOLVER_H
#define VNCS_NBESOLVER_H

#include <sofa/core/behavior/OdeSolver.h>
#include <VNCS/Types.h>

namespace VNCS
{
class NBESolver : public sofa::core::behavior::OdeSolver
{
public:
    SOFA_CLASS(NBESolver, sofa::core::behavior::OdeSolver);
    NBESolver();
    void init() override;

    void cleanup() override;

    void solve(const sofa::core::ExecParams *params,
               SReal dt,
               sofa::core::MultiVecCoordId xResult,
               sofa::core::MultiVecDerivId vResult) override;

    VNCS::Real rayleighStiffness() const;
    void setRayleighStiffness(VNCS::Real newRayleighStiffness);

    VNCS::Real rayleighMass() const;
    void setRayleighMass(VNCS::Real newRayleighMass);

    VNCS::Real newtonThreshold() const;
    void setNewtonThreshold(VNCS::Real newNewtonThreshold);

    int maxNewtonIterations() const;
    void setMaxNewtonIterations(int newMaxNewtonIterations);

    int lineSearchMaxIterations() const;
    void setLineSearchMaxIterations(int newLineSearchMaxIterations);

private:
    VNCS::Real m_rayleighStiffness;
    VNCS::Real m_rayleighMass;
    VNCS::Real m_newtonThreshold;
    int m_maxNewtonIterations;
    int m_lineSearchMaxIterations;
};
}  // namespace VNCS

#endif  // VNCS_NBESOLVER_H
