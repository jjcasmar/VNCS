#ifndef VNCS_SIM2D_MECHANICALOBJECT_H
#define VNCS_SIM2D_MECHANICALOBJECT_H

#include <VNCS/Spaces.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/State.h>
#include <sofa/core/State.inl>
#include <VNCS/DataExtensions.h>

namespace sofa
{
namespace component
{
namespace container
{
template <>
inline void MechanicalObject<VNCS::Space2D::F22>::vThreshold(sofa::core::VecId v, double t)
{
}

template <>
inline double MechanicalObject<VNCS::Space2D::F22>::vDot(const core::ExecParams *,
                                                         core::ConstVecId a,
                                                         core::ConstVecId b)
{
    Real r = 0.0;
    return r;
}

template <>
inline double MechanicalObject<VNCS::Space2D::F22>::getConstraintJacobianTimesVecDeriv(unsigned int line,
                                                                                       core::ConstVecId id)
{
    SReal result = 0;
    return result;
}

template <>
inline void MechanicalObject<VNCS::Space2D::F22>::applyScale(const double sx, const double sy, const double sz)
{
}

template <>
inline bool MechanicalObject<VNCS::Space2D::F22>::isIndependent() const
{
    return false;
}

template <>
inline void MechanicalObject<VNCS::Space2D::F22>::storeResetState()
{
    // store a reset state only for independent dofs (mapped dofs are deduced from independent dofs)
    if (!isIndependent())
        return;
}

template <>
inline void MechanicalObject<VNCS::Space2D::F22>::buildIdentityBlocksInJacobian(
    const sofa::helper::vector<unsigned int> &list_n,
    core::MatrixDerivId &mID)
{
}

template <>
inline SReal MechanicalObject<VNCS::Space2D::F22>::vSum(const core::ExecParams *params, core::ConstVecId a, unsigned l)
{
    Real r = 0.0;
    return r;
}

template <>
inline SReal MechanicalObject<VNCS::Space2D::F22>::vMax(const core::ExecParams *params, core::ConstVecId a)
{
    Real r = 0.0;
    return r;
}

template <>
inline void MechanicalObject<VNCS::Space2D::F22>::getConstraintJacobian(const core::ConstraintParams *cParams,
                                                                        sofa::defaulttype::BaseMatrix *J,
                                                                        unsigned int &off)
{
}

template <>
inline void MechanicalObject<VNCS::Space2D::F21>::vThreshold(sofa::core::VecId v, double t)
{
}

template <>
inline double MechanicalObject<VNCS::Space2D::F21>::vDot(const core::ExecParams *,
                                                         core::ConstVecId a,
                                                         core::ConstVecId b)
{
    Real r = 0.0;
    return r;
}

template <>
inline double MechanicalObject<VNCS::Space2D::F21>::getConstraintJacobianTimesVecDeriv(unsigned int line,
                                                                                       core::ConstVecId id)
{
    SReal result = 0;
    return result;
}

template <>
inline void MechanicalObject<VNCS::Space2D::F21>::applyScale(const double sx, const double sy, const double sz)
{
}

template <>
inline bool MechanicalObject<VNCS::Space2D::F21>::isIndependent() const
{
    return false;
}

template <>
inline void MechanicalObject<VNCS::Space2D::F21>::storeResetState()
{
    // store a reset state only for independent dofs (mapped dofs are deduced from independent dofs)
    if (!isIndependent())
        return;
}

template <>
inline void MechanicalObject<VNCS::Space2D::F21>::buildIdentityBlocksInJacobian(
    const sofa::helper::vector<unsigned int> &list_n,
    core::MatrixDerivId &mID)
{
}

template <>
inline SReal MechanicalObject<VNCS::Space2D::F21>::vSum(const core::ExecParams *params, core::ConstVecId a, unsigned l)
{
    Real r = 0.0;
    return r;
}

template <>
inline SReal MechanicalObject<VNCS::Space2D::F21>::vMax(const core::ExecParams *params, core::ConstVecId a)
{
    Real r = 0.0;
    return r;
}

template <>
inline void MechanicalObject<VNCS::Space2D::F21>::getConstraintJacobian(const core::ConstraintParams *cParams,
                                                                        sofa::defaulttype::BaseMatrix *J,
                                                                        unsigned int &off)
{
}

template <>
inline void MechanicalObject<VNCS::Space3D::F33>::vThreshold(sofa::core::VecId v, double t)
{
}

template <>
inline double MechanicalObject<VNCS::Space3D::F33>::vDot(const core::ExecParams *,
                                                         core::ConstVecId a,
                                                         core::ConstVecId b)
{
    Real r = 0.0;
    return r;
}

template <>
inline double MechanicalObject<VNCS::Space3D::F33>::getConstraintJacobianTimesVecDeriv(unsigned int line,
                                                                                       core::ConstVecId id)
{
    SReal result = 0;
    return result;
}

template <>
inline void MechanicalObject<VNCS::Space3D::F33>::applyScale(const double sx, const double sy, const double sz)
{
}

template <>
inline bool MechanicalObject<VNCS::Space3D::F33>::isIndependent() const
{
    return false;
}

template <>
inline void MechanicalObject<VNCS::Space3D::F33>::storeResetState()
{
    // store a reset state only for independent dofs (mapped dofs are deduced from independent dofs)
    if (!isIndependent())
        return;
}

template <>
inline void MechanicalObject<VNCS::Space3D::F33>::buildIdentityBlocksInJacobian(
    const sofa::helper::vector<unsigned int> &list_n,
    core::MatrixDerivId &mID)
{
}

template <>
inline SReal MechanicalObject<VNCS::Space3D::F33>::vSum(const core::ExecParams *params, core::ConstVecId a, unsigned l)
{
    Real r = 0.0;
    return r;
}

template <>
inline SReal MechanicalObject<VNCS::Space3D::F33>::vMax(const core::ExecParams *params, core::ConstVecId a)
{
    Real r = 0.0;
    return r;
}

template <>
inline void MechanicalObject<VNCS::Space3D::F33>::getConstraintJacobian(const core::ConstraintParams *cParams,
                                                                        sofa::defaulttype::BaseMatrix *J,
                                                                        unsigned int &off)
{
}

template <>
inline void MechanicalObject<VNCS::Space3D::F32>::vThreshold(sofa::core::VecId v, double t)
{
}

template <>
inline double MechanicalObject<VNCS::Space3D::F32>::vDot(const core::ExecParams *,
                                                         core::ConstVecId a,
                                                         core::ConstVecId b)
{
    Real r = 0.0;
    return r;
}

template <>
inline double MechanicalObject<VNCS::Space3D::F32>::getConstraintJacobianTimesVecDeriv(unsigned int line,
                                                                                       core::ConstVecId id)
{
    SReal result = 0;
    return result;
}

template <>
inline void MechanicalObject<VNCS::Space3D::F32>::applyScale(const double sx, const double sy, const double sz)
{
}

template <>
inline bool MechanicalObject<VNCS::Space3D::F32>::isIndependent() const
{
    return false;
}

template <>
inline void MechanicalObject<VNCS::Space3D::F32>::storeResetState()
{
    // store a reset state only for independent dofs (mapped dofs are deduced from independent dofs)
    if (!isIndependent())
        return;
}

template <>
inline void MechanicalObject<VNCS::Space3D::F32>::buildIdentityBlocksInJacobian(
    const sofa::helper::vector<unsigned int> &list_n,
    core::MatrixDerivId &mID)
{
}

template <>
inline SReal MechanicalObject<VNCS::Space3D::F32>::vSum(const core::ExecParams *params, core::ConstVecId a, unsigned l)
{
    Real r = 0.0;
    return r;
}

template <>
inline SReal MechanicalObject<VNCS::Space3D::F32>::vMax(const core::ExecParams *params, core::ConstVecId a)
{
    Real r = 0.0;
    return r;
}

template <>
inline void MechanicalObject<VNCS::Space3D::F32>::getConstraintJacobian(const core::ConstraintParams *cParams,
                                                                        sofa::defaulttype::BaseMatrix *J,
                                                                        unsigned int &off)
{
}

template <>
inline void MechanicalObject<VNCS::Space3D::F31>::vThreshold(sofa::core::VecId v, double t)
{
}

template <>
inline double MechanicalObject<VNCS::Space3D::F31>::vDot(const core::ExecParams *,
                                                         core::ConstVecId a,
                                                         core::ConstVecId b)
{
    Real r = 0.0;
    return r;
}

template <>
inline double MechanicalObject<VNCS::Space3D::F31>::getConstraintJacobianTimesVecDeriv(unsigned int line,
                                                                                       core::ConstVecId id)
{
    SReal result = 0;
    return result;
}

template <>
inline void MechanicalObject<VNCS::Space3D::F31>::applyScale(const double sx, const double sy, const double sz)
{
}

template <>
inline bool MechanicalObject<VNCS::Space3D::F31>::isIndependent() const
{
    return false;
}

template <>
inline void MechanicalObject<VNCS::Space3D::F31>::storeResetState()
{
    // store a reset state only for independent dofs (mapped dofs are deduced from independent dofs)
    if (!isIndependent())
        return;
}

template <>
inline void MechanicalObject<VNCS::Space3D::F31>::buildIdentityBlocksInJacobian(
    const sofa::helper::vector<unsigned int> &list_n,
    core::MatrixDerivId &mID)
{
}

template <>
inline SReal MechanicalObject<VNCS::Space3D::F31>::vSum(const core::ExecParams *params, core::ConstVecId a, unsigned l)
{
    Real r = 0.0;
    return r;
}

template <>
inline SReal MechanicalObject<VNCS::Space3D::F31>::vMax(const core::ExecParams *params, core::ConstVecId a)
{
    Real r = 0.0;
    return r;
}

template <>
inline void MechanicalObject<VNCS::Space3D::F31>::getConstraintJacobian(const core::ConstraintParams *cParams,
                                                                        sofa::defaulttype::BaseMatrix *J,
                                                                        unsigned int &off)
{
}
}  // namespace container
}  // namespace component
}  // namespace sofa

#include <SofaBaseMechanics/MechanicalObject.inl>

namespace VNCS
{
namespace Sim2D
{
using MO = sofa::component::container::MechanicalObject<VNCS::Space2D::VecType>;
using FMO21 = sofa::component::container::MechanicalObject<VNCS::Space2D::F21>;
using FMO22 = sofa::component::container::MechanicalObject<VNCS::Space2D::F22>;

}  // namespace Sim2D
namespace Sim3D
{
using MO = sofa::component::container::MechanicalObject<VNCS::Space3D::VecType>;
using FMO33 = sofa::component::container::MechanicalObject<VNCS::Space3D::F33>;
using FMO32 = sofa::component::container::MechanicalObject<VNCS::Space3D::F32>;
using FMO31 = sofa::component::container::MechanicalObject<VNCS::Space3D::F31>;

}  // namespace Sim3D
}  // namespace VNCS

template class sofa::component::container::MechanicalObject<VNCS::Space2D::VecType>;
template class sofa::component::container::MechanicalObject<VNCS::Space3D::VecType>;
template class sofa::core::State<VNCS::Space2D::F22>;
template class sofa::component::container::MechanicalObject<VNCS::Space2D::F22>;
template class sofa::core::State<VNCS::Space3D::F33>;
template class sofa::component::container::MechanicalObject<VNCS::Space3D::F33>;
template class sofa::core::State<VNCS::Space3D::F32>;
template class sofa::component::container::MechanicalObject<VNCS::Space3D::F32>;
template class sofa::core::State<VNCS::Space3D::F31>;
template class sofa::component::container::MechanicalObject<VNCS::Space3D::F31>;
template class sofa::core::State<VNCS::Space2D::F21>;
template class sofa::component::container::MechanicalObject<VNCS::Space2D::F21>;
#endif  // VNCS_SIM2D_MECHANICALOBJECT_H
