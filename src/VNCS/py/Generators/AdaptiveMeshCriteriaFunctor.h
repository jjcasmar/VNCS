#ifndef VNCS_PY_ADAPTIVEMESHCRITERIAFUNCTOR_H
#define VNCS_PY_ADAPTIVEMESHCRITERIAFUNCTOR_H

#include <VNCS/Generators/AdaptiveMeshCriteria.h>
#include <pybind11/pybind11.h>

#include <VNCS/Spaces.h>

namespace VNCS
{
namespace Generators
{
namespace py
{
class AdaptiveMeshCriteriaFunctor : public VNCS::Generators::AdaptiveMeshCriteriaFunctor
{
    using K = VNCS::Space2D::K;
    using Real = VNCS::Space2D::Real;

public:
    std::pair<Real, Real> operator()(const VNCS::Space2D::Point &p) const;
};

void adaptiveMeshCriteriaFunctor(pybind11::module &m);
}  // namespace py
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_PY_ADAPTIVEMESHCRITERIAFUNCTOR_H
