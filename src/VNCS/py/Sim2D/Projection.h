#ifndef VNCS_PY_PROJECTION_H
#define VNCS_PY_PROJECTION_H

#include <pybind11/pybind11.h>

namespace VNCS
{
namespace Sim2D
{
namespace py
{
void projection(pybind11::module &m);
}  // namespace py
}  // namespace Sim2D
}  // namespace VNCS

#endif  //  VNCS_PY_PROJECTION_H
