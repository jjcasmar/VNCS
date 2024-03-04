#ifndef VNCS_PY_ENHACERELATIONMAP_H
#define VNCS_PY_ENHACERELATIONMAP_H

#include <VNCS/Spaces.h>
#include <vector>
#include <array>
#include <pybind11/pybind11.h>

namespace VNCS
{
namespace Sim2D
{
namespace py
{
void enhaceMap(pybind11::module &m);

}  // namespace py
}  // namespace Sim2D
}  // namespace VNCS

#endif  // ENHACERELATIONMAP_H
