#ifndef VNCS_PY_PYMECHANICALOBJECT_H
#define VNCS_PY_PYMECHANICALOBJECT_H

#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <VNCS/DeformationGradientTypes.h>
#include <VNCS/Types.h>

namespace VNCS
{
namespace py
{
std::vector<VNCS::Sim3D::Coord> restPositions(pybind11::handle h);

std::vector<VNCS::Sim3D::Coord> positions(pybind11::handle h);

std::vector<VNCS::Sim3D::Deriv> velocities(pybind11::handle h);
std::vector<VNCS::Sim3D::Deriv> forces(pybind11::handle h);

std::vector<VNCS::F33D::Coord> deformationGradients(pybind11::handle h);

std::vector<VNCS::F33D::Deriv> stresses(pybind11::handle h);

namespace module
{
void mechanicalObject(pybind11::module &m);
}

}  // namespace py
}  // namespace VNCS

#endif  // PYMECHANICALOBJECT_H
