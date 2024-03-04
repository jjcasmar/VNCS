#include "PyMechanicalObject.h"

#include <sofa/core/objectmodel/BaseObject.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <SofaBaseMechanics/MechanicalObject.inl>
#include <sofa/core/State.h>
#include <sofa/core/State.inl>

#include <VNCS/DataExtensions.h>

template class sofa::component::container::MechanicalObject<VNCS::Sim3D>;
template class sofa::component::container::MechanicalObject<VNCS::F33D>;

std::vector<VNCS::Sim3D::Coord> VNCS::py::restPositions(pybind11::handle h)
{
    auto bo = h.cast<sofa::core::objectmodel::BaseObject *>();
    auto mo = dynamic_cast<sofa::component::container::MechanicalObject<VNCS::Sim3D> *>(bo);

    auto p = mo->readRestPositions();
    std::vector<VNCS::Sim3D::Coord> v;
    std::copy(std::begin(p), std::end(p), std::back_inserter(v));
    return v;
}

std::vector<VNCS::Sim3D::Coord> VNCS::py::positions(pybind11::handle h)
{
    auto bo = h.cast<sofa::core::objectmodel::BaseObject *>();
    auto mo = dynamic_cast<sofa::component::container::MechanicalObject<VNCS::Sim3D> *>(bo);

    auto p = mo->readPositions();
    std::vector<VNCS::Sim3D::Coord> v;
    std::copy(std::begin(p), std::end(p), std::back_inserter(v));
    return v;
}

std::vector<VNCS::Sim3D::Deriv> VNCS::py::velocities(pybind11::handle h)
{
    auto bo = h.cast<sofa::core::objectmodel::BaseObject *>();
    auto mo = dynamic_cast<sofa::component::container::MechanicalObject<VNCS::Sim3D> *>(bo);

    auto p = mo->readVelocities();
    std::vector<VNCS::Sim3D::Deriv> v;
    std::copy(std::begin(p), std::end(p), std::back_inserter(v));
    return v;
}

std::vector<VNCS::Sim3D::Deriv> VNCS::py::forces(pybind11::handle h)
{
    auto bo = h.cast<sofa::core::objectmodel::BaseObject *>();
    auto mo = dynamic_cast<sofa::component::container::MechanicalObject<VNCS::Sim3D> *>(bo);

    auto p = mo->readForces();
    std::vector<VNCS::Sim3D::Deriv> v;
    std::copy(std::begin(p), std::end(p), std::back_inserter(v));
    return v;
}

std::vector<VNCS::F33D::Coord> VNCS::py::deformationGradients(pybind11::handle h)
{
    auto bo = h.cast<sofa::core::objectmodel::BaseObject *>();
    auto mo = dynamic_cast<sofa::component::container::MechanicalObject<VNCS::F33D> *>(bo);

    auto p = mo->readPositions();
    std::vector<VNCS::F33D::Coord> v;
    std::copy(std::begin(p), std::end(p), std::back_inserter(v));
    return v;
}

std::vector<VNCS::F33D::Deriv> VNCS::py::stresses(pybind11::handle h)
{
    auto bo = h.cast<sofa::core::objectmodel::BaseObject *>();
    auto mo = dynamic_cast<sofa::component::container::MechanicalObject<VNCS::F33D> *>(bo);

    auto p = mo->readForces();
    std::vector<VNCS::F33D::Deriv> v;
    std::copy(std::begin(p), std::end(p), std::back_inserter(v));
    return v;
}

void VNCS::py::module::mechanicalObject(pybind11::module &m)
{
    m.def("restPositions", &restPositions);
    m.def("positions", &positions);
    m.def("velocities", &velocities);
    m.def("forces", &forces);
    m.def("deformationGradients", &deformationGradients);
    m.def("stresses", &stresses);
}
