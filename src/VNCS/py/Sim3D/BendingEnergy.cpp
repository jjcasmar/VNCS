#include "BendingEnergy.h"

#include <VNCS/Sim3D/BendingForceField.h>
#include <sofa/core/sptr.h>

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim3D::py::bendingEnergy(pybind11::module &m)
{
    pybind11::class_<VNCS::Sim3D::BendingForceField,
                     sofa::core::objectmodel::BaseObject,
                     sofa::core::sptr<VNCS::Sim3D::BendingForceField>>(m, "Bending")
        .def(pybind11::init<>())
        .def_property("edgeMeshPath",
                      nullptr,
                      [](VNCS::Sim3D::BendingForceField &f, const std::string &s) { f.setEdgeMeshPath(s); })
        .def_property("stiffness",
                      &VNCS::Sim3D::BendingForceField::bendingStiffness,
                      &VNCS::Sim3D::BendingForceField::setBendingStiffness);
}
