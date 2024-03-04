#include "BendingEnergy.h"

#include <VNCS/Sim2D/BendingForceField.h>
#include <sofa/core/sptr.h>

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim2D::py::bendingEnergy(pybind11::module &m)
{
    pybind11::class_<VNCS::Sim2D::BendingForceField,
                     sofa::core::objectmodel::BaseObject,
                     sofa::core::sptr<VNCS::Sim2D::BendingForceField>>(m, "Bending")
        .def(pybind11::init<>())
        .def_property("edgeMeshPath", nullptr, [](VNCS::Sim2D::BendingForceField &f, const std::string &s) {f.setEdgeMeshPath(s);})
        .def_property("stiffness",
                      &VNCS::Sim2D::BendingForceField::bendingStiffness,
                      &VNCS::Sim2D::BendingForceField::setBendingStiffness);
}
