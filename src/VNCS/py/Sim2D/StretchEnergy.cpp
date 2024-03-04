#include "StretchEnergy.h"
#include <sofa/core/sptr.h>

#include <VNCS/Sim2D/StretchForceField.h>

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);
void VNCS::Sim2D::py::stretchEnergy(pybind11::module &m)
{
    pybind11::class_<VNCS::Sim2D::StretchForceField,
                     sofa::core::objectmodel::BaseObject,
                     sofa::core::sptr<VNCS::Sim2D::StretchForceField>>(m, "Stretch")
        .def(pybind11::init<>())
        .def_property("edgeMeshPath", nullptr, [](VNCS::Sim2D::StretchForceField &f, const std::string &s) {f.setEdgeMeshPath(s);})
        .def_property("stiffness",
                      &VNCS::Sim2D::StretchForceField::stretchStiffness,
                      &VNCS::Sim2D::StretchForceField::setStretchStiffness);
}
