#include "StretchEnergy.h"

#include <VNCS/Sim3D/StretchForceField.h>
#include <sofa/core/sptr.h>

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::Sim3D::py::stretchEnergy(pybind11::module &m)
{
    pybind11::class_<VNCS::Sim3D::StretchForceField,
                     sofa::core::objectmodel::BaseObject,
                     sofa::core::sptr<VNCS::Sim3D::StretchForceField>>(m, "Stretch")
        .def(pybind11::init<>())
        .def_property("edgeMeshPath",
                      nullptr,
                      [](VNCS::Sim3D::StretchForceField &f, const std::string &s) { f.setEdgeMeshPath(s); })
        .def_property("stiffness",
                      &VNCS::Sim3D::StretchForceField::stretchStiffness,
                      &VNCS::Sim3D::StretchForceField::setStretchStiffness);
}
