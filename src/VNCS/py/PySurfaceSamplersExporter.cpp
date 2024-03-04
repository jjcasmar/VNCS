#include "PySurfaceSamplersExporter.h"

#include <VNCS/SurfaceSamplersExporter.h>
#include <sofa/core/sptr.h>

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);

void VNCS::py::module::surfaceSamplersExporter(pybind11::module &m)
{
    pybind11::class_<VNCS::SurfaceSamplersExporter, sofa::core::sptr<VNCS::SurfaceSamplersExporter>>(
        m, "SurfaceSamplersExporter")
        .def(pybind11::init(
            [](pybind11::args &args, pybind11::kwargs &kwargs) { return new VNCS::SurfaceSamplersExporter(); }))
        .def_property(
            "path",
            [](VNCS::SurfaceSamplersExporter *exporter) { return ""; },
            [](VNCS::SurfaceSamplersExporter *exporter, std::string path) { exporter->setPath(path); })
        .def_property("blendingField",
                      &VNCS::SurfaceSamplersExporter::blendingField,
                      &VNCS::SurfaceSamplersExporter::setBlendingField);
}
