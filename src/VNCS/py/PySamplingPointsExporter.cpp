#include "PySamplingPointsExporter.h"

#include <VNCS/SamplingPointsExporter.h>
#include <sofa/core/sptr.h>

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);
void VNCS::py::module::samplingPointsExporter(pybind11::module &m)
{
    pybind11::class_<VNCS::SamplingPointsExporter, sofa::core::sptr<VNCS::SamplingPointsExporter>>(
        m, "SamplingPointsExporter")
        .def(pybind11::init(
            [](pybind11::args &args, pybind11::kwargs &kwargs) { return new VNCS::SamplingPointsExporter(); }))
        .def_property(
            "finePath",
            [](VNCS::SamplingPointsExporter *exporter) { return ""; },
            [](VNCS::SamplingPointsExporter *exporter, std::string path) { exporter->setFineFilePath(path); })
        .def_property(
            "coarsePath",
            [](VNCS::SamplingPointsExporter *exporter) { return ""; },
            [](VNCS::SamplingPointsExporter *exporter, std::string path) { exporter->setCoarseFilePath(path); })
        .def_property("blendingField",
                      &VNCS::SamplingPointsExporter::blendingField,
                      &VNCS::SamplingPointsExporter::setBlendingField);
}
