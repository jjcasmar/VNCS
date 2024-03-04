#include "SamplingPoints.h"
#include <VNCS/Spaces.h>
#include <VNCS/SamplingPoints.h>

#include <pybind11/stl.h>

// Customization point for the Coord type
namespace sofa
{
namespace defaulttype
{
void from_json(const json &j, VNCS::Space1D::VecType::Coord &coord)
{
    j[0].get_to(coord.x());
}

void from_json(const json &j, VNCS::Space2D::VecType::Coord &coord)
{
    j[0].get_to(coord.x());
    j[1].get_to(coord.y());
}

void from_json(const json &j, VNCS::Space3D::VecType::Coord &coord)
{
    j[0].get_to(coord.x());
    j[1].get_to(coord.y());
    j[2].get_to(coord.y());
}
}  // namespace defaulttype
}  // namespace sofa

void VNCS::py::samplingPoints(pybind11::module &m)
{
    {
        using Sampler = VNCS::Sampler<VNCS::Space3D>;
        using ShapeFunction = VNCS::ShapeFunction<VNCS::Space3D>;

        pybind11::class_<ShapeFunction>(m, "ShapeFunction3")
            .def_readonly("nodeIndex", &ShapeFunction::nodeIndex)
            .def_readonly("v", &ShapeFunction::v)
            .def_readonly("dv", &ShapeFunction::dv);

        pybind11::class_<Sampler>(m, "Sampler3")
            .def_readonly("x", &Sampler::x)
            .def_readonly("w", &Sampler::w)
            .def_readonly("a", &Sampler::a)
            .def_readonly("coarseFunctions", &Sampler::coarseShapeFunctions)
            .def_readonly("fineFunctions", &Sampler::fineShapeFunctions);
    }

    {
        using Sampler = VNCS::Sampler<VNCS::Space2D>;
        using ShapeFunction = VNCS::ShapeFunction<VNCS::Space2D>;

        pybind11::class_<ShapeFunction>(m, "ShapeFunction2")
            .def_readonly("nodeIndex", &ShapeFunction::nodeIndex)
            .def_readonly("v", &ShapeFunction::v)
            .def_readonly("dv", &ShapeFunction::dv);

        pybind11::class_<Sampler>(m, "Sampler2")
            .def_readonly("x", &Sampler::x)
            .def_readonly("w", &Sampler::w)
            .def_readonly("a", &Sampler::a)
            .def_readonly("coarseFunctions", &Sampler::coarseShapeFunctions)
            .def_readonly("fineFunctions", &Sampler::fineShapeFunctions);
    }

    {
        using Sampler = VNCS::Sampler<VNCS::Space1D>;
        using ShapeFunction = VNCS::ShapeFunction<VNCS::Space1D>;

        pybind11::class_<ShapeFunction>(m, "ShapeFunction1")
            .def_readonly("nodeIndex", &ShapeFunction::nodeIndex)
            .def_readonly("v", &ShapeFunction::v)
            .def_readonly("dv", &ShapeFunction::dv);

        pybind11::class_<Sampler>(m, "Sampler1")
            .def_readonly("x", &Sampler::x)
            .def_readonly("w", &Sampler::w)
            .def_readonly("a", &Sampler::a)
            .def_readonly("coarseFunctions", &Sampler::coarseShapeFunctions)
            .def_readonly("fineFunctions", &Sampler::fineShapeFunctions);
    }

    pybind11::class_<SamplingPointsHolder<VNCS::Space3D>>(m, "SamplingPoints3")
        .def(pybind11::init<std::string>())
        .def_property_readonly(
            "samplers", [](SamplingPointsHolder<VNCS::Space3D> &holder) -> const VNCS::SamplingPoints<VNCS::Space3D> & {
                return *holder.samplingPoints;
            });

    pybind11::class_<SamplingPointsHolder<VNCS::Space2D>>(m, "SamplingPoints2")
        .def(pybind11::init<std::string>())
        .def_property_readonly(
            "samplers", [](SamplingPointsHolder<VNCS::Space2D> &holder) -> const VNCS::SamplingPoints<VNCS::Space2D> & {
                return *holder.samplingPoints;
            });

    pybind11::class_<SamplingPointsHolder<VNCS::Space1D>>(m, "SamplingPoints1")
        .def(pybind11::init<std::string>())
        .def_property_readonly(
            "samplers", [](SamplingPointsHolder<VNCS::Space1D> &holder) -> const VNCS::SamplingPoints<VNCS::Space1D> & {
                return *holder.samplingPoints;
            });
}
