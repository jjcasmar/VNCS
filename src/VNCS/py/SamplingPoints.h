#ifndef VNCS_PY_SAMPLINGPOINTS_H
#define VNCS_PY_SAMPLINGPOINTS_H

#include <filesystem>
#include <memory>
#include <VNCS/SamplingPoints.h>
#include <pybind11/pybind11.h>

namespace VNCS
{
namespace py
{
template <typename MaterialSpace>
struct SamplingPointsHolder {
    SamplingPointsHolder(const std::filesystem::path &path);

    std::shared_ptr<VNCS::SamplingPoints<MaterialSpace>> samplingPoints;
};

template <typename MaterialSpace>
SamplingPointsHolder<MaterialSpace>::SamplingPointsHolder(const std::filesystem::path &path)
{
    samplingPoints =
        std::make_shared<VNCS::SamplingPoints<MaterialSpace>>(VNCS::loadSamplingPoints<MaterialSpace>(path));
}

void samplingPoints(pybind11::module &m);
}  // namespace py
}  // namespace VNCS

#endif  // SAMPLINGPOINTS_H
