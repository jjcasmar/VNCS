#ifndef VNCS_SAMPLINGPOINTS_H
#define VNCS_SAMPLINGPOINTS_H

#include <memory>
#include <optional>
#include <filesystem>
#include <fstream>

#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>

using json = nlohmann::json;

namespace
{
const auto KEY_POINT = std::string("point");
const auto KEY_WEIGHT = std::string("weight");
const auto KEY_ALPHA = std::string("alpha");
const auto KEY_SHAPE_FUNCTIONS = std::string("shapeFunctions");
const auto KEY_COARSE = std::string("coarse");
const auto KEY_FINE = std::string("fine");
const auto KEY_NODE_INDEX = std::string("nodeIndex");
const auto KEY_SHAPE_FUNCTION_VALUE = std::string("v");
const auto KEY_SHAPE_FUNCTION_D_VALUE = std::string("dv");
}  // namespace

// Add a couple of customization points to support std::optional in nlohmann::json
namespace nlohmann
{
template <typename T>
void from_json(const nlohmann::json &j, std::optional<T> &t)
{
    if (!j.is_null())
        t = j.get<T>();
}

template <typename T>
void to_json(nlohmann::json &j, const std::optional<T> &t)
{
    if (t.has_value())
        j = t.value();
}
}  // namespace nlohmann

namespace VNCS
{
template <typename MaterialSpace>
struct ShapeFunction {
    using Real = typename MaterialSpace::VecType::Real;
    int nodeIndex;
    Real v;
    std::optional<std::array<Real, MaterialSpace::dim>> dv;
};

template <typename MaterialSpace>
struct Sampler {
    using Real = typename MaterialSpace::VecType::Real;
    using Coord = typename MaterialSpace::VecType::Coord;

    Coord x;
    Real w;
    Real a;
    std::optional<int> coarseFaceIdx;
    std::optional<int> fineFaceIdx;
    std::vector<ShapeFunction<MaterialSpace>> coarseShapeFunctions;
    std::vector<ShapeFunction<MaterialSpace>> fineShapeFunctions;
};

template <typename Space>
using SamplingPoints = std::vector<Sampler<Space>>;

template <typename T>
void from_json(const json &j, ShapeFunction<T> &shapeFunction)
{
    j.at("nodeIndex").get_to(shapeFunction.nodeIndex);
    j.at("v").get_to(shapeFunction.v);
    if (j.contains("dv"))
        j.at("dv").get_to(shapeFunction.dv);
}

template <typename T>
void to_json(json &j, const ShapeFunction<T> &shapeFunction)
{
    j["nodeIndex"] = shapeFunction.nodeIndex;
    j["v"] = shapeFunction.v;
    if (shapeFunction.dv)
        j["dv"] = shapeFunction.dv;
}

template <typename T>
void from_json(const json &j, Sampler<T> &sampler)
{
    j.at("x").get_to(sampler.x);
    if (j.contains("w"))
        j.at("w").get_to(sampler.w);
    j.at("a").get_to(sampler.a);
    if (j.contains("coarseFaceIndex"))
        j.at("coarseFaceIndex").get_to(sampler.coarseFaceIdx);
    if (j.contains("fineFaceIndex"))
        j.at("fineFaceIndex").get_to(sampler.fineFaceIdx);
    const auto shapeFunctions = j.at("shapeFunctions");
    j.at("shapeFunctions").at("coarse").get_to(sampler.coarseShapeFunctions);
    j.at("shapeFunctions").at("fine").get_to(sampler.fineShapeFunctions);
}

template <typename T>
void to_json(json &j, const Sampler<T> &sampler)
{
    j["x"] = sampler.x;
    j["w"] = sampler.w;
    j["a"] = sampler.a;
    if (sampler.coarseFaceIdx)
        j["coarseFaceIndex"] = sampler.coarseFaceIdx;
    if (sampler.fineFaceIdx)
        j["fineFaceIndex"] = sampler.fineFaceIdx;
    j["shapeFunctions"]["coarse"] = sampler.coarseShapeFunctions;
    j["shapeFunctions"]["fine"] = sampler.fineShapeFunctions;
}

template <typename T>
SamplingPoints<T> loadSamplingPoints(const std::filesystem::path &path)
{
    if (std::filesystem::exists(path)) {
        std::ifstream jsonData(path);

        json j;
        jsonData >> j;

        SamplingPoints<T> samplingPoints = j.get<SamplingPoints<T>>();
        return samplingPoints;
    }
    return {};
}

template <typename T>
void saveSamplingPoints(const std::filesystem::path &path, const SamplingPoints<T> &samplingPoints)
{
    std::ofstream jsonData(path);

    json j = samplingPoints;
    jsonData << j;
}

}  // namespace VNCS

#endif  // SAMPLINGPOINTS_H
