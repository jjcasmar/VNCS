#include "FileSamplingPoints.h"
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace
{
const auto KEY_POINT = std::string("point");
const auto KEY_WEIGHT = std::string("weight");
const auto KEY_ALPHA = std::string("alpha");
const auto KEY_D_ALPHA = std::string("da");
const auto KEY_SHAPE_FUNCTIONS = std::string("shapeFunctions");
const auto KEY_COARSE = std::string("coarse");
const auto KEY_FINE = std::string("fine");
const auto KEY_NODE_INDEX = std::string("nodeIndex");
const auto KEY_SHAPE_FUNCTION_VALUE = std::string("v");
const auto KEY_SHAPE_FUNCTION_D_VALUE = std::string("dv");
}  // namespace

VNCS::FileSamplingPoints::FileSamplingPoints()
    : m_source(initData(&m_source, "source", "Source file"))
{
}

void VNCS::FileSamplingPoints::init()
{
    std::ifstream in(m_source.getAbsolutePath());
    json samplers;
    in >> samplers;

    if (samplers.is_array()) {
        std::transform(
            std::begin(samplers), std::end(samplers), std::back_inserter(m_samplers), [](const auto &jSampler) {
                Sampler sampler;
                auto x = jSampler.value(KEY_POINT, std::array<double, 3>{});
                sampler.x[0] = x[0];
                sampler.x[1] = x[1];
                sampler.x[2] = x[2];
                sampler.w = jSampler.value(KEY_WEIGHT, 0.0);
                sampler.a = jSampler.value(KEY_ALPHA, 0.0);
                sampler.da = jSampler.value(KEY_D_ALPHA, std::array<double, 3>{});
                const auto jCoarseShapeFunctions =
                    jSampler.value(KEY_SHAPE_FUNCTIONS, json{}).value(KEY_COARSE, json{});
                const auto jFineShapeFunctions = jSampler.value(KEY_SHAPE_FUNCTIONS, json{}).value(KEY_FINE, json{});

                std::transform(std::begin(jCoarseShapeFunctions),
                               std::end(jCoarseShapeFunctions),
                               std::back_inserter(sampler.coarseShapeFunctions),
                               [](const auto &jShapeFunction) {
                                   Sampler::ShapeFunction shapeFunction;
                                   shapeFunction.nodeIndex = jShapeFunction.value(KEY_NODE_INDEX, 0.0);
                                   shapeFunction.v = jShapeFunction.value(KEY_SHAPE_FUNCTION_VALUE, 0.0);
                                   shapeFunction.dv =
                                       jShapeFunction.value(KEY_SHAPE_FUNCTION_D_VALUE, std::array<double, 3>{});
                                   return shapeFunction;
                               });

                std::transform(std::begin(jFineShapeFunctions),
                               std::end(jFineShapeFunctions),
                               std::back_inserter(sampler.fineShapeFunctions),
                               [](const auto &jShapeFunction) {
                                   Sampler::ShapeFunction shapeFunction;
                                   shapeFunction.nodeIndex = jShapeFunction.value(KEY_NODE_INDEX, 0.0);
                                   shapeFunction.v = jShapeFunction.value(KEY_SHAPE_FUNCTION_VALUE, 0.0);
                                   shapeFunction.dv =
                                       jShapeFunction.value(KEY_SHAPE_FUNCTION_D_VALUE, std::array<double, 3>{});
                                   return shapeFunction;
                               });
                return sampler;
            });
    }
    Inherit::init();
}
