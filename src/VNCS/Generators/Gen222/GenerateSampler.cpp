#include "GenerateSampler.h"

#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>
#include <range/v3/to_container.hpp>

#include <VNCS/Generators/IntegrationScheme.h>

VNCS::Generators::Gen222::GenerateSampler::GenerateSampler(
    const VNCS::Generators::PMPTriangleIntersection &coarseIntersectionFunc,
    const VNCS::Generators::PMPTriangleIntersection &fineIntersectionFunc)
    : m_coarseIntersectionFunc(coarseIntersectionFunc)
    , m_fineIntersectionFunc(fineIntersectionFunc)
{
}

VNCS::Sampler<VNCS::Space2D> VNCS::Generators::Gen222::GenerateSampler::operator()(
    const SamplingPointsExporter::SamplerDefinition &samplerDefinition,
    const SamplingPointsExporter::Mesh &coarseMesh,
    CoarseSampler)
{
    using SamplerDefinition = SamplingPointsExporter::SamplerDefinition;

    // Get the indices and points of the face
    const auto vIndices =
        coarseMesh.vertices_around_face(coarseMesh.halfedge(samplerDefinition.baseFace)) | ranges::to_vector;
    const auto v = vIndices | ranges::views::transform([&coarseMesh](const auto &vd) { return coarseMesh.point(vd); }) |
                   ranges::to_vector;

    // Create the J matrix. This matrix allows to convert the coordinates of the sampler point to an actual position
    // and to compute derivatives
    Eigen::Matrix3d J;
    J << v[1][0] - v[0][0], v[2][0] - v[0][0], v[0][0],  //
        v[1][1] - v[0][1], v[2][1] - v[0][1], v[0][1],   //
        0, 0, 1;
    const Eigen::Matrix3d Jinv = J.inverse();

    const auto baseInterpolationValues = VNCS::Sim2D::Tri_3::interpolation(samplerDefinition.coordinates);
    const auto baseInterpolationPartialValues = VNCS::Sim2D::Tri_3::interpolationPartial(samplerDefinition.coordinates);

    const Eigen::Vector2d basePoint = (J * (samplerDefinition.coordinates.homogeneous())).head<2>();

    VNCS::Sampler<VNCS::Space2D> sampler;
    sampler.x = {basePoint[0], basePoint[1]};
    sampler.w = samplerDefinition.weight;
    sampler.a = 1.0;
    sampler.coarseFaceIdx = samplerDefinition.baseFace.idx();

    auto &baseShapeFunctions = sampler.coarseShapeFunctions;

    for (const auto &[vIndex, phi, dPhi] : ranges::views::zip(vIndices,  //
                                                              baseInterpolationValues,
                                                              baseInterpolationPartialValues)) {
        VNCS::ShapeFunction<VNCS::Space2D> shapeFunction;
        shapeFunction.nodeIndex = vIndex;
        shapeFunction.v = phi;
        const Eigen::Vector2d dv = Jinv.transpose().block<2, 2>(0, 0) * dPhi;
        shapeFunction.dv = {dv[0], dv[1]};
        baseShapeFunctions.push_back(shapeFunction);
    }
    return sampler;
}

VNCS::Sampler<VNCS::Space2D> VNCS::Generators::Gen222::GenerateSampler::operator()(
    const SamplingPointsExporter::SamplerDefinition &samplerDefinition,
    const SamplingPointsExporter::Mesh &mesh,
    FineSampler)
{
    using SamplerDefinition = SamplingPointsExporter::SamplerDefinition;

    // Get the indices and points of the face
    const auto vIndices = mesh.vertices_around_face(mesh.halfedge(samplerDefinition.baseFace)) | ranges::to_vector;
    const auto v =
        vIndices | ranges::views::transform([&mesh](const auto &vd) { return mesh.point(vd); }) | ranges::to_vector;

    // Create the J matrix. This matrix allows to convert the coordinates of the sampler point to an actual position
    // and to compute derivatives
    Eigen::Matrix3d J;
    J << v[1][0] - v[0][0], v[2][0] - v[0][0], v[0][0],  //
        v[1][1] - v[0][1], v[2][1] - v[0][1], v[0][1],   //
        0, 0, 1;
    const Eigen::Matrix3d Jinv = J.inverse();

    const auto baseInterpolationValues = VNCS::Sim2D::Tri_3::interpolation(samplerDefinition.coordinates);
    const auto baseInterpolationPartialValues = VNCS::Sim2D::Tri_3::interpolationPartial(samplerDefinition.coordinates);

    const Eigen::Vector2d basePoint = (J * (samplerDefinition.coordinates.homogeneous())).head<2>();

    VNCS::Sampler<VNCS::Space2D> sampler;
    sampler.x = {basePoint[0], basePoint[1]};
    sampler.w = samplerDefinition.weight;
    sampler.a = 0.0;
    sampler.fineFaceIdx = samplerDefinition.baseFace.idx();

    auto &baseShapeFunctions = sampler.fineShapeFunctions;

    for (const auto &[vIndex, phi, dPhi] : ranges::views::zip(vIndices,  //
                                                              baseInterpolationValues,
                                                              baseInterpolationPartialValues)) {
        VNCS::ShapeFunction<VNCS::Space2D> shapeFunction;
        shapeFunction.nodeIndex = vIndex;
        shapeFunction.v = phi;
        const Eigen::Vector2d dv = Jinv.transpose().block<2, 2>(0, 0) * dPhi;
        shapeFunction.dv = {dv[0], dv[1]};
        baseShapeFunctions.push_back(shapeFunction);
    }
    return sampler;
}
