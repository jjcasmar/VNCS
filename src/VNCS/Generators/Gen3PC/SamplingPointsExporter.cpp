#include "SamplingPointsExporter.h"

#include <spdlog/spdlog.h>
#include <spdlog/pattern_formatter.h>

#include <VNCS/Generators/TetraMeshIntersection.h>
#include <VNCS/Generators/KMeans.h>
#include <VNCS/Generators/IntegrationScheme.h>

#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/action/sort.hpp>
#include <range/v3/view/group_by.hpp>

namespace
{
struct SamplerDefinition {
    Eigen::Vector3d coordinates;
    VNCS::Real weight;
    int tetraIdx;
};

    VNCS::Sampler<VNCS::Space3D> generateSampler(
    const SamplerDefinition &samplerDefinition,
    const VNCS::Space3D::TetraMesh &mesh)
{
    // Get the indices and points of the face
    const auto vIndices = mesh.tetras[samplerDefinition.tetraIdx];
    const auto v = vIndices |
                   ranges::views::transform([&mesh](const auto &vd) { return mesh.points[vd]; }) |
                   ranges::to_vector;

    // Create the J matrix. This matrix allows to convert the coordinates of the sampler point to an actual position
    // and to compute derivatives
    Eigen::Matrix4d J;
    J << v[1][0] - v[0][0], v[2][0] - v[0][0], v[3][0] - v[0][0],
        v[0][0],  //
        v[1][1] - v[0][1], v[2][1] - v[0][1], v[3][1] - v[0][1],
        v[0][1],  //
        v[1][2] - v[0][2], v[2][2] - v[0][2], v[3][2] - v[0][2],
        v[0][2],  //
        0, 0, 0, 1;
    const Eigen::Matrix4d Jinv = J.inverse();

    const auto baseInterpolationValues = VNCS::Sim3D::Tetra_4::interpolation(samplerDefinition.coordinates);
    const auto baseInterpolationPartialValues =
        VNCS::Sim3D::Tetra_4::interpolationPartial(samplerDefinition.coordinates);

    const Eigen::Vector3d basePoint = (J * (samplerDefinition.coordinates.homogeneous())).head<3>();

    VNCS::Sampler<VNCS::Space3D> sampler;
    sampler.x = {basePoint[0], basePoint[1], basePoint[2]};
    sampler.w = samplerDefinition.weight;
    sampler.a = 1.0;

    auto &baseShapeFunctions = sampler.coarseShapeFunctions;

    for (const auto &[vIndex, phi, dPhi] : ranges::views::zip(vIndices,  //
                                                              baseInterpolationValues,
                                                              baseInterpolationPartialValues)) {
        VNCS::ShapeFunction<VNCS::Space3D> shapeFunction;
        shapeFunction.nodeIndex = vIndex;
        shapeFunction.v = phi;
        const Eigen::Vector3d dv = Jinv.transpose().block<3, 3>(0, 0) * dPhi;
        shapeFunction.dv = {dv[0], dv[1], dv[2]};
        baseShapeFunctions.push_back(shapeFunction);
    }
    return sampler;
}
}  // namespace

void VNCS::Generators::Gen333::SamplingPointsExporter::operator()(const TetraMesh &coarseMesh,
                                                                  const TetraMesh &fineMesh)
{
    const auto logger = spdlog::get("Gen333");
    const VNCS::Generators::TetraMeshIntersection coarseIntersectionFunc(coarseMesh);
    const VNCS::Generators::TetraMeshIntersection fineIntersectionFunc(fineMesh);

    for (const auto &[coarseFaceId, coarseFace] : coarseMesh.tetras | ranges::views::enumerate) {
        // Get the indices and points of the face
        const auto vIndices = coarseFace;
        const auto v = vIndices |
                       ranges::views::transform([&coarseMesh](const auto &vd) { return coarseMesh.points[vd]; }) |
                       ranges::to_vector;

        const Space3D::Tetra coarseTetra = {v[0], v[1], v[2], v[3]};

        ::SamplerDefinition samplerDefinition;
        samplerDefinition.coordinates = {1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0};
        samplerDefinition.tetraIdx = coarseFaceId;
        samplerDefinition.weight = coarseTetra.volume();
        m_coarseSamplingPoints.push_back(generateSampler(samplerDefinition,  //
                                                         coarseMesh));
    }

    for (const auto [fineFaceId, fineFace] : fineMesh.tetras | ranges::views::enumerate) {
        const auto vIndices = fineFace;
        const auto v = vIndices | ranges::views::transform([&fineMesh](const auto &vd) { return fineMesh.points[vd]; }) |
                       ranges::to_vector;

        const Space3D::Tetra fineTriangle{v[0], v[1], v[2], v[3]};
        ::SamplerDefinition samplerDefinition;
        samplerDefinition.coordinates = {1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0};
        samplerDefinition.tetraIdx = fineFaceId;
        samplerDefinition.weight = fineTriangle.volume();
        auto sampler = generateSampler(samplerDefinition, fineMesh);
        sampler.fineShapeFunctions = sampler.coarseShapeFunctions;
        sampler.coarseShapeFunctions.clear();
        sampler.a = 0.0;
        m_fineSamplingPoints.push_back(sampler);
    }

    {
        json coarseSamplers = m_coarseSamplingPoints;
        std::ofstream oCoarse(m_coarseSamplersFilePath);
        oCoarse << std::setw(4) << coarseSamplers << std::endl;
    }

    {
        json fineSamplers = m_fineSamplingPoints;
        std::ofstream oFine(m_fineSamplersFilePath);
        oFine << std::setw(4) << fineSamplers << std::endl;
    }
}
