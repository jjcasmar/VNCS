#include "GenerateSampler.h"

#include <VNCS/Generators/IntegrationScheme.h>
#include <CGAL/intersections.h>
#include <spdlog/spdlog.h>

#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>
#include <range/v3/to_container.hpp>

VNCS::Sampler<VNCS::Space2D> VNCS::Generators::Gen223::GenerateSampler::operator()(
    const VNCS::Space3D::Mesh::Face_index &faceHandle,
    const VNCS::Space3D::Mesh &fineMesh)
{
    // Get the indices and points of the face
    const auto vIndices = fineMesh.vertices_around_face(fineMesh.halfedge(faceHandle)) | ranges::to_vector;
    const auto v = vIndices | ranges::views::transform([&fineMesh](const auto &vd) { return fineMesh.point(vd); }) |
                   ranges::to_vector;

    VNCS::Space3D::Triangle t3d{v[0], v[1], v[2]};

    // Convert to a 2D triangle
    const auto L10squared = (v[1] - v[0]).squared_length();
    const auto L10 = std::sqrt(L10squared);
    const auto L20squared = (v[2] - v[0]).squared_length();
    const auto L21squared = (v[2] - v[1]).squared_length();

    const auto v0 = VNCS::Space2D::Point(0, 0);
    const auto v1 = VNCS::Space2D::Point(L10, 0);

    const auto v2x = (L10squared + L20squared - L21squared) / (2.0 * L10);
    const auto v2 = VNCS::Space2D::Point(v2x, std::sqrt(L20squared - v2x * v2x));
    VNCS::Space2D::Triangle t2d{v0, v1, v2};

    // Create the J matrix. This matrix allows to convert the coordinates of the sampler point to an actual position
    // and to compute derivatives
    Eigen::Matrix3d J;
    J << v1[0] - v0[0], v2[0] - v0[0], v0[0],  //
        v1[1] - v0[1], v2[1] - v0[1], v0[1],   //
        0, 0, 1;
    const Eigen::Matrix3d Jinv = J.inverse();

    const auto baseInterpolationValues = VNCS::Sim2D::Tri_3::interpolation({1.0 / 3.0, 1.0 / 3.0});
    const auto baseInterpolationPartialValues = VNCS::Sim2D::Tri_3::interpolationPartial({1.0 / 3.0, 1.0 / 3.0});

    VNCS::Sampler<VNCS::Space2D> sampler;
    sampler.w = t2d.area();
    sampler.a = 0.0;
    sampler.fineFaceIdx = faceHandle.idx();

    auto &shapeFunctions = sampler.fineShapeFunctions;

    for (const auto [vIndex, phi, dPhi] : ranges::views::zip(vIndices,  //
                                                             baseInterpolationValues,
                                                             baseInterpolationPartialValues)) {
        VNCS::ShapeFunction<VNCS::Space2D> shapeFunction;
        shapeFunction.nodeIndex = vIndex;
        shapeFunction.v = phi;
        const Eigen::Vector2d dv = Jinv.transpose().block<2, 2>(0, 0) * dPhi;
        shapeFunction.dv = {dv[0], dv[1]};
        shapeFunctions.push_back(shapeFunction);
    }
    return sampler;
}
