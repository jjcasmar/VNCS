#include "SamplingPointsExporter.h"

#include <VNCS/Generators/PMPTriangleIntersection.h>
#include <VNCS/Generators/Gen212/GenerateSampler.h>

#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/view/subrange.hpp>
#include <VNCS/EdgeMesh.h>

namespace VNCS
{
namespace Generators
{
namespace Gen212
{
void SamplingPointsExporter::operator()(const VNCS::Space2D::Mesh &coarseMesh, const VNCS::EdgeMesh<VNCS::Space2D> &fineMesh)
{
    GenerateSampler generateSampler;
    for (const auto [coarseFaceId, coarseFace] : coarseMesh.faces() | ranges::views::enumerate) {
        // Get the indices and points of the face
        const auto v = coarseMesh.vertices_around_face(coarseMesh.halfedge(coarseFace)) |
                       ranges::views::transform([&coarseMesh](const auto &vd) { return coarseMesh.point(vd); }) |
                       ranges::to_vector;

        const VNCS::Space2D::Triangle coarseTriangle = {v[0], v[1], v[2]};

        CoarseSamplerDefinition samplerDefinition;
        samplerDefinition.coordinates = {1.0 / 3.0, 1.0 / 3.0};
        samplerDefinition.baseFace = coarseFace;
        samplerDefinition.weight = coarseTriangle.area();
        m_coarseSamplingPoints.push_back(generateSampler(samplerDefinition,  //
                                                         coarseMesh));
    }

    for (const auto &vertexHandle : boost::make_iterator_range(boost::vertices(fineMesh))) {

        FineSamplerDefinition samplerDefinition;
        samplerDefinition.pointHandle = vertexHandle;
        samplerDefinition.weight = 0;

        const auto &pointsMap = boost::get(PointVertexTag{}, fineMesh);
        const auto &point0 = boost::get(pointsMap, vertexHandle);

        for (const auto &edge : boost::make_iterator_range(boost::out_edges(vertexHandle, fineMesh)))
        {
            const auto targetVertexHandle = boost::target(edge, fineMesh);
            const auto &point1 = boost::get(pointsMap, targetVertexHandle);
            const VNCS::Space2D::Segment segment = VNCS::Space2D::Segment{point0, point1};
            samplerDefinition.weight += 1.0 / 2.0 * std::sqrt(segment.squared_length());
        }

        m_fineSamplingPoints.push_back(generateSampler(samplerDefinition));
    }

    json coarseSamplers = m_coarseSamplingPoints;
    json fineSamplers = m_fineSamplingPoints;

    std::ofstream oCoarse(m_coarseSamplersFilePath);
    oCoarse << std::setw(4) << coarseSamplers << std::endl;

    std::ofstream oFineNode(m_fineSamplersFilePath);
    oFineNode << std::setw(4) << fineSamplers << std::endl;
}

std::filesystem::path SamplingPointsExporter::coarseSamplersFilePath() const
{
    return m_coarseSamplersFilePath;
}

void SamplingPointsExporter::setCoarseSamplersFilePath(const std::filesystem::path &coarseSamplersFilePath)
{
    m_coarseSamplersFilePath = coarseSamplersFilePath;
}

std::filesystem::path SamplingPointsExporter::fineSamplersFilePath() const
{
    return m_fineSamplersFilePath;
}

void SamplingPointsExporter::setFineSamplersFilePath(const std::filesystem::path &fineSamplersFilePath)
{
    m_fineSamplersFilePath = fineSamplersFilePath;
}
}  // namespace Gen212
}  // namespace Generators
}  // namespace VNCS
