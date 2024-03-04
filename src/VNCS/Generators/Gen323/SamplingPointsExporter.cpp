#include "SamplingPointsExporter.h"

#include <spdlog/spdlog.h>
#include <spdlog/pattern_formatter.h>

#include <VNCS/Generators/TetraMeshIntersection.h>
#include <VNCS/Generators/Gen323/GenerateSampler.h>
#include <VNCS/Generators/KMeans.h>
#include <VNCS/Generators/IntegrationScheme.h>

#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/action/sort.hpp>
#include <range/v3/view/group_by.hpp>

void VNCS::Generators::Gen323::SamplingPointsExporter::operator()(const TetraMesh &coarseMesh,
                                                                  const VNCS::Space3D::Mesh &fineMesh)
{
    GenerateSampler generateSampler;
    for (const auto &[coarseTetraId, coarseTetra] : coarseMesh.tetras | ranges::views::enumerate) {
        // Get the indices and points of the face
        const auto vIndices = coarseTetra;
        const auto v = vIndices |
                       ranges::views::transform([&coarseMesh](const auto &vd) { return coarseMesh.points[vd]; }) |
                       ranges::to_vector;

        const VNCS::Space3D::Tetra coarseTetrahedron(v[0], v[1], v[2], v[3]);

        for (const auto &sample : VNCS::Sim3D::Tetra_4::Gauss1::Scheme) {
            CoarseSamplerDefinition samplerDefinition;
            samplerDefinition.coordinates = sample.point;
            samplerDefinition.baseFace = coarseTetraId;
            samplerDefinition.weight = sample.weight * coarseTetrahedron.volume();
            m_coarseSamplingPoints.push_back(generateSampler(samplerDefinition,  //
                                                             coarseMesh));
        }
    }

    for (const auto &[fineFaceId, fineFace] : fineMesh.faces() | ranges::views::enumerate) {
        m_fineSamplingPoints.push_back(generateSampler(fineFace,  //
                                                       fineMesh));
    }

    json coarseSamplers = m_coarseSamplingPoints;
    json fineSamplers = m_fineSamplingPoints;

    std::ofstream oCoarse(m_coarseSamplersFilePath);
    oCoarse << std::setw(4) << coarseSamplers << std::endl;

    std::ofstream oNodalFine(m_fineSamplersFilePath);
    oNodalFine << std::setw(4) << fineSamplers << std::endl;
}

std::filesystem::path VNCS::Generators::Gen323::SamplingPointsExporter::coarseSamplersFilePath() const
{
    return m_coarseSamplersFilePath;
}

void VNCS::Generators::Gen323::SamplingPointsExporter::setCoarseSamplersFilePath(
    const std::filesystem::path &coarseSamplersFilePath)
{
    m_coarseSamplersFilePath = coarseSamplersFilePath;
}

std::filesystem::path VNCS::Generators::Gen323::SamplingPointsExporter::fineSamplersFilePath() const
{
    return m_fineSamplersFilePath;
}

void VNCS::Generators::Gen323::SamplingPointsExporter::setFineSamplersFilePath(
    const std::filesystem::path &fineSamplersFilePath)
{
    m_fineSamplersFilePath = fineSamplersFilePath;
}
