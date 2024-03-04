#include "SamplingPointsExporter.h"

#include <spdlog/spdlog.h>
#include <spdlog/pattern_formatter.h>

#include <VNCS/Generators/Gen223Barycentric/GenerateSampler.h>
#include <VNCS/Generators/IntegrationScheme.h>

#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

#include <range/v3/view/filter.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/action/sort.hpp>
#include <range/v3/view/group_by.hpp>

void VNCS::Generators::Gen223Barycentric::SamplingPointsExporter::operator()(const VNCS::Space3D::Mesh &coarseMesh,
                                                                             const VNCS::Space3D::Mesh &fineMesh)
{
    GenerateSampler generateSampler;

    for (const auto &[coarseFaceId, coarseFace] : coarseMesh.faces() | ranges::views::enumerate) {
        m_coarseSamplingPoints.push_back(generateSampler(coarseFace,  //
                                                         coarseMesh));
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

std::filesystem::path VNCS::Generators::Gen223Barycentric::SamplingPointsExporter::coarseSamplersFilePath() const
{
    return m_coarseSamplersFilePath;
}

void VNCS::Generators::Gen223Barycentric::SamplingPointsExporter::setCoarseSamplersFilePath(
    const std::filesystem::path &coarseSamplersFilePath)
{
    m_coarseSamplersFilePath = coarseSamplersFilePath;
}

std::filesystem::path VNCS::Generators::Gen223Barycentric::SamplingPointsExporter::fineSamplersFilePath() const
{
    return m_fineSamplersFilePath;
}

void VNCS::Generators::Gen223Barycentric::SamplingPointsExporter::setFineSamplersFilePath(
    const std::filesystem::path &fineSamplersFilePath)
{
    m_fineSamplersFilePath = fineSamplersFilePath;
}
