#include "SamplingPointsExporter.h"
#include <VNCS/DataExtensions.h>

#include <VNCS/IntegrationScheme.h>

#include <nlohmann/json.hpp>
#include <fstream>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>
#include <random>

#include "TetraAABBTree.h"

//#include <range/v3/view/iota.hpp>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/view/take.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/equal.hpp>
#include <range/v3/algorithm/equal_range.hpp>
#include <range/v3/view/group_by.hpp>
#include <range/v3/action/sort.hpp>

#include "KMeans.h"

#include <numeric>
#include <Eigen/Dense>

using namespace VNCS;
using json = nlohmann::json;

SamplingPointsExporter::SamplingPointsExporter()
    : sofa::core::objectmodel::BaseObject()
    , m_coarsePoints(initData(&m_coarsePoints, "coarsePoints", "Coarse points"))
    , m_finePoints(initData(&m_finePoints, "finePoints", "Fine points"))
    , m_lowResTetrahedra(initData(&m_lowResTetrahedra, "lowResTetrahedra", "Tetrahedra in the coarse region"))
    , m_highResTetrahedra(initData(&m_highResTetrahedra, "highResTetrahedra", "Tetrahedra in the fine region"))
{
}

void SamplingPointsExporter::init()
{
    // Sample a default tetrahedron
    const VNCS::Tetrahedron standardTetrahedron{VNCS::Point_3{0, 0, 0},  //
                                                VNCS::Point_3{1, 0, 0},
                                                VNCS::Point_3{0, 1, 0},
                                                VNCS::Point_3{0, 0, 1}};

    std::random_device rd;
    auto seed = 941263991;  // rd();
    spdlog::get("VNCS")->info("seed = {}", seed);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> r0Generator(0, 1.0);
    std::uniform_real_distribution<> r1Generator(0, 1.0);
    std::uniform_real_distribution<> r2Generator(0, 1.0);

    std::vector<VNCS::Point_3> mcSamplers;
    for (int i = 0; i < 10000; ++i) {
        const auto r0 = r0Generator(gen);
        const auto r1 = r1Generator(gen);
        const auto r2 = r2Generator(gen);
        auto p = VNCS::Point_3{r0, r1, r2};
        while (!standardTetrahedron.has_on_bounded_side(p)) {
            const auto r0 = r0Generator(gen);
            const auto r1 = r1Generator(gen);
            const auto r2 = r2Generator(gen);
            p = VNCS::Point_3{r0, r1, r2};
        }
        mcSamplers.push_back(p);
    }

    //    const std::vector<VNCS::Point_3> mcSamplers =
    //        ranges::views::iota(0) |
    //        ranges::views::transform([&gen, &r0Generator, &r1Generator, &r2Generator](const auto i) {
    //            const auto r0 = r0Generator(gen);
    //            const auto r1 = r1Generator(gen);
    //            const auto r2 = r2Generator(gen);
    //            return VNCS::Point_3(r0, r1, r2);
    //        }) |
    //        ranges::views::filter(
    //            [standardTetrahedron](const auto &p) { return standardTetrahedron.has_on_bounded_side(p); }) |
    //        ranges::views::take(10000) | ranges::to_vector;

    // Now I need to export information about integration points
    // The scheme is the following:
    /*
     * [
     *     {
     *         "point": [x,y,z] // Integration point in 3D. Actually this is unneeded
     *         "weight": 0.5 // Weight of the integration point
     *         "alpha": 0.5 // Blending value at the integration point
     *         "da": [dax,day,daz] // Alpha derivative at the integration point
     *         "shapeFunctions": {
     *             "coarse" [
     *                 "nodeIndex": 0 // Index of a coarse node
     *                 "v": 0.5 // Coarse shape function value of the ith coarse node at the integration point
     *                 "dv": [dvx, dvy, dvz] // Coarse shape function derivative of the ith coarse node at the
     * integration point
     *             ],
     *
     *             "fine": [
     *                 // At most
     *                 "nodeIndex": 0 // Index of a fine node
     *                 "v": 0.5 // Fine shape function value of the ith fine node at the integration point
     *                 "dv": [dvx, dvy, dvz] // Fine shape function derivative of the ith fine node at the
     * integration point
     *             ]
     *         }
     *     }
     * ]
     */

    const auto &lowResTetrahedra = make_read_accessor(m_lowResTetrahedra);
    const auto &coarseVertices = make_read_accessor(m_coarsePoints);

    const auto &highResTetrahedra = make_read_accessor(m_highResTetrahedra);
    const auto &fineVertices = make_read_accessor(m_finePoints);

    // Build an AABB tree with the coarse tetrahedra
    const auto sPoint2Point = [](const VNCS::SofaTypes::Point &p) { return VNCS::Point_3{p.x(), p.y(), p.z()}; };
    const std::vector<VNCS::Tetrahedron> coarseResCGALTetrahedra =
        lowResTetrahedra |
        ranges::views::transform([&coarseVertices, sPoint2Point](const VNCS::SofaTypes::Tetra &sTetra) {
            auto pA = sPoint2Point(coarseVertices[sTetra[0]]);
            auto pB = sPoint2Point(coarseVertices[sTetra[1]]);
            auto pC = sPoint2Point(coarseVertices[sTetra[2]]);
            auto pD = sPoint2Point(coarseVertices[sTetra[3]]);
            VNCS::Tetrahedron tetra{pA, pB, pC, pD};
            return tetra;
        }) |
        ranges::to_vector;

    TetraAABBTree coarseAABBTree(std::begin(coarseResCGALTetrahedra), std::end(coarseResCGALTetrahedra));
    coarseAABBTree.build();

    // Build an AABB tree with the fine tetrahedra
    const std::vector<VNCS::Tetrahedron> fineResCGALTetrahedra =
        highResTetrahedra |
        ranges::views::transform([&fineVertices, sPoint2Point](const VNCS::SofaTypes::Tetra &sTetra) {
            auto pA = sPoint2Point(fineVertices[sTetra[0]]);
            auto pB = sPoint2Point(fineVertices[sTetra[1]]);
            auto pC = sPoint2Point(fineVertices[sTetra[2]]);
            auto pD = sPoint2Point(fineVertices[sTetra[3]]);
            VNCS::Tetrahedron tetra{pA, pB, pC, pD};
            return tetra;
        }) |
        ranges::to_vector;

    TetraAABBTree fineAABBTree(std::begin(fineResCGALTetrahedra), std::end(fineResCGALTetrahedra));
    fineAABBTree.build();

    // Blending factor for each coarse and fine point
    const std::vector<double> coarseBlendingFactor = coarseVertices |
                                                     ranges::views::transform([this, sPoint2Point](const auto &p) {
                                                         return m_blendingField->blending(sPoint2Point(p));
                                                     }) |
                                                     ranges::to_vector;
    const std::vector<double> fineBlendingFactor = fineVertices |
                                                   ranges::views::transform([this, sPoint2Point](const auto &p) {
                                                       return m_blendingField->blending(sPoint2Point(p));
                                                   }) |
                                                   ranges::to_vector;

    {
        json coarseSamplers;
        for (const auto &tetra : lowResTetrahedra) {
            // Get the four vertices of the coarse tetrahedron and check if they are all in the coarse region

            // Check where are the four vertices
            const auto v =
                tetra | ranges::views::transform([&coarseVertices](const auto i) {
                    return Eigen::Vector3d(coarseVertices[i].x(), coarseVertices[i].y(), coarseVertices[i].z());
                }) |
                ranges::to_vector;

            const auto coarseBlendingFactorsForCoarseTetra =
                tetra |
                ranges::views::transform([&coarseBlendingFactor](const auto i) { return coarseBlendingFactor[i]; }) |
                ranges::to_vector;

            if (ranges::any_of(coarseBlendingFactorsForCoarseTetra, [](const auto &v) { return v == 1.0; })) {
                Eigen::Matrix4d J;
                J << v[1][0] - v[0][0], v[2][0] - v[0][0], v[3][0] - v[0][0], v[0][0],  //
                    v[1][1] - v[0][1], v[2][1] - v[0][1], v[3][1] - v[0][1], v[0][1],   //
                    v[1][2] - v[0][2], v[2][2] - v[0][2], v[3][2] - v[0][2], v[0][2],   //
                    0, 0, 0, 1;
                const Eigen::Matrix4d Jinv = J.inverse();

                for (const IntegrationPoint &integrationPoint : Tet_4_1) {
                    json sampler;
                    Eigen::Vector4d locationPoint = J * integrationPoint.point.homogeneous();

                    const auto pointsOutsideFineTetrahedra =
                        mcSamplers | ranges::views::transform([&J](const VNCS::Point_3 &mcSampler) {
                            const Eigen::Vector3d mcSamplerEigen{mcSampler[0], mcSampler[1], mcSampler[2]};
                            const Eigen::Vector3d pEigen = (J * mcSamplerEigen.homogeneous()).head<3>();
                            return VNCS::Point_3(pEigen[0], pEigen[1], pEigen[2]);
                        }) |
                        ranges::views::filter([&fineAABBTree](const VNCS::Point_3 &sampler) {
                            return !fineAABBTree.do_intersect(sampler);
                        }) |
                        ranges::to_vector;

                    if (pointsOutsideFineTetrahedra.size()) {
                        sampler["point"] = {locationPoint[0], locationPoint[1], locationPoint[2]};
                        sampler["weight"] = J.determinant() * integrationPoint.weight *
                                            static_cast<VNCS::Real>(pointsOutsideFineTetrahedra.size()) /
                                            static_cast<VNCS::Real>(mcSamplers.size());

                        sampler["alpha"] = 1.0;
                        sampler["da"] = {0, 0, 0};

                        sampler["shapeFunctions"] = json::object();
                        sampler["shapeFunctions"]["coarse"] = json::array();

                        const auto interpolationValues = Tet_4::interpolation(integrationPoint.point);
                        const auto interpolationPartialValues = Tet_4::interpolationPartial(integrationPoint.point);
                        for (int i = 0; i < 4; ++i) {
                            json coarseShapeFunction;
                            coarseShapeFunction["nodeIndex"] = tetra[i];
                            coarseShapeFunction["v"] = interpolationValues[i];
                            const Eigen::Vector3d dv =
                                (Jinv.transpose().block<3, 3>(0, 0) * interpolationPartialValues[i]);
                            coarseShapeFunction["dv"] = {dv[0], dv[1], dv[2]};
                            sampler["shapeFunctions"]["coarse"].push_back(coarseShapeFunction);
                        }

                        coarseSamplers.push_back(sampler);
                    }
                }
            }
        }
        std::ofstream oCoarse(m_coarseFilePath);
        oCoarse << std::setw(4) << coarseSamplers << std::endl;
    }

    {
        json fineSamplers;
        for (const auto &[tetraId, tetra] : ranges::views::enumerate(highResTetrahedra)) {
            spdlog::get("VNCS")->info("Generating samplers for fine tetrahedron #{}", tetraId);
            const auto v = tetra | ranges::views::transform([&fineVertices](const auto i) {
                               return Eigen::Vector3d(fineVertices[i].x(), fineVertices[i].y(), fineVertices[i].z());
                           }) |
                           ranges::to_vector;

            const auto fineBlendingFactorsForFineTetra =
                tetra |
                ranges::views::transform([&fineBlendingFactor](const auto i) { return fineBlendingFactor[i]; }) |
                ranges::to_vector;

            Eigen::Matrix4d J;
            J << v[1][0] - v[0][0], v[2][0] - v[0][0], v[3][0] - v[0][0], v[0][0],  //
                v[1][1] - v[0][1], v[2][1] - v[0][1], v[3][1] - v[0][1], v[0][1],   //
                v[1][2] - v[0][2], v[2][2] - v[0][2], v[3][2] - v[0][2], v[0][2],   //
                0, 0, 0, 1;

            Eigen::Matrix4d Jinv = J.inverse();

            if (ranges::all_of(fineBlendingFactorsForFineTetra, [](const auto v) { return v == 0.0; })) {
                spdlog::get("VNCS")->info(
                    "\tAll vertices have a blending value of 0. Consider it full and add only 1 sampler");
                for (const IntegrationPoint &integrationPoint : Tet_4_1) {
                    json sampler;
                    Eigen::Vector4d locationPoint = J * integrationPoint.point.homogeneous();

                    sampler["point"] = {locationPoint[0], locationPoint[1], locationPoint[2]};
                    sampler["weight"] = J.determinant() * integrationPoint.weight;

                    sampler["alpha"] = 0.0;
                    sampler["da"] = {0, 0, 0};

                    const auto interpolationValues = Tet_4::interpolation(integrationPoint.point);
                    const auto interpolationPartialValues = Tet_4::interpolationPartial(integrationPoint.point);

                    sampler["shapeFunctions"] = json::object();
                    sampler["shapeFunctions"]["fine"] = json::array();

                    for (unsigned int i = 0u; i < 4; ++i) {
                        json fineShapeFunction;
                        fineShapeFunction["nodeIndex"] = tetra[i];
                        fineShapeFunction["v"] = interpolationValues[i];
                        const Eigen::Vector3d dv = (Jinv.transpose().block<3, 3>(0, 0) * interpolationPartialValues[i]);
                        fineShapeFunction["dv"] = {dv[0], dv[1], dv[2]};
                        sampler["shapeFunctions"]["fine"].push_back(fineShapeFunction);
                    }

                    sampler["shapeFunctions"]["coarse"] = json::array();
                    fineSamplers.push_back(sampler);
                }
            } else {
                // This is an hybrid tetrahedron, so we need to check if the integration is normal scheme or k-means
                // If the 4 vertices of the tetrahedron lay in the same coarse tetrahedron, then is a normal tetrahedron
                auto coarseTetras =
                    v | ranges::views::transform([&coarseAABBTree](const auto &p) {
                        return coarseAABBTree.any_intersected_primitive(VNCS::Point_3{p[0], p[1], p[2]});
                    }) |
                    ranges::to_vector;

                if (std::equal(std::begin(coarseTetras) + 1, std::end(coarseTetras), std::begin(coarseTetras))) {
                    if (!coarseTetras[0].has_value()) {
                        spdlog::get("VNCS")->info(
                            "\tAll vertices are outside coarse tetrahedra. Dont add samplers on this fine tetrahedron "
                            "as it's not clear if its hybrid or full");
                    } else {
                        spdlog::get("VNCS")->info(
                            "\tAll vertices share the same coarse tetrahedra. Consider it hybrid and add 4 sampler");
                        // Compute J for the coarse tetrahedron
                        const auto &coarseTetraIterator = coarseTetras[0].value();
                        const auto &coarseTetra = *coarseTetraIterator;
                        auto coarseTetraId = std::distance(std::begin(coarseResCGALTetrahedra), coarseTetraIterator);
                        const auto &coarseSofaTetra = lowResTetrahedra[coarseTetraId];

                        Eigen::Matrix4d JCoarse;
                        JCoarse << coarseTetra[1][0] - coarseTetra[0][0], coarseTetra[2][0] - coarseTetra[0][0],
                            coarseTetra[3][0] - coarseTetra[0][0], coarseTetra[0][0],  //
                            coarseTetra[1][1] - coarseTetra[0][1], coarseTetra[2][1] - coarseTetra[0][1],
                            coarseTetra[3][1] - coarseTetra[0][1], coarseTetra[0][1],  //
                            coarseTetra[1][2] - coarseTetra[0][2], coarseTetra[2][2] - coarseTetra[0][2],
                            coarseTetra[3][2] - coarseTetra[0][2], coarseTetra[0][2],  //
                            0, 0, 0, 1;
                        // All the fine vertices are in the same tetrahedron. Add 4 samplers
                        for (const IntegrationPoint &integrationPoint : Tet_4_4) {
                            json sampler;
                            Eigen::Vector4d locationPoint = J * integrationPoint.point.homogeneous();

                            sampler["point"] = {locationPoint[0], locationPoint[1], locationPoint[2]};
                            sampler["weight"] = J.determinant() * integrationPoint.weight;

                            const auto interpolationValues = Tet_4::interpolation(integrationPoint.point);
                            const auto interpolationPartialValues = Tet_4::interpolationPartial(integrationPoint.point);

                            sampler["alpha"] = std::inner_product(std::begin(interpolationValues),
                                                                  std::end(interpolationValues),
                                                                  std::begin(fineBlendingFactorsForFineTetra),
                                                                  0.0);
                            const Eigen::Vector3d da = Jinv.transpose().block<3, 3>(0, 0) *
                                                       std::inner_product(std::begin(interpolationPartialValues),
                                                                          std::end(interpolationPartialValues),
                                                                          std::begin(fineBlendingFactorsForFineTetra),
                                                                          Eigen::Vector3d{0, 0, 0});
                            sampler["da"] = {da[0], da[1], da[2]};

                            sampler["shapeFunctions"] = json::object();
                            sampler["shapeFunctions"]["fine"] = json::array();

                            for (unsigned int i = 0u; i < 4; ++i) {
                                json fineShapeFunction;
                                fineShapeFunction["nodeIndex"] = tetra[i];
                                fineShapeFunction["v"] = interpolationValues[i];
                                const Eigen::Vector3d dv =
                                    (Jinv.transpose().block<3, 3>(0, 0) * interpolationPartialValues[i]);
                                fineShapeFunction["dv"] = {dv[0], dv[1], dv[2]};
                                sampler["shapeFunctions"]["fine"].push_back(fineShapeFunction);
                            }

                            const auto JCoarseInv = JCoarse.inverse();
                            const Eigen::Vector3d coarseIntegrationPoint = (JCoarseInv * locationPoint).head<3>();
                            const auto coarseInterpolationValues = Tet_4::interpolation(coarseIntegrationPoint);
                            const auto coarseInterpolationPartialValues =
                                Tet_4::interpolationPartial(coarseIntegrationPoint);
                            sampler["shapeFunctions"]["coarse"] = json::array();
                            for (unsigned int i = 0u; i < 4; ++i) {
                                json fineShapeFunction;
                                fineShapeFunction["nodeIndex"] = coarseSofaTetra[i];
                                fineShapeFunction["v"] = coarseInterpolationValues[i];
                                const Eigen::Vector3d dv =
                                    (JCoarseInv.transpose().block<3, 3>(0, 0) * coarseInterpolationPartialValues[i]);
                                fineShapeFunction["dv"] = {dv[0], dv[1], dv[2]};
                                sampler["shapeFunctions"]["coarse"].push_back(fineShapeFunction);
                            }
                            fineSamplers.push_back(sampler);
                        }
                    }
                } else {
                    // The fine tetrahedron is in several coarse tetrahedra. k-means!
                    spdlog::get("VNCS")->info(
                        "\tVertices are in different coarse tetrahedra. Compute regions and k-means");

                    struct PointAndTetraIdx {
                        VNCS::Point_3 point;
                        long tetrahedron;
                    };

                    auto sortedPointAndCoarseTetra =
                        mcSamplers | ranges::views::transform([](const auto &mcSampler) {
                            return Eigen::Vector3d(mcSampler[0], mcSampler[1], mcSampler[2]);
                        }) |
                        ranges::views::transform(
                            [&J](const auto &mcSampler) { return (J * mcSampler.homogeneous()).eval(); }) |
                        ranges::views::transform([](const auto &mcSampler) {
                            return VNCS::Point_3(mcSampler[0], mcSampler[1], mcSampler[2]);
                        }) |
                        ranges::views::transform([&coarseAABBTree](const auto &mcSampler) {
                            return std::make_tuple(mcSampler, coarseAABBTree.any_intersected_primitive(mcSampler));
                        }) |
                        ranges::views::filter([this](const auto &tetrahedron) -> bool {
                            return std::get<1>(tetrahedron).has_value() ||
                                   m_blendingField->blending(std::get<0>(tetrahedron)) == 0.0;
                        }) |
                        ranges::views::transform([&coarseResCGALTetrahedra](const auto &tetrahedron) {
                            return PointAndTetraIdx{
                                std::get<0>(tetrahedron),
                                std::distance(std::begin(coarseResCGALTetrahedra), std::get<1>(tetrahedron).value())};
                        }) |
                        ranges::to_vector | ranges::actions::sort(std::less<long>{}, &PointAndTetraIdx::tetrahedron);

                    auto groupedSortedPointCoarseTetra =
                        sortedPointAndCoarseTetra | ranges::views::group_by([](const auto &lhs, const auto &rhs) {
                            return lhs.tetrahedron == rhs.tetrahedron;
                        }) |
                        ranges::views::transform([](const auto &a) { return a | ranges::to_vector; }) |
                        ranges::to_vector;

                    spdlog::get("VNCS")->info("\t\tFound {} coarse tetrahedra for fine tetrahedron #{}",
                                              groupedSortedPointCoarseTetra.size(),
                                              tetraId);
                    // Compute k-means of the different sampled regions
                    for (const auto &[groupId, pointsForTetra] :
                         ranges::views::enumerate(groupedSortedPointCoarseTetra)) {
                        const std::vector<VNCS::Point_3> points =
                            pointsForTetra | ranges::views::transform(&PointAndTetraIdx::point) | ranges::to_vector;

                        if (points.size() >= 4) {
                            auto [centers, numberOfPointsPerRegion, pointsPerRegion, iterations] =
                                kMeans<4>(points, 500000);
                            const auto coarseTetraId = pointsForTetra[0].tetrahedron;
                            const auto &coarseTetra = coarseResCGALTetrahedra[coarseTetraId];
                            const auto &coarseSofaTetra = lowResTetrahedra[coarseTetraId];

                            spdlog::get("VNCS")->info("\t\tCompute k-means with {} clusters {}/{}",
                                                      4,
                                                      groupId + 1,
                                                      groupedSortedPointCoarseTetra.size());
                            for (const auto &[center, nbPoints] :
                                 ranges::views::zip(centers, numberOfPointsPerRegion)) {
                                spdlog::get("VNCS")->info("\t\t\tk-means region with {} vertices", nbPoints);
                                // Compute J for the coarse tetrahedron where this sampler is
                                json sampler;
                                sampler["point"] = {center[0], center[1], center[2]};
                                sampler["weight"] = J.determinant() * static_cast<VNCS::Real>(nbPoints) /
                                                    static_cast<VNCS::Real>(mcSamplers.size());

                                const Eigen::Vector3d centerInStandardTetrahedron =
                                    (Jinv * Eigen::Vector3d(center[0], center[1], center[2]).homogeneous()).head<3>();

                                const auto interpolationValues = Tet_4::interpolation(centerInStandardTetrahedron);
                                const auto interpolationPartialValues =
                                    Tet_4::interpolationPartial(centerInStandardTetrahedron);

                                sampler["alpha"] = std::inner_product(std::begin(interpolationValues),
                                                                      std::end(interpolationValues),
                                                                      std::begin(fineBlendingFactorsForFineTetra),
                                                                      0.0);
                                const Eigen::Vector3d da =
                                    Jinv.transpose().block<3, 3>(0, 0) *
                                    std::inner_product(std::begin(interpolationPartialValues),
                                                       std::end(interpolationPartialValues),
                                                       std::begin(fineBlendingFactorsForFineTetra),
                                                       Eigen::Vector3d{0, 0, 0});
                                sampler["da"] = {da[0], da[1], da[2]};

                                sampler["shapeFunctions"] = json::object();
                                sampler["shapeFunctions"]["fine"] = json::array();

                                for (unsigned int i = 0u; i < 4; ++i) {
                                    json fineShapeFunction;
                                    fineShapeFunction["nodeIndex"] = tetra[i];
                                    fineShapeFunction["v"] = interpolationValues[i];
                                    const Eigen::Vector3d dv =
                                        (Jinv.transpose().block<3, 3>(0, 0) * interpolationPartialValues[i]);
                                    fineShapeFunction["dv"] = {dv[0], dv[1], dv[2]};
                                    sampler["shapeFunctions"]["fine"].push_back(fineShapeFunction);
                                }

                                // Compute J for the coarse tetrahedron where this sampler is

                                Eigen::Matrix4d JCoarse;
                                JCoarse << coarseTetra[1][0] - coarseTetra[0][0], coarseTetra[2][0] - coarseTetra[0][0],
                                    coarseTetra[3][0] - coarseTetra[0][0], coarseTetra[0][0],  //
                                    coarseTetra[1][1] - coarseTetra[0][1], coarseTetra[2][1] - coarseTetra[0][1],
                                    coarseTetra[3][1] - coarseTetra[0][1], coarseTetra[0][1],  //
                                    coarseTetra[1][2] - coarseTetra[0][2], coarseTetra[2][2] - coarseTetra[0][2],
                                    coarseTetra[3][2] - coarseTetra[0][2], coarseTetra[0][2],  //
                                    0, 0, 0, 1;
                                const auto JCoarseInv = JCoarse.inverse();
                                const Eigen::Vector3d coarseIntegrationPoint =
                                    (JCoarseInv * Eigen::Vector3d(center[0], center[1], center[2]).homogeneous())
                                        .head<3>();
                                const auto coarseInterpolationValues = Tet_4::interpolation(coarseIntegrationPoint);
                                const auto coarseInterpolationPartialValues =
                                    Tet_4::interpolationPartial(coarseIntegrationPoint);
                                sampler["shapeFunctions"]["coarse"] = json::array();
                                for (unsigned int i = 0u; i < 4; ++i) {
                                    json coarseShapeFunction;
                                    coarseShapeFunction["nodeIndex"] = coarseSofaTetra[i];
                                    coarseShapeFunction["v"] = coarseInterpolationValues[i];
                                    const Eigen::Vector3d dv = (JCoarseInv.transpose().block<3, 3>(0, 0) *
                                                                coarseInterpolationPartialValues[i]);
                                    coarseShapeFunction["dv"] = {dv[0], dv[1], dv[2]};
                                    sampler["shapeFunctions"]["coarse"].push_back(coarseShapeFunction);
                                }
                                if (fineSamplers.size() == 14301) {
                                    spdlog::get("VNCS")->info("sampler = {}", sampler.dump(2));
                                }
                                fineSamplers.push_back(sampler);
                            }
                        }
                    }
                }
            }
        }
        // Write file
        std::ofstream oFine(m_fineFilePath);
        oFine << std::setw(4) << fineSamplers << std::endl;
    }
}

void SamplingPointsExporter::setCoarseFilePath(const std::filesystem::path &coarseFilePath)
{
    m_coarseFilePath = coarseFilePath;
}

void SamplingPointsExporter::setFineFilePath(const std::filesystem::path &fineFilePath)
{
    m_fineFilePath = fineFilePath;
}

void SamplingPointsExporter::setBlendingField(std::shared_ptr<BlendingField> blendingField)
{
    m_blendingField = blendingField;
}

std::shared_ptr<BlendingField> SamplingPointsExporter::blendingField() const
{
    return m_blendingField;
}
