#ifndef VNCS_SIM2D_KMEANSP_H
#define VNCS_SIM2D_KMEANSP_H

#include <numeric>
#include <range/v3/view/take.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/zip.hpp>
#include <range/v3/action/shuffle.hpp>
#include <range/v3/action/transform.hpp>
#include <range/v3/numeric/accumulate.hpp>

#include <spdlog/fmt/ostr.h>
#include <Eigen/Dense>

#include <spdlog/spdlog.h>

#include <omp.h>

namespace VNCS
{
namespace Generators
{
template <typename P, int k, typename Bootstrap, typename Closest, typename Center>
struct KMeans {
    struct KMeansResult {
        std::array<P, k> centers;
        std::array<std::vector<int>, k> associatedPoints;
        int iterations;
    };

    std::optional<KMeansResult> operator()(const std::vector<P> &points,
                                           int maxIter,
                                           Bootstrap bootstrap = Bootstrap(),
                                           Closest closest = Closest(),
                                           Center center = Center())
    {
        regions.clear();
        regions.resize(points.size());

        std::array<P, k> centers = bootstrap(points);

        auto converges = false;
        auto iter = 0;
        std::array<int, k> numberOfPointsPerCenter;

        while (!converges && iter < maxIter) {
            converges = true;

            for (const auto &[pointIndex, p] : ranges::views::enumerate(points)) {
                // Given the current centers and a point, computes the closest center
                auto centerIndex = closest(p, centers);
                converges = converges && regions[pointIndex] == centerIndex;
                regions[pointIndex] = centerIndex;
            }

            // Compute new centers
            centers = center(regions, points);

            iter++;
        }

        if (converges) {
            std::array<std::vector<int>, k> pointsPerCenter;
            for (const auto [pointId, region] : ranges::view::enumerate(regions)) {
                pointsPerCenter[region].push_back(pointId);
            }

            return KMeansResult{centers, pointsPerCenter, iter};
        }

        return {};
    }

private:
    std::vector<int> regions;
};

template <typename P, int k>
struct RandomKMeansBootstrap {
    std::array<P, k> operator()(const std::vector<P> &points) const
    {
        auto current = 0;
        auto insertInArray = [&current](auto currentArray, auto &&newValue) {
            currentArray[current] = std::move(newValue);
            current++;
            return currentArray;
        };

        std::vector<int> indices;
        indices.resize(points.size());
        std::iota(indices.begin(), indices.end(), 0);

        std::default_random_engine gen;
        indices |= ranges::action::shuffle(gen);

        const auto centers =
            ranges::accumulate(indices | ranges::views::take(k) |
                                   ranges::views::transform([&points](const auto index) { return points[index]; }),
                               std::array<P, k>{},
                               insertInArray);

        return centers;
    }
};

template <typename P, int k>
struct CGALCenter {
    std::array<P, k> operator()(const std::vector<int> &regions, const std::vector<P> &points)
    {
        std::array<P, k> centers;
        centers.fill(CGAL::ORIGIN);

        std::array<int, k> pointsPerCenter;
        pointsPerCenter.fill(0);
        for (int i = 0; i < points.size(); ++i) {
            const auto center = centers[regions[i]];
            const auto point = points[i];
            centers[regions[i]] = center + (point - CGAL::ORIGIN);
            pointsPerCenter[regions[i]]++;
        }

        for (int i = 0; i < k; ++i)
            centers[i] =
                CGAL::ORIGIN + (centers[i] - CGAL::ORIGIN) / static_cast<VNCS::Space3D::Real>(pointsPerCenter[i]);

        return centers;
    }
};
}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_SIM2D_PROJECTION_H
