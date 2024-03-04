#ifndef VNCS_INTEGRATIONSCHEME_H
#define VNCS_INTEGRATIONSCHEME_H

#include <array>
#include <Eigen/Dense>

struct Tet_4 {
    static std::array<double, 4> interpolation(const Eigen::Vector3d &p)
    {
        return {1 - p[0] - p[1] - p[2], p[0], p[1], p[2]};
    }

    static std::array<Eigen::Vector3d, 4> interpolationPartial(const Eigen::Vector3d &p)
    {
        return {
            Eigen::Vector3d{-1, -1, -1}, Eigen::Vector3d{1, 0, 0}, Eigen::Vector3d{0, 1, 0}, Eigen::Vector3d{0, 0, 1}};
    }
};

struct Tri_3 {
    static std::array<double, 3> interpolation(const Eigen::Vector2d &p) { return {1 - p[0] - p[1], p[0], p[1]}; }

    static std::array<Eigen::Vector2d, 3> interpolationPartial(const Eigen::Vector2d &p)
    {
        return {Eigen::Vector2d{-1, -1}, Eigen::Vector2d{1, 0}, Eigen::Vector2d{0, 1}};
    }
};

struct IntegrationPoint {
    Eigen::Vector3d point;
    double weight;
};

template <int IntegrationPoints>
using IntegrationScheme = std::array<IntegrationPoint, IntegrationPoints>;

const static Eigen::Vector3d p0{0.5854101966249685, 0.1381966011250105, 0.1381966011250105};
const static Eigen::Vector3d p1{0.1381966011250105, 0.1381966011250105, 0.1381966011250105};
const static Eigen::Vector3d p2{0.1381966011250105, 0.1381966011250105, 0.5854101966249685};
const static Eigen::Vector3d p3{0.1381966011250105, 0.5854101966249685, 0.1381966011250105};

const static Eigen::Vector3d p00{1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0};

const static IntegrationScheme<4> Tet_4_4 = {IntegrationPoint{p0, 0.25 / 6.0},
                                             IntegrationPoint{p1, 0.25 / 6.0},
                                             IntegrationPoint{p2, 0.25 / 6.0},
                                             IntegrationPoint{p3, 0.25 / 6.0}};

const static IntegrationScheme<1> Tet_4_1 = {IntegrationPoint{p00, 1.0 / 6.0}};
const static IntegrationScheme<1> Tet_3_1 = {IntegrationPoint{p00, 1.0 / 6.0}};
#endif  // VNCS_INTEGRATIONSCHEME_H
