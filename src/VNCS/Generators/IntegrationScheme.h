#ifndef VNCS_SIM2D_INTEGRATIONSCHEME_H
#define VNCS_SIM2D_INTEGRATIONSCHEME_H

#include <array>
#include <Eigen/Dense>
#include <VNCS/Spaces.h>

namespace VNCS
{
namespace Sim2D
{
namespace Tri_3
{
struct IntegrationPoint {
    Eigen::Vector2d point;
    VNCS::Space2D::Real weight;
};  // namespace IntegrationPoint

template <int IntegrationPoints>
using IntegrationScheme = std::array<IntegrationPoint, IntegrationPoints>;

inline std::array<VNCS::Space2D::Real, 3> interpolation(const Eigen::Vector2d &p)
{
    return {1 - p[0] - p[1], p[0], p[1]};
}

inline std::array<Eigen::Vector2d, 3> interpolationPartial(const Eigen::Vector2d &p)
{
    return {Eigen::Vector2d{-1, -1}, Eigen::Vector2d{1, 0}, Eigen::Vector2d{0, 1}};
}

namespace Gauss1
{
const static Eigen::Vector2d p0{1.0 / 3.0, 1.0 / 3.0};
const static IntegrationScheme<1> Scheme = {IntegrationPoint{p0, 1.0 / 2.0}};
}  // namespace Gauss1

namespace Gauss3
{
const static Eigen::Vector2d p0{1.0 / 4.0, 1.0 / 2.0};
const static Eigen::Vector2d p1{1.0 / 2.0, 1.0 / 4.0};
const static Eigen::Vector2d p2{1.0 / 4.0, 1.0 / 4.0};
const static IntegrationScheme<3> Scheme = {IntegrationPoint{p0, 1.0 / 6.0},
                                            IntegrationPoint{p1, 1.0 / 6.0},
                                            IntegrationPoint{p2, 1.0 / 6.0}};
}  // namespace Gauss3

namespace Gauss4
{
const static Eigen::Vector2d p0{1.0 / 3.0, 1.0 / 3.0};
const static Eigen::Vector2d p1{0.6, 0.2};
const static Eigen::Vector2d p2{0.2, 0.6};
const static Eigen::Vector2d p3{0.2, 0.2};
const static IntegrationScheme<4> Scheme = {IntegrationPoint{p0, -9.0 / 32.0},
                                            IntegrationPoint{p1, 25.0 / 96.0},
                                            IntegrationPoint{p2, 25.0 / 96.0},
                                            IntegrationPoint{p3, 25.0 / 96.0}};
}  // namespace Gauss4
}  // namespace Tri_3

namespace Segment_2
{
inline std::array<VNCS::Space2D::Real, 2> interpolation(const VNCS::Space2D::Real &p)
{
    return {1.0 / 2.0 - 1.0 / 2.0 * p, 1.0 / 2.0 + 1.0 / 2.0 * p};
}

inline std::array<VNCS::Space2D::Real, 2> interpolationPartial(const VNCS::Space2D::Real &p)
{
    return {-1.0 / 2.0, 1.0 / 2.0};
}

}  // namespace Segment_2
}  // namespace Sim2D

namespace Sim3D
{
struct IntegrationPoint {
    Eigen::Vector3d point;
    VNCS::Space2D::Real weight;
};  // namespace IntegrationPoint

namespace Tetra_4
{
template <int IntegrationPoints>
using IntegrationScheme = std::array<IntegrationPoint, IntegrationPoints>;

inline std::array<VNCS::Space3D::Real, 4> interpolation(const Eigen::Vector3d &p)
{
    return {1 - p[0] - p[1] - p[2], p[0], p[1], p[2]};
}

inline std::array<Eigen::Vector3d, 4> interpolationPartial(const Eigen::Vector3d &p)
{
    return {Eigen::Vector3d{-1, -1, -1}, Eigen::Vector3d{1, 0, 0}, Eigen::Vector3d{0, 1, 0}, Eigen::Vector3d{0, 0, 1}};
}

namespace Gauss1
{
const static Eigen::Vector3d p0{0.25, 0.25, 0.25};

const static IntegrationScheme<1> Scheme = {IntegrationPoint{p0, 1.0}};
}  // namespace Gauss1

namespace Gauss4
{
const static Eigen::Vector3d p0{0.1381966011250105, 0.1381966011250105, 0.1381966011250105};
const static Eigen::Vector3d p1{0.5854101966249685, 0.1381966011250105, 0.1381966011250105};
const static Eigen::Vector3d p2{0.1381966011250105, 0.5854101966249685, 0.1381966011250105};
const static Eigen::Vector3d p3{0.1381966011250105, 0.1381966011250105, 0.5854101966249685};

const static IntegrationScheme<4> Scheme = {IntegrationPoint{p0, 0.25},
                                            IntegrationPoint{p1, 0.25},
                                            IntegrationPoint{p2, 0.25},
                                            IntegrationPoint{p3, 0.25}};
}  // namespace Gauss4
}  // namespace Tetra_4
namespace Tri_3
{
struct IntegrationPoint {
    Eigen::Vector2d point;
    VNCS::Space2D::Real weight;
};  // namespace IntegrationPoint

template <int IntegrationPoints>
using IntegrationScheme = std::array<IntegrationPoint, IntegrationPoints>;

inline std::array<VNCS::Space2D::Real, 3> interpolation(const Eigen::Vector2d &p)
{
    return {1 - p[0] - p[1], p[0], p[1]};
}

inline std::array<Eigen::Vector2d, 3> interpolationPartial(const Eigen::Vector2d &p)
{
    return {Eigen::Vector2d{-1, -1}, Eigen::Vector2d{1, 0}, Eigen::Vector2d{0, 1}};
}

namespace Gauss1
{
const static Eigen::Vector2d p0{1.0 / 3.0, 1.0 / 3.0};
const static IntegrationScheme<1> Scheme = {IntegrationPoint{p0, 1.0}};
}  // namespace Gauss1

namespace Gauss3
{
const static Eigen::Vector2d p0{1.0 / 4.0, 1.0 / 2.0};
const static Eigen::Vector2d p1{1.0 / 2.0, 1.0 / 4.0};
const static Eigen::Vector2d p2{1.0 / 4.0, 1.0 / 4.0};
const static IntegrationScheme<3> Scheme = {IntegrationPoint{p0, 1.0 / 3.0},
                                            IntegrationPoint{p1, 1.0 / 3.0},
                                            IntegrationPoint{p2, 1.0 / 3.0}};
}  // namespace Gauss3

}  // namespace Tri_3
}  // namespace Sim3D
}  // namespace VNCS

#endif  // VNCS_INTEGRATIONSCHEME_H
