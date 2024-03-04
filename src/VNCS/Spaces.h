#ifndef VNCS_SPACES_H
#define VNCS_SPACES_H

#include <VNCS/Types.h>

#include <sofa/defaulttype/VecTypes.h>
#include <VNCS/DeformationGradientTypes.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Iterator_range.h>

#include <range/v3/view/adaptor.hpp>

template <typename T>
inline constexpr bool ranges::enable_borrowed_range<CGAL::Iterator_range<T>> = true;

template <typename T>
inline constexpr bool ranges::enable_view<CGAL::Iterator_range<T>> = true;

namespace VNCS
{
struct Space1D {
    constexpr static int dim = 1;
    using VecType = sofa::defaulttype::Vec1dTypes;
};

struct Space2D {
    constexpr static int dim = 2;
    using VecType = sofa::defaulttype::Vec2dTypes;
    using Real = VNCS::Real;

    using F22 = VNCS::F<Real, VNCS::Space2D, VNCS::Space2D>;
    using F21 = VNCS::F<Real, VNCS::Space2D, VNCS::Space1D>;

    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point = K::Point_2;
    using Mesh = CGAL::Surface_mesh<Point>;
    using Polyline = std::vector<Point>;
    using Segment = K::Segment_2;
    using Triangle = K::Triangle_2;
};

struct Space3D {
    constexpr static int dim = 3;
    using VecType = sofa::defaulttype::Vec3dTypes;
    using Real = VNCS::Real;

    using F33 = VNCS::F<Real, VNCS::Space3D, VNCS::Space3D>;
    using F32 = VNCS::F<Real, VNCS::Space3D, VNCS::Space2D>;
    using F31 = VNCS::F<Real, VNCS::Space3D, VNCS::Space1D>;

    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point = K::Point_3;
    using Mesh = CGAL::Surface_mesh<Point>;
    using Polyline = std::vector<Point>;
    using Segment = K::Segment_3;
    using Triangle = K::Triangle_3;
    using Tetra = K::Tetrahedron_3;

    struct TetraMesh {
        using TetraIndices = std::array<int, 4>;
        std::vector<VNCS::Space3D::Point> points;
        std::vector<TetraIndices> tetras;
    };
};
}  // namespace VNCS

#endif  // VNCS_SPACES_H
