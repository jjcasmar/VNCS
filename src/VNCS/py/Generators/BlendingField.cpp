#include "BlendingField.h"
#include <pybind11/pybind11.h>
#include <spdlog/spdlog.h>
#include <pybind11/stl.h>
#include <memory>

#include <VNCS/Generators/BlendingField.h>
#include <VNCS/Spaces.h>

#include <VNCS/Generators/PMPTriangleIntersection.h>

#include <range/v3/view/transform.hpp>
#include <range/v3/to_container.hpp>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

namespace
{
template <typename WorldSpace>
class BlendingField : public VNCS::Generators::BlendingField<typename WorldSpace::Real, typename WorldSpace::Point>
{
public:
    BlendingField() = default;
    ~BlendingField() = default;

    typename WorldSpace::Real blending(const typename WorldSpace::Point &p) const final
    {
        using Inherit = VNCS::Generators::BlendingField<typename WorldSpace::Real, typename WorldSpace::Point>;
        std::array<typename WorldSpace::Real, WorldSpace::dim> pyp;
        for (int i = 0; i < WorldSpace::dim; ++i)
            pyp[i] = p[i];
        PYBIND11_OVERLOAD_PURE(typename WorldSpace::Real, /* Return type */
                               Inherit,                   /* Parent class */
                               blending,                  /* Name of function in C++ (must match Python name) */
                               pyp                        /* Argument(s) */
        );
    }
};

}  // namespace

void VNCS::Generators::py::blendingField(pybind11::module &m)
{
    using Inherit2 = VNCS::Generators::BlendingField<VNCS::Space2D::Real, VNCS::Space2D::Point>;
    pybind11::class_<Inherit2, ::BlendingField<VNCS::Space2D>, std::shared_ptr<Inherit2>>(m, "BlendingField2")
        .def(pybind11::init<>())
        .def("blending", &Inherit2::blending);

    using Inherit3 = VNCS::Generators::BlendingField<VNCS::Space3D::Real, VNCS::Space3D::Point>;
    pybind11::class_<Inherit3, ::BlendingField<VNCS::Space3D>, std::shared_ptr<Inherit3>>(m, "BlendingField3")
        .def(pybind11::init<>())
        .def("blending", &Inherit3::blending);

    struct PMPTriangleIntersectionBinding {
        PMPTriangleIntersectionBinding(VNCS::Space2D::Mesh mesh_)
            : mesh(mesh_)
            , intersection(mesh)
        {
        }

        PMPTriangleIntersectionBinding(const PMPTriangleIntersectionBinding &other)
            : mesh(other.mesh)
            , intersection(mesh)
        {
        }

        PMPTriangleIntersectionBinding(PMPTriangleIntersectionBinding &&other)
            : mesh(std::move(other.mesh))
            , intersection(mesh)
        {
        }

        VNCS::Space2D::Mesh mesh;
        VNCS::Generators::PMPTriangleIntersection intersection;
    };

    pybind11::class_<PMPTriangleIntersectionBinding>(m, "TriangleIntersection")
        .def(pybind11::init([](const std::vector<std::array<VNCS::Space2D::Real, 3>> &points,
                               const std::vector<std::array<int, 3>> &triangles) {
            const auto pointsCGAL =
                points | ranges::views::transform([](const auto &p) { return VNCS::Space2D::Point(p[0], p[1]); }) |
                ranges::to_vector;

            VNCS::Space2D::Mesh mesh;

            namespace PMP = CGAL::Polygon_mesh_processing;
            PMP::polygon_soup_to_polygon_mesh(pointsCGAL, triangles, mesh);

            PMPTriangleIntersectionBinding intersector(std::move(mesh));
            return intersector;
        }))
        .def("__call__",
             [](const PMPTriangleIntersectionBinding &intersector, const std::array<VNCS::Space2D::Real, 2> &p)
                 -> std::tuple<unsigned int, std::array<VNCS::Space2D::Real, 3>> {
                 const auto point = VNCS::Space2D::Point(p[0], p[1]);
                 const auto closestTriangle = intersector.intersection.query(point, true);
                 if (closestTriangle) {
                     const auto points =
                         intersector.mesh.vertices_around_face(intersector.mesh.halfedge(closestTriangle.value())) |
                         ranges::views::transform(
                             [&intersector](const auto vId) { return intersector.mesh.point(vId); }) |
                         ranges::to_vector;

                     using TriangleCoordinates =
                         CGAL::Barycentric_coordinates::Triangle_coordinates_2<VNCS::Space2D::K>;
                     TriangleCoordinates triangleCoordinates(points[0], points[1], points[2]);

                     std::array<VNCS::Space2D::Real, 3> coordinates;
                     triangleCoordinates(point, std::begin(coordinates));

                     return std::make_tuple(closestTriangle->idx(), coordinates);
                 }
                 return {};
             });
}
