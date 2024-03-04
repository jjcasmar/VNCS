#ifndef VNCS_COARSESTATESETTER_H
#define VNCS_COARSESTATESETTER_H

#include <filesystem>
#include <VNCS/Spaces.h>
#include <fmt/format.h>

#include <VNCS/Generators/PMPTriangleIntersection.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/IO/print_wavefront.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

namespace sofa
{
namespace component
{
namespace container
{
template <typename T>
class MechanicalObject;
}
}  // namespace component
}  // namespace sofa

namespace VNCS
{
template <typename T>
class CoarseStateSetter
{
    template <typename T>
    struct BarycentricCoordinates {
    };

    template <>
    struct BarycentricCoordinates<VNCS::Space2D::VecType> : public std::array<VNCS::Space2D::Real, 3> {
        VNCS::Space2D::Mesh::Face_index face;
    };

    template <>
    struct BarycentricCoordinates<VNCS::Space3D::VecType> : public std::array<VNCS::Space2D::Real, 4> {
    };

public:
    CoarseStateSetter();

    void setXNObject(sofa::component::container::MechanicalObject<T> *xNObject) { m_xNObject = xNObject; }
    void setXCObject(sofa::component::container::MechanicalObject<T> *xCObject) { m_xCObject = xCOject; }

    void init()
    {
        // Given the coarse simulation, the xN and the xC positions, I need to find the relation between each simulated
        // point (xN or xC) with the triangles in the coarse simulation

        // For that, we compute the barycentric coordinates of each xN or xC wrt to the elements of xN. Obviously, we
        // first need to know in which element each xN or xC lays.

        // We load the coarse simulation
        {
            std::vector<VNCS::Space2D::Point> coarsePoints;
            std::vector<std::vector<std::size_t>> coarseFaces;
            std::ifstream coarseMeshIn(coarseMeshPath);
            CGAL::read_OBJ(coarseMeshIn, coarsePoints, coarseFaces);

            namespace PMP = CGAL::Polygon_mesh_processing;
            PMP::polygon_soup_to_polygon_mesh(coarsePoints, coarseFaces, coarseInputMesh);
        }

        VNCS::Generators::PMPTriangleIntersection intersection()
    }

private:
    std::filesystem::path m_statePath;

    sofa::component::container::MechanicalObject<T> m_xNObject;
    sofa::component::container::MechanicalObject<T> m_xCObject;

    std::vector<BarycentricCoordinates<T>> m_xNBarycentricCoordinates;
    std::vector<BarycentricCoordinates<T>> m_xCBarycentricCoordinates;
};
}  // namespace VNCS

#endif  // COARSESTATESETTER_H
