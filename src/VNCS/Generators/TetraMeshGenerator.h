#ifndef VNCS_GENERATORS_TETRAMESHGENERATOR_H
#define VNCS_GENERATORS_TETRAMESHGENERATOR_H

#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <VNCS/Spaces.h>

namespace VNCS
{
namespace Generators
{
namespace detail
{
typedef typename CGAL::Mesh_3::Robust_intersection_traits_3<VNCS::Space3D::K> Geom_traits;
// typedef K Geom_traits;
typedef typename CGAL::Mesh_polyhedron_3<Geom_traits>::type Polyhedron;
typedef typename CGAL::Polyhedral_mesh_domain_with_features_3<Geom_traits, Polyhedron> Mesh_domain;
typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef typename CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
}  // namespace detail

struct C3t3Criteria {
    VNCS::Space3D::Real facet_angle;
    VNCS::Space3D::Real facet_size;
    VNCS::Space3D::Real facet_distance;
    VNCS::Space3D::Real cell_radius_edge_ratio;
    VNCS::Space3D::Real cell_size;

    VNCS::Space3D::Real sharpAngle;
    VNCS::Space3D::Real sharpEdgeSize;
};

struct RemeshCriteria {
    VNCS::Space3D::Real target_edge_length;
    VNCS::Space3D::Real sharpAngle;
    int iterations;
};

// C3t3 defines a terahedral mesh embedded in a triangulation
typedef typename CGAL::Mesh_complex_3_in_triangulation_3<VNCS::Generators::detail::Tr,
                                                         VNCS::Generators::detail::Mesh_domain::Corner_index,
                                                         VNCS::Generators::detail::Mesh_domain::Curve_segment_index>
    C3t3;

C3t3 makeC3T3FromMesh(const VNCS::Space3D::Mesh &mesh, C3t3Criteria criteria);

std::pair<std::vector<VNCS::Space3D::Point>, std::vector<std::array<int, 4>>> c3t3ToTetrahedra(const C3t3 &c3t3);

}  // namespace Generators
}  // namespace VNCS

#endif  // TETRAMESHGENERATOR_H
