#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/constructions/kernel_ftC3.h>
#include <pybind11/pybind11.h>

#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <VNCS/Spaces.h>
#include <VNCS/SamplingPoints.h>

#include <VNCS/py/SamplingPoints.h>
#include <VNCS/py/ConjugateGradient.h>
#include <VNCS/py/MechanicalOperations.h>

#include <VNCS/py/Sim2D/MechanicalObject.h>
#include <VNCS/py/Sim2D/Projection.h>
#include <VNCS/py/Sim2D/UfMap.h>
#include <VNCS/py/Sim2D/ConcatMap.h>
#include <VNCS/py/Sim2D/DisplacementMap.h>
#include <VNCS/py/Sim2D/Mass.h>
#include <VNCS/py/Sim2D/DeformationGradientMO.h>
#include <VNCS/py/Sim2D/DeformationGradientMap.h>
#include <VNCS/py/Sim2D/StVK.h>
#include <VNCS/py/Sim2D/BendingEnergy.h>
#include <VNCS/py/Sim2D/StretchEnergy.h>
#include <VNCS/py/Sim2D/FinalPositionMap.h>
#include <VNCS/py/Sim2D/Vec2D23D.h>
#include <VNCS/py/Sim2D/EnhaceRelationMap.h>

#include <VNCS/py/Sim3D/MechanicalObject.h>
#include <VNCS/py/Sim3D/Projection.h>
#include <VNCS/py/Sim3D/UfMap.h>
#include <VNCS/py/Sim3D/ConcatMap.h>
#include <VNCS/py/Sim3D/DisplacementMap.h>
#include <VNCS/py/Sim3D/Mass.h>
#include <VNCS/py/Sim3D/DeformationGradientMO.h>
#include <VNCS/py/Sim3D/DeformationGradientMap.h>
#include <VNCS/py/Sim3D/StVK.h>
#include <VNCS/py/Sim3D/BendingEnergy.h>
#include <VNCS/py/Sim3D/StretchEnergy.h>
#include <VNCS/py/Sim3D/FinalPositionMap.h>

#include <VNCS/py/Generators/AdaptiveMeshCriteriaFunctor.h>
#include <VNCS/py/Generators/BlendingField.h>
#include <VNCS/py/Generators/Gen222.h>
#include <VNCS/py/Generators/Gen212.h>
#include <VNCS/py/Generators/Gen202.h>
#include <VNCS/py/Generators/Gen333.h>
#include <VNCS/py/Generators/Gen323.h>
#include <VNCS/py/Generators/Gen223.h>
#include <VNCS/py/Generators/Gen223Barycentric.h>

#include <sofa/core/behavior/MultiVec.h>
#include <sofa/core/ExecParams.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <SofaEigen2Solver/EigenVectorWrapper.h>
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>

#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <spdlog_setup/conf.h>
#include <filesystem>

#include <CGAL/Exact_spherical_kernel_3.h>

using Spherical_K = CGAL::Exact_spherical_kernel_3;
using Point_3 = CGAL::Point_3<Spherical_K>;
using Sphere_3 = CGAL::Sphere_3<Spherical_K>;

PYBIND11_MODULE(PyVNCS, m)
{
    m.doc() = "VNCS python bindings";
    m.attr("__version__") = "dev";

    if (std::filesystem::exists("spdlog_setup.toml"))
        spdlog_setup::from_file("spdlog_setup.toml");

    Eigen::setNbThreads(omp_get_max_threads() / 2);

    auto module2D = m.def_submodule("Sim2D", "2D bindings for VNCS");
    auto module3D = m.def_submodule("Sim3D", "3D bindings for VNCS");
    auto moduleGenerators = m.def_submodule("Generators", "Generators for VNCS");

    VNCS::py::samplingPoints(m);
    VNCS::py::conjugateGradient(m);
    VNCS::py::mechanicalOperations(m);

    VNCS::Generators::py::blendingField(moduleGenerators);
    VNCS::Generators::py::gen222(moduleGenerators);
    VNCS::Generators::py::gen212(moduleGenerators);
    VNCS::Generators::py::gen202(moduleGenerators);
    VNCS::Generators::py::gen333(moduleGenerators);
    VNCS::Generators::py::gen323(moduleGenerators);
    VNCS::Generators::py::gen223(moduleGenerators);
    VNCS::Generators::py::gen223Barycentric(moduleGenerators);
    VNCS::Generators::py::adaptiveMeshCriteriaFunctor(moduleGenerators);

    VNCS::Sim2D::py::mechanicalObject(module2D);
    VNCS::Sim2D::py::projection(module2D);
    VNCS::Sim2D::py::ufMap(module2D);
    VNCS::Sim2D::py::concatMap(module2D);
    VNCS::Sim2D::py::displacementMap(module2D);
    VNCS::Sim2D::py::mass(module2D);
    VNCS::Sim2D::py::deformationGradientMO(module2D);
    VNCS::Sim2D::py::deformationGradientMap(module2D);
    VNCS::Sim2D::py::stVK(module2D);
    VNCS::Sim2D::py::bendingEnergy(module2D);
    VNCS::Sim2D::py::stretchEnergy(module2D);
    VNCS::Sim2D::py::finalPositionMap(module2D);
    VNCS::Sim2D::py::vec2D23D(module2D);
    VNCS::Sim2D::py::enhaceMap(module2D);

    VNCS::Sim3D::py::mechanicalObject(module3D);
    VNCS::Sim3D::py::projection(module3D);
    VNCS::Sim3D::py::ufMap(module3D);
    VNCS::Sim3D::py::concatMap(module3D);
    VNCS::Sim3D::py::displacementMap(module3D);
    VNCS::Sim3D::py::mass(module3D);
    VNCS::Sim3D::py::deformationGradientMO(module3D);
    VNCS::Sim3D::py::deformationGradientMap(module3D);
    VNCS::Sim3D::py::stVK(module3D);
    VNCS::Sim3D::py::stretchEnergy(module3D);
    VNCS::Sim3D::py::bendingEnergy(module3D);
    VNCS::Sim3D::py::finalPositionMap(module3D);

    // We use this for the hand demo
    m.def("pointsInsideMesh", [](std::vector<std::array<VNCS::Space3D::Real, 3>> points, std::string meshPath) {
        VNCS::Space3D::Mesh mesh;

        {
            std::vector<VNCS::Space3D::Point> meshPoints;
            std::vector<std::vector<std::size_t>> meshFaces;
            std::ifstream coarseMeshIn(meshPath);
            CGAL::read_OBJ(coarseMeshIn, meshPoints, meshFaces);

            namespace PMP = CGAL::Polygon_mesh_processing;
            PMP::polygon_soup_to_polygon_mesh(meshPoints, meshFaces, mesh);
        }

        std::vector<VNCS::Space3D::Mesh> ccMeshes;

        CGAL::Polygon_mesh_processing::split_connected_components(mesh, ccMeshes);

        std::vector<std::pair<int, int>> pointsInside;
        for (int i = 0; i < ccMeshes.size(); ++i) {
            CGAL::Side_of_triangle_mesh<decltype(mesh), VNCS::Space3D::K> inside(mesh);

            for (int j = 0; j < points.size(); ++j) {
                VNCS::Space3D::Point p(points[i][0], points[i][1], points[i][2]);
                if (inside(p) == CGAL::ON_BOUNDED_SIDE)
                    pointsInside.push_back(std::pair(j, i));
            }
        }
        return pointsInside;
    });

    // m.def("intersection",
    //       [](const Eigen::Vector3d p0,
    //          double d0,
    //          const Eigen::Vector3d p1,
    //          double d1,
    //          const Eigen::Vector3d p2,
    //          double d2) {
    //           Sphere_3 s0(Point_3(p0[0], p0[1], p0[2]), d0);
    //           Sphere_3 s1(Point_3(p1[0], p1[1], p1[2]), d1);
    //           Sphere_3 s2(Point_3(p2[0], p2[1], p2[2]), d2);

    //           std::vector<Point_3> intersecs;
    //           CGAL::intersection(s0, s1, s2, std::back_inserter(intersecs));
    //           return intersecs;
    //       });

    m.def("isPointInsideMesh",
          [](std::array<VNCS::Space3D::Real, 3> point,
             const std::vector<std::array<VNCS::Space3D::Real, 3>> &points,
             std::vector<std::array<int, 3>> &triangles) {
              VNCS::Space3D::Mesh mesh;

              {
                  namespace PMP = CGAL::Polygon_mesh_processing;
                  PMP::polygon_soup_to_polygon_mesh(points, triangles, mesh);
              }

              CGAL::Side_of_triangle_mesh<decltype(mesh), VNCS::Space3D::K> inside(mesh);
              VNCS::Space3D::Point p(point[0], point[1], point[2]);
              if (inside(p) == CGAL::ON_BOUNDED_SIDE)
                  return true;
              return false;
          });

    m.def("computeBarycentricMatrix2",
          [](std::vector<std::array<VNCS::Real, 2>> points,
             std::vector<std::array<VNCS::Real, 3>> triMeshPoints,
             std::vector<std::array<int, 3>> tris) {
              std::vector<Eigen::Triplet<VNCS::Real>> triplets;

              const auto toCGAL3D = [](const std::array<VNCS::Real, 3> &v) { return VNCS::Space2D::Point(v[0], v[1]); };
              const auto toCGAL2D = [](const std::array<VNCS::Real, 2> &v) { return VNCS::Space2D::Point(v[0], v[1]); };
              const auto toCGALTri = [&triMeshPoints, toCGAL2D, toCGAL3D](const std::array<int, 3> &v) {
                  VNCS::Space2D::Triangle t(toCGAL3D(triMeshPoints[v[0]]),  //
                                            toCGAL3D(triMeshPoints[v[1]]),
                                            toCGAL3D(triMeshPoints[v[2]]));
                  return t;
              };

              std::vector<VNCS::Space2D::Triangle> triangles;
              for (const auto &tri : tris) {
                  const auto t = toCGALTri(tri);
                  triangles.push_back(t);
              }

              for (int pId = 0; pId < points.size(); ++pId) {
                  VNCS::Space2D::Point p = toCGAL2D(points[pId]);
                  for (int i = 0; i < triangles.size(); ++i) {
                      const auto &triangle = triangles[i];
                      if (!triangle.has_on_unbounded_side(p)) {
                          const auto &indices = tris[i];
                          const auto p0 = toCGAL3D(triMeshPoints[indices[0]]);
                          const auto p1 = toCGAL3D(triMeshPoints[indices[1]]);
                          const auto p2 = toCGAL3D(triMeshPoints[indices[2]]);

                          Eigen::Matrix3d A({{p0[0], p1[0], p2[0]},  //
                                             {p0[1], p1[1], p2[1]},
                                             {1.0, 1.0, 1.0}});
                          Eigen::Vector3d b = A.inverse() * Eigen::Vector3d({p[0], p[1], 1.0});
                          triplets.emplace_back(2 * pId + 0, 2 * indices[0] + 0, b[0]);
                          triplets.emplace_back(2 * pId + 0, 2 * indices[1] + 0, b[1]);
                          triplets.emplace_back(2 * pId + 0, 2 * indices[2] + 0, b[2]);

                          triplets.emplace_back(2 * pId + 1, 2 * indices[0] + 1, b[0]);
                          triplets.emplace_back(2 * pId + 1, 2 * indices[1] + 1, b[1]);
                          triplets.emplace_back(2 * pId + 1, 2 * indices[2] + 1, b[2]);
                          break;
                      }
                  }
              }
              Eigen::SparseMatrix<VNCS::Real> matrix;
              matrix.resize(2 * points.size(), 2 * triMeshPoints.size());
              matrix.setFromTriplets(std::begin(triplets), std::end(triplets));
              return matrix;
          });

    m.def("computeBarycentricMatrix",
          [](std::vector<std::array<VNCS::Real, 3>> points,
             std::vector<std::array<VNCS::Real, 3>> tetraMeshPoints,
             std::vector<std::array<int, 4>> tetras) {
              std::vector<Eigen::Triplet<VNCS::Real>> triplets;

              const auto toCGAL = [](const std::array<VNCS::Real, 3> &v) {
                  return VNCS::Space3D::Point(v[0], v[1], v[2]);
              };
              const auto toCGALTetra = [&tetraMeshPoints, toCGAL](const std::array<int, 4> &v) {
                  VNCS::Space3D::Tetra t(toCGAL(tetraMeshPoints[v[0]]),
                                         toCGAL(tetraMeshPoints[v[1]]),
                                         toCGAL(tetraMeshPoints[v[2]]),
                                         toCGAL(tetraMeshPoints[v[3]]));
                  return t;
              };

              std::vector<VNCS::Space3D::Tetra> tetrahedra;
              for (const auto &tetra : tetras) {
                  const VNCS::Space3D::Tetra t = toCGALTetra(tetra);
                  tetrahedra.push_back(t);
              }

              for (int pId = 0; pId < points.size(); ++pId) {
                  VNCS::Space3D::Point p = toCGAL(points[pId]);
                  for (int i = 0; i < tetrahedra.size(); ++i) {
                      const auto &tetra = tetrahedra[i];
                      if (!tetra.has_on_unbounded_side(p)) {
                          const auto &indices = tetras[i];
                          const auto p0 = toCGAL(tetraMeshPoints[indices[0]]);
                          const auto p1 = toCGAL(tetraMeshPoints[indices[1]]);
                          const auto p2 = toCGAL(tetraMeshPoints[indices[2]]);
                          const auto p3 = toCGAL(tetraMeshPoints[indices[3]]);

                          Eigen::Matrix4d A({{p0[0], p1[0], p2[0], p3[0]},
                                             {p0[1], p1[1], p2[1], p3[1]},
                                             {p0[2], p1[2], p2[2], p3[2]},
                                             {1.0, 1.0, 1.0, 1.0}});
                          Eigen::Vector4d b = A.inverse() * Eigen::Vector4d({p[0], p[1], p[2], 1.0});
                          triplets.emplace_back(3 * pId + 0, 3 * indices[0] + 0, b[0]);
                          triplets.emplace_back(3 * pId + 0, 3 * indices[1] + 0, b[1]);
                          triplets.emplace_back(3 * pId + 0, 3 * indices[2] + 0, b[2]);
                          triplets.emplace_back(3 * pId + 0, 3 * indices[3] + 0, b[3]);

                          triplets.emplace_back(3 * pId + 1, 3 * indices[0] + 1, b[0]);
                          triplets.emplace_back(3 * pId + 1, 3 * indices[1] + 1, b[1]);
                          triplets.emplace_back(3 * pId + 1, 3 * indices[2] + 1, b[2]);
                          triplets.emplace_back(3 * pId + 1, 3 * indices[3] + 1, b[3]);

                          triplets.emplace_back(3 * pId + 2, 3 * indices[0] + 2, b[0]);
                          triplets.emplace_back(3 * pId + 2, 3 * indices[1] + 2, b[1]);
                          triplets.emplace_back(3 * pId + 2, 3 * indices[2] + 2, b[2]);
                          triplets.emplace_back(3 * pId + 2, 3 * indices[3] + 2, b[3]);
                          break;
                      }
                  }
              }
              Eigen::SparseMatrix<VNCS::Real> matrix;
              matrix.resize(3 * points.size(), 3 * tetraMeshPoints.size());
              matrix.setFromTriplets(std::begin(triplets), std::end(triplets));
              return matrix;
          });

    m.def("computeForce", [](Eigen::VectorXd p, sofa::core::objectmodel::BaseContext *ctx) {
        sofa::core::MechanicalParams params;
        sofa::simulation::common::MechanicalOperations mop(&params, ctx);
        sofa::simulation::common::VectorOperations vop(&params, ctx);

        sofa::core::behavior::MultiVecCoord pos(&vop, sofa::core::VecCoordId::position());
        sofa::core::behavior::MultiVecDeriv vel(&vop, sofa::core::VecDerivId::velocity());
        sofa::core::behavior::MultiVecDeriv f(&vop, sofa::core::VecDerivId::force());
        sofa::core::behavior::MultiVecDeriv dx(&vop, sofa::core::VecDerivId::dx());
        dx.realloc(&vop, false, true);

        vel.clear();
        f.clear();

        // Assign the p value
        sofa::component::linearsolver::DefaultMultiMatrixAccessor accessor;
        mop.getMatrixDimension(nullptr, nullptr, &accessor);
        const auto n = static_cast<Eigen::Index>(accessor.getGlobalDimension());

        sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> rEigenWrapped(p);
        mop.baseVector2MultiVector(&rEigenWrapped, pos.id(), &accessor);

        mop.propagateXAndV(pos, vel);
        mop.computeForce(f);

        sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> pEigenWrapped(p);
        mop.multiVector2BaseVector(f.id(), &pEigenWrapped, &accessor);
        return p;
    });

    m.def("computeKdx", [](Eigen::VectorXd pEig, Eigen::VectorXd dxEig, sofa::core::objectmodel::BaseContext *ctx) {
        sofa::core::MechanicalParams params;
        sofa::simulation::common::MechanicalOperations mop(&params, ctx);
        sofa::simulation::common::VectorOperations vop(&params, ctx);

        sofa::core::behavior::MultiVecCoord pos(&vop, sofa::core::VecCoordId::position());
        sofa::core::behavior::MultiVecDeriv vel(&vop, sofa::core::VecDerivId::velocity());
        sofa::core::behavior::MultiVecDeriv f(&vop, sofa::core::VecDerivId::force());

        sofa::core::behavior::MultiVecDeriv dx(&vop, sofa::core::VecDerivId::dx());
        dx.realloc(&vop, false, true);

        sofa::core::behavior::MultiVecDeriv df(&vop);
        dx.realloc(&vop, false, true);

        vel.clear();
        f.clear();

        // Assign the p to pos
        sofa::component::linearsolver::DefaultMultiMatrixAccessor accessor;
        mop.getMatrixDimension(nullptr, nullptr, &accessor);
        const auto n = static_cast<Eigen::Index>(accessor.getGlobalDimension());

        {
            sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> rEigenWrapped(pEig);
            mop.baseVector2MultiVector(&rEigenWrapped, pos.id(), &accessor);
        }

        mop.propagateXAndV(pos, vel);

        // Compute forces, as that can precompute stuff for the hessian
        mop.computeForce(f);

        {
            sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> rEigenWrapped(dxEig);
            mop.baseVector2MultiVector(&rEigenWrapped, dx.id(), &accessor);
        }

        mop.propagateDxAndResetDf(dx, df);
        mop.addMBKdx(df, 0.0, 0.0, 1.0);

        sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> pEigenWrapped(pEig);
        mop.multiVector2BaseVector(df.id(), &pEigenWrapped, &accessor);
        return pEig;
    });
}
