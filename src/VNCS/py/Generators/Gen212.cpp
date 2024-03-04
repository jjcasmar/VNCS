#include "Gen212.h"

#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <VNCS/Generators/Gen212/Gen212.h>
#include "AdaptiveMeshCriteriaFunctor.h"
#include "VNCS/EdgeMesh.h"
#include <Eigen/Dense>

#include <boost/range/iterator_range_core.hpp>
#include <string>
#include <filesystem>

namespace py = pybind11;

void VNCS::Generators::py::gen212(::py::module &m)
{
    class EdgeMeshBind
    {
    public:
        EdgeMeshBind(const std::string &path)
        {
            m_edgeMesh = [path]() {
                std::vector<VNCS::Space2D::Point> fineMeshPoints;
                std::vector<std::vector<std::size_t>> coarseLines;
                std::ifstream fineMeshIn(path);
                VNCS::read_OBJ(fineMeshIn, fineMeshPoints, coarseLines);

                return VNCS::lines_soup_to_edge_mesh<VNCS::Space2D>(fineMeshPoints, coarseLines);
            }();
        }

        const std::vector<Eigen::Vector2d> vertices() const
        {
            std::vector<Eigen::Vector2d> r;
            auto pMap = boost::get(VNCS::PointVertexTag{}, m_edgeMesh);
            for (const auto vertex : boost::make_iterator_range(boost::vertices(m_edgeMesh))) {
                const auto &p = pMap[vertex];
                r.push_back({p[0], p[1]});
            }
            return r;
        }

        const std::vector<std::pair<int, int>> edges() const
        {
            std::vector<std::pair<int, int>> r;
            for (const auto edge : boost::make_iterator_range(boost::edges(m_edgeMesh))) {
                r.push_back({boost::source(edge, m_edgeMesh), boost::target(edge, m_edgeMesh)});
            }
            return r;
        }

        VNCS::EdgeMesh<VNCS::Space2D> m_edgeMesh;
    };
    ::py::class_<EdgeMeshBind>(m, "EdgeMesh")
        .def(pybind11::init<std::string>())
        .def_property_readonly("vertices", &EdgeMeshBind::vertices)
        .def_property_readonly("edges", &EdgeMeshBind::edges);

    ::py::class_<VNCS::Generators::Gen212::Generator>(m, "Gen212")
        .def(pybind11::init<>())
        .def_property(
            "coarseInputMeshPath",
            [](const VNCS::Generators::Gen212::Generator &mo) {},
            [](VNCS::Generators::Gen212::Generator &gen212, const std::string &meshPath) {
                gen212.setCoarseMeshPath(meshPath);
            })
        .def_property(
            "fineInputMeshPath",
            [](const VNCS::Generators::Gen212::Generator &mo) {},
            [](VNCS::Generators::Gen212::Generator &gen212, const std::string &meshPath) {
                gen212.setEdgeMeshPath(meshPath);
            })
        .def_property("remeshCoarseInputMesh",
                      &VNCS::Generators::Gen212::Generator::remeshCoarseMesh,
                      &VNCS::Generators::Gen212::Generator::setRemeshCoarseMesh)
        .def_property(
            "coarseCriteria",
            [](const VNCS::Generators::Gen212::Generator &mo) { return mo.coarseCriteria(); },
            [](VNCS::Generators::Gen212::Generator &gen212,
               std::shared_ptr<VNCS::Generators::py::AdaptiveMeshCriteriaFunctor> criteria) {
                gen212.setCoarseCriteria(criteria);
            })
        .def_property(
            "coarseSamplersFilePath",
            [](const VNCS::Generators::Gen212::Generator &mo) {},
            [](VNCS::Generators::Gen212::Generator &gen212, std::string &coarseSamplersFilePath) {
                gen212.setCoarseSamplersFilePath(coarseSamplersFilePath);
            })
        .def_property(
            "fineSamplersFilePath",
            [](const VNCS::Generators::Gen212::Generator &mo) {},
            [](VNCS::Generators::Gen212::Generator &gen212, std::string &fineSamplersFilePath) {
                gen212.setFineSamplersFilePath(fineSamplersFilePath);
            })
        .def_property(
            "clusterMatrixFilePath",
            [](const VNCS::Generators::Gen212::Generator &mo) {},
            [](VNCS::Generators::Gen212::Generator &gen212, std::string &clusterMatrixFilePath) {
                gen212.setClusterMatrixFilePath(clusterMatrixFilePath);
            })
        .def_property(
            "dofFilePath",
            [](const VNCS::Generators::Gen212::Generator &mo) {},
            [](VNCS::Generators::Gen212::Generator &gen212, std::string &dofFilePath) {
                gen212.setDofFilePath(dofFilePath);
            })
        .def("create", &VNCS::Generators::Gen212::Generator::operator());
}
