#include <boost/range/iterator_range_core.hpp>
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <type_traits>
#include <VNCS/Spaces.h>
#include <VNCS/EdgeMesh.h>

#include <boost/range/iterator_range.hpp>

TEST_CASE("EdgeMesh")
{
    SECTION("Sim3D")
    {
        using Mesh3D = VNCS::EdgeMesh<VNCS::Space3D>;

        SECTION("Traits")
        {
            CHECK(std::is_default_constructible<Mesh3D>::value);
            CHECK(std::is_copy_constructible<Mesh3D>::value);
            CHECK(std::is_move_constructible<Mesh3D>::value);
            CHECK(std::is_copy_assignable<Mesh3D>::value);
            CHECK(std::is_move_assignable<Mesh3D>::value);
        }

        SECTION("HasPointMap")
        {
            Mesh3D mesh;
            const auto &pMap = boost::get(VNCS::PointVertexTag{}, mesh);

            const auto v0Handle = boost::add_vertex(mesh);
            boost::put(pMap, v0Handle, VNCS::Space3D::Point(1, 2, 3));

            const auto p0 = boost::get(pMap, v0Handle);
            CHECK(p0 == VNCS::Space3D::Point(1, 2, 3));
        }

        SECTION("EdgeTraversal")
        {
            Mesh3D mesh;
            boost::add_edge(0, 2, mesh);
            boost::add_edge(2, 3, mesh);

            CHECK(boost::num_vertices(mesh) == 4);
            CHECK(boost::num_edges(mesh) == 2);

            // Print (so I understand the API basically)
            const auto edgesBegin = boost::edges(mesh).first;
            const auto edgesEnd = boost::edges(mesh).second;
            for (const auto &edge : boost::make_iterator_range(edgesBegin, edgesEnd))
                std::cout << "(" << boost::source(edge, mesh) << ", " << boost::target(edge, mesh) << ")\n";
            std::cout << std::endl;

            const auto verticesBegin = boost::vertices(mesh).first;
            const auto verticesEnd = boost::vertices(mesh).second;
            for (const auto &vertex : boost::make_iterator_range(verticesBegin, verticesEnd)) {
                std::cout << "Vertex: " << vertex << "\n";
                const auto vOutEdgesBegin = boost::out_edges(vertex, mesh).first;
                const auto vOutEdgesEnd = boost::out_edges(vertex, mesh).second;

                for (const auto &edge : boost::make_iterator_range(vOutEdgesBegin, vOutEdgesEnd)) {
                    std::cout << edge << "\n";
                    std::cout << "(" << boost::source(edge, mesh) << ", " << boost::target(edge, mesh) << ")\n";
                }

                const auto vInEdgesBegin = boost::in_edges(vertex, mesh).first;
                const auto vInEdgesEnd = boost::in_edges(vertex, mesh).second;
                for (const auto &edge : boost::make_iterator_range(vInEdgesBegin, vInEdgesEnd)) {
                    std::cout << edge << "\n";
                    std::cout << "(" << boost::source(edge, mesh) << ", " << boost::target(edge, mesh) << ")\n";
                }
            }
        }
    }
}
