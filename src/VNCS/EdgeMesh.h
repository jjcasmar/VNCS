#ifndef VNCS_EDGEMESH_H
#define VNCS_EDGEMESH_H

#include <boost/graph/properties.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <vector>

namespace VNCS
{
struct PointVertexTag {
    using kind = boost::vertex_property_tag;
};

template <typename Space>
using PointVertexProperty = boost::property<PointVertexTag, typename Space::Point>;

template <typename Space>
using EdgeMesh = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, PointVertexProperty<Space>>;

template <typename Space>
EdgeMesh<Space> lines_soup_to_edge_mesh(const std::vector<typename Space::Point> &points,
                                        const std::vector<std::vector<std::size_t>> &lines)
{
    EdgeMesh<Space> mesh;
    const auto pMap = boost::get(PointVertexTag{}, mesh);
    for (const auto &p : points) {
        const auto vertex = boost::add_vertex(mesh);
        boost::put(pMap, vertex, p);
    }

    for (const auto &line : lines) {
        for (int i = 1; i < line.size(); ++i)
            boost::add_edge(line[i - 1], line[i], mesh);
    }

    return mesh;
}

template <typename P>
bool read_OBJ(std::istream &input, std::vector<P> &points, std::vector<std::vector<std::size_t>> &segments)
{
    P p;
    std::string line;
    while (getline(input, line)) {
        if (line[0] == 'v' && line[1] == ' ') {
            std::istringstream iss(line.substr(1));
            iss >> p;
            if (!iss)
                return false;
            points.push_back(p);
        } else if (line[0] == 'l') {
            std::istringstream iss(line.substr(1));
            int i;
            segments.push_back(std::vector<std::size_t>());
            while (iss >> i) {
                if (i < 1) {
                    segments.back().push_back(points.size() + i);  // negative indices are relative references
                } else {
                    segments.back().push_back(i - 1);
                }
                iss.ignore(256, ' ');
            }
        } else {
            // std::cerr<<"ERROR : Cannnot read line beginning with "<<line[0]<<std::endl;
            continue;
        }
    }
    return true;
}
}  // namespace VNCS

#endif  // VNCS_EDGEMESH_H