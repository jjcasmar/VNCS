#ifndef VNCS_GENERATORS_ADAPTIVECRITERIA_H
#define VNCS_GENERATORS_ADAPTIVECRITERIA_H

#include <CGAL/Mesh_2/Face_badness.h>
#include <CGAL/Delaunay_mesh_criteria_2.h>
#include <utility>
#include <ostream>
#include <VNCS/Spaces.h>

namespace VNCS
{
namespace Generators
{
class AdaptiveMeshCriteriaFunctor
{
public:
    virtual ~AdaptiveMeshCriteriaFunctor() = default;

    virtual std::pair<VNCS::Space2D::Real, VNCS::Space2D::Real> operator()(const VNCS::Space2D::Point &p) const = 0;
};

template <class CDT>
class AdaptiveMeshCriteria
{
protected:
    typedef typename CDT::Geom_traits Geom_traits;
    Geom_traits m_geomTraits;
    std::shared_ptr<AdaptiveMeshCriteriaFunctor> m_sizeFunctor;

public:
    AdaptiveMeshCriteria(const Geom_traits &traits = Geom_traits())
        : m_geomTraits(traits)
    {
    }

    AdaptiveMeshCriteria(std::shared_ptr<AdaptiveMeshCriteriaFunctor> sizeFunctor,
                         const Geom_traits &traits = Geom_traits())
        : m_geomTraits(traits)
        , m_sizeFunctor(sizeFunctor)
    {
    }

    // first: squared_minimum_sine
    // second: size
    struct Quality : public std::pair<double, double> {
        typedef std::pair<double, double> Base;
        double B;

        Quality()
            : Base(){};
        Quality(double _sine, double _size)
            : Base(_sine, _size)
        {
        }

        const double &size() const { return second; }
        const double &sine() const { return first; }

        // q1<q2 means q1 is prioritised over q2
        // ( q1 == *this, q2 == q )
        bool operator<(const Quality &q) const
        {
            if (size() > 1)
                if (q.size() > 1)
                    return (size() > q.size());
                else
                    return true;  // *this is big but not q
            else if (q.size() > 1)
                return false;  // q is big but not *this
            return (sine() < q.sine());
        }

        std::ostream &operator<<(std::ostream &out) const
        {
            return out << "(size=" << size() << ", sine=" << sine() << ")";
        }
    };

    class Is_bad
    {
    protected:
        Geom_traits traits;
        std::shared_ptr<AdaptiveMeshCriteriaFunctor> sizeFunctor;

    public:
        Is_bad(std::shared_ptr<AdaptiveMeshCriteriaFunctor> sizeFunctor, Geom_traits traits = Geom_traits())
            : traits(traits)
            , sizeFunctor(sizeFunctor)
        {
        }

        CGAL::Mesh_2::Face_badness operator()(const Quality q) const
        {
            if (q.size() > 1)
                return CGAL::Mesh_2::IMPERATIVELY_BAD;
            if (q.sine() < q.B)
                return CGAL::Mesh_2::BAD;
            else
                return CGAL::Mesh_2::NOT_BAD;
        }

        CGAL::Mesh_2::Face_badness operator()(const typename CDT::Face_handle &fh, Quality &q) const
        {
            typedef typename CDT::Geom_traits Geom_traits;
            using Point_2 = typename Geom_traits::Point_2;
            typedef typename Geom_traits::Compute_area_2 Compute_area_2;
            typedef typename Geom_traits::Compute_squared_distance_2 Compute_squared_distance_2;

            Geom_traits traits; /** @warning traits with data!! */

            Compute_squared_distance_2 squared_distance = traits.compute_squared_distance_2_object();

            const Point_2 &pa = fh->vertex(0)->point();
            const Point_2 &pb = fh->vertex(1)->point();
            const Point_2 &pc = fh->vertex(2)->point();

            const auto blend_a = (*sizeFunctor)(pa);
            const auto blend_b = (*sizeFunctor)(pb);
            const auto blend_c = (*sizeFunctor)(pc);

            const auto sizeBound = (blend_a.first + blend_b.first + blend_c.first) / 3.0;
            const auto squared_size_bound = sizeBound * sizeBound;
            const auto B = (blend_a.second + blend_b.second + blend_c.second) / 3.0;
            q.B = B;

            double a = CGAL::to_double(squared_distance(pb, pc));
            double b = CGAL::to_double(squared_distance(pc, pa));
            double c = CGAL::to_double(squared_distance(pa, pb));

            double max_sq_length;  // squared max edge length
            double second_max_sq_length;

            if (a < b) {
                if (b < c) {
                    max_sq_length = c;
                    second_max_sq_length = b;
                } else {  // c<=b
                    max_sq_length = b;
                    second_max_sq_length = (a < c ? c : a);
                }
            } else  // b<=a
            {
                if (a < c) {
                    max_sq_length = c;
                    second_max_sq_length = a;
                } else {  // c<=a
                    max_sq_length = a;
                    second_max_sq_length = (b < c ? c : b);
                }
            }

            q.second = 0;
            if (squared_size_bound != 0) {
                //	  std::cerr << squared_size_bound << std::endl;
                q.second = max_sq_length / squared_size_bound;
                // normalized by size bound to deal
                // with size field
                if (q.size() > 1) {
                    q.first = 1;  // (do not compute sine)
                    return CGAL::Mesh_2::IMPERATIVELY_BAD;
                }
            }

            Compute_area_2 area_2 = traits.compute_area_2_object();

            double area = 2 * CGAL::to_double(area_2(pa, pb, pc));

            q.first = (area * area) / (max_sq_length * second_max_sq_length);  // (sine)

            if (q.sine() < B)
                return CGAL::Mesh_2::BAD;
            else
                return CGAL::Mesh_2::NOT_BAD;
        }
    };

    Is_bad is_bad_object() const { return Is_bad(m_sizeFunctor); }
};

}  // namespace Generators
}  // namespace VNCS

#endif  //  VNCS_GENERATORS_ADAPTIVECRITERIA_H
