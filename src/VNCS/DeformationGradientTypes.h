#ifndef VNCS_DEFORMATIONGRADIENTTYPES_H
#define VNCS_DEFORMATIONGRADIENTTYPES_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/MapMapSparseMatrix.h>

#include <sofa/defaulttype/DataTypeInfo.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <SofaBaseMechanics/MechanicalObject.inl>

namespace VNCS
{
template <typename R, typename WorldSpace_, typename MaterialSpace_>
class F
{
public:
    using Real = R;
    using VecReal = sofa::helper::vector<Real>;
    using WorldSpace = WorldSpace_;
    using MaterialSpace = MaterialSpace_;
    using Matrix = Eigen::Matrix<Real, WorldSpace::dim, MaterialSpace::dim>;

    struct Coord : public Matrix {
        enum { total_size = WorldSpace::dim * MaterialSpace::dim };

        Coord() { Matrix::setZero(); }

        template <typename OtherDerived>
        Coord(const Eigen::MatrixBase<OtherDerived> &other)
            : Matrix(other)
        {
        }

        template <typename OtherDerived>
        Coord &operator=(const Eigen::MatrixBase<OtherDerived> &other)
        {
            this->Matrix::operator=(other);
            return *this;
        }

        void clear() { Matrix::setZero(); }
        constexpr static int size() { return total_size; }
    };

    enum { coord_total_size = MaterialSpace::dim };
    enum { deriv_total_size = MaterialSpace::dim };
    enum { spatial_dimensions = WorldSpace::dim };

    using Deriv = Coord;
    using VecCoord = sofa::helper::vector<Coord>;
    using VecDeriv = sofa::helper::vector<Deriv>;
    using MatrixDeriv = sofa::defaulttype::MapMapSparseMatrix<Deriv>;

    static void get(Real &x, Real &y, Real &z, const Deriv &v) {}
    static void set(Coord &v, Real x, Real y, Real z) {}
    static void add(Coord &v, Real x, Real y, Real z) {}
    static Deriv interpolate(const sofa::helper::vector<Deriv> &ancestors, const sofa::helper::vector<Real> &coefs)
    {
        Deriv c = Deriv::Zero();
        for (int i = 0; i < ancestors.size(); ++i)
            c += ancestors[i] * coefs[i];
        return c;
    }
};

}  // namespace VNCS

template <typename R, typename WorldSpace, typename MaterialSpace>
std::istream &operator>>(std::istream &s, typename VNCS::F<R, WorldSpace, MaterialSpace>::Coord &d)
{
    return s;
}

#endif  // DEFORMATIONGRADIENTTYPES_H
