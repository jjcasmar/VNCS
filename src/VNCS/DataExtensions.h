#ifndef VNCS_DATA_EXTENSIONS_H
#define VNCS_DATA_EXTENSIONS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <istream>
#include <sofa/core/objectmodel/Data.h>

namespace VNCS
{
template <typename... Ts>
struct overload : Ts... {
    using Ts::operator()...;
};

template <typename... Ts>
overload(Ts...) -> overload<Ts...>;

template <typename T>
T &make_write_accessor(T &t)
{
    return sofa::helper::WriteAccessor<T>(t).wref();
}

template <typename T>
const T &make_read_accessor(const T &t)
{
    return sofa::helper::ReadAccessor<T>(t).ref();
}

template <typename T>
T &make_write_accessor(sofa::Data<T> &t)
{
    return sofa::helper::WriteAccessor<sofa::Data<T>>(t).wref();
}

template <typename T>
const T &make_read_accessor(const sofa::Data<T> &t)
{
    return sofa::helper::ReadAccessor<sofa::Data<T>>(t).ref();
}

}  // namespace VNCS

namespace Eigen
{
template <typename Scalar,
          int RowsAtCompileTime,
          int ColsAtCompileTime,
          int Options,
          int MaxRowsAtCompileTime,
          int MaxColsAtCompileTime>
std::istream &operator>>(
    std::istream &is,
    Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime> &)
{
    return is;
}

template <typename Scalar>
std::istream &operator<<(std::istream &is, const Eigen::SparseMatrix<Scalar> &)
{
    return is;
}

template <typename Scalar>
std::istream &operator>>(std::istream &is, Eigen::SparseMatrix<Scalar> &)
{
    return is;
}

template <typename Scalar, int RowsAtCompileTime>
std::istream &operator>>(std::istream &is, Eigen::DiagonalMatrix<Scalar, RowsAtCompileTime> &)
{
    return is;
}

template <typename Scalar, int RowsAtCompileTime>
std::ostream &operator<<(std::ostream &is, const Eigen::DiagonalMatrix<Scalar, RowsAtCompileTime> &m)
{
    is << m.diagonal();
    return is;
}

template <typename Scalar, int RowsAtCompileTime = Eigen::Dynamic>
Scalar *begin(Eigen::Matrix<Scalar, RowsAtCompileTime, 1> &m)
{
    return m.data();
}

template <typename Scalar, int RowsAtCompileTime = Eigen::Dynamic>
Scalar *end(Eigen::Matrix<Scalar, RowsAtCompileTime, 1> &m)
{
    return m.data() + m.rows();
}

template <typename Scalar, int RowsAtCompileTime = Eigen::Dynamic>
const Scalar *begin(const Eigen::Matrix<Scalar, RowsAtCompileTime, 1> &m)
{
    return m.data();
}

template <typename Scalar, int RowsAtCompileTime = Eigen::Dynamic>
const Scalar *end(const Eigen::Matrix<Scalar, RowsAtCompileTime, 1> &m)
{
    return m.data() + m.rows();
}

}  // namespace Eigen

#endif  // VNCS_DATA_EXTENSIONS_H
