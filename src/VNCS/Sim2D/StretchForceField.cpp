#include "StretchForceField.h"

#include <VNCS/DataExtensions.h>
#include <boost/range/iterator_range_core.hpp>
#include <range/v3/view/zip.hpp>
#include <range/v3/numeric/accumulate.hpp>

#include <VNCS/Logger.h>

namespace
{
VNCS::Sim2D::StretchForceField::Deriv make_deriv(const Eigen::Vector2d &v)
{
    return {v[0], v[1]};
}
}  // namespace

VNCS::Sim2D::StretchForceField::StretchForceField()
    : Inherit()
{
}

void VNCS::Sim2D::StretchForceField::init()
{
    Inherit::init();
    const auto msState = getMState();
    const auto x = msState->readRestPositions();

    for (const auto &edge : boost::make_iterator_range(boost::edges(m_edgeMesh))) {
        const auto source = boost::source(edge, m_edgeMesh);
        const auto target = boost::target(edge, m_edgeMesh);

        const auto x0 = x[source];
        const auto x1 = x[target];

        const Eigen::Vector2d L0 =
            (Eigen::Map<const Eigen::Vector2d>(x1.data()) - Eigen::Map<const Eigen::Vector2d>(x0.data()));

        m_elements.push_back({static_cast<int>(source), static_cast<int>(target), L0});
    }
}

VNCS::Space3D::Real VNCS::Sim2D::StretchForceField::getPotentialEnergy(const sofa::core::MechanicalParams *,
                                                                       const DataVecCoord &xVec) const
{
    const auto &u = make_read_accessor(xVec);

    auto energyFunctor = [this, &u](const VNCS::Sim2D::StretchForceField::StretchElement &element) -> VNCS::Real {
        const auto leftIndex = element.i0;
        const auto rightIndex = element.i1;

        const Eigen::Vector2d u0 = Eigen::Map<const Eigen::Vector2d>(&u[leftIndex][0]);
        const Eigen::Vector2d u1 = Eigen::Map<const Eigen::Vector2d>(&u[rightIndex][0]);

        const VNCS::Real L0norm = element.restLength.norm();
        const Eigen::Vector2d L = element.restLength + u1 - u0;
        const VNCS::Real Lnorm = L.norm();

        return 0.5 * m_stretchStiffness * L0norm * pow(((Lnorm / L0norm) - 1.0), 2.0);
        // }
    };

    return ranges::accumulate(m_elements, 0.0, std::plus<>(), energyFunctor);
}

void VNCS::Sim2D::StretchForceField::addForce(const sofa::core::MechanicalParams *mparams,
                                              DataVecDeriv &fVec,
                                              const DataVecCoord &xVec,
                                              const DataVecDeriv &vVec)
{
    const auto msState = getMState();
    const auto &u = make_read_accessor(xVec);
    auto &f = make_write_accessor(fVec);
    m_hessians.clear();
    for (const auto &element : m_elements) {
        const auto i0 = element.i0;
        const auto i1 = element.i1;

        const Eigen::Vector2d u0 = Eigen::Map<const Eigen::Vector2d>(&u[i0][0]);
        const Eigen::Vector2d u1 = Eigen::Map<const Eigen::Vector2d>(&u[i1][0]);

        {
            const auto A = []() {
                Eigen::Matrix<VNCS::Real, 2, 4> A;
                A.block<2, 2>(0, 0) = Eigen::Matrix2d::Identity();
                A.block<2, 2>(0, 2) = -Eigen::Matrix2d::Identity();
                return A;
            }();

            const auto u = [&u1, &u0]() {
                Eigen::Matrix<VNCS::Real, 4, 1> u;
                u.block<2, 1>(0, 0) = u1;
                u.block<2, 1>(2, 0) = u0;
                return u;
            }();

            const auto b = element.restLength;

            const auto L0 = element.restLength;
            const auto L0norm = L0.norm();

            const auto L = (A * u + b).eval();
            const auto Lnorm = L.norm();

            const auto eps = Lnorm / L0norm - 1;

            const auto depsdu = (1.0 / (L0norm * L0norm * (eps + 1)) * (A * u + b).transpose() * A).eval();

            const auto d2epsdu2 = 1.0 / (L0norm * L0norm * (eps + 1)) *
                                  (A.transpose() * A - L0norm * L0norm * depsdu.transpose() * depsdu);

            const auto dEdu_MA = (-(m_stretchStiffness * L0norm * eps * depsdu).transpose()).eval();
            const auto d2Edu = (m_stretchStiffness * L0norm * eps * d2epsdu2 +
                                m_stretchStiffness * L0norm * depsdu.transpose() * depsdu)
                                   .eval();
            f[element.i0] -= ::make_deriv(dEdu_MA.block<2, 1>(0, 0));
            f[element.i1] -= ::make_deriv(dEdu_MA.block<2, 1>(2, 0));

            m_hessians.push_back(-d2Edu);
        }
    }
}

void VNCS::Sim2D::StretchForceField::addDForce(const sofa::core::MechanicalParams *mparams,
                                               DataVecDeriv &fDiffVec,
                                               const DataVecDeriv &xDiffVec)
{
    const auto &dxis = make_read_accessor(xDiffVec);
    auto &dfs = make_write_accessor(fDiffVec);
    for (const auto &[element, hessian] : ranges::views::zip(m_elements, m_hessians)) {
        const auto leftIndex = element.i0;
        const auto rightIndex = element.i1;

        const auto xLeft = Eigen::Map<const Eigen::Vector2d>(&dxis[leftIndex][0]);
        const auto xRight = Eigen::Map<const Eigen::Vector2d>(&dxis[rightIndex][0]);

        Eigen::Matrix<VNCS::Real, 4, 1> dx;
        dx << xLeft, xRight;

        Eigen::Matrix<VNCS::Real, 4, 1> df = mparams->kFactor() * hessian * dx;

        dfs[leftIndex][0] += df[0];
        dfs[leftIndex][1] += df[1];

        dfs[rightIndex][0] += df[2];
        dfs[rightIndex][1] += df[3];
    }
}

namespace VNCS
{
namespace Sim2D
{
Real StretchForceField::stretchStiffness() const
{
    return m_stretchStiffness;
}

void StretchForceField::setStretchStiffness(Real newStretchStiffness)
{
    m_stretchStiffness = newStretchStiffness;
}

const std::vector<Eigen::Matrix<Real, 4, 4>> &StretchForceField::hessians() const
{
    return m_hessians;
}

void StretchForceField::setEdgeMeshPath(const std::filesystem::path &meshPath)
{
    m_edgeMesh = [meshPath]() {
        std::vector<VNCS::Space2D::Point> fineMeshPoints;
        std::vector<std::vector<std::size_t>> coarseLines;
        std::ifstream fineMeshIn(meshPath);
        VNCS::read_OBJ(fineMeshIn, fineMeshPoints, coarseLines);

        return VNCS::lines_soup_to_edge_mesh<VNCS::Space2D>(fineMeshPoints, coarseLines);
    }();
}
}  // namespace Sim2D
}  // namespace VNCS
