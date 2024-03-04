#include "BendingForceField.h"

#include <VNCS/DataExtensions.h>
#include <boost/range/iterator_range_core.hpp>
#include <range/v3/view/zip.hpp>

namespace
{
auto make_eigen = [](const auto &v) -> Eigen::Vector2d { return {v[0], v[1]}; };
auto view_eigen = [](const auto &v) -> Eigen::Map<const Eigen::Vector2d> {
    return Eigen::Map<const Eigen::Vector2d>(v.data());
};
}  // namespace

VNCS::Sim2D::BendingForceField::BendingForceField()
    : Inherit()
{
}

void VNCS::Sim2D::BendingForceField::init()
{
    Inherit::init();
    const auto msState = getMState();
    const auto x = msState->readRestPositions();

    for (const auto &vertex : boost::make_iterator_range(boost::vertices(m_edgeMesh))) {
        const auto &[outEdgesBegin, outEdgesEnd] = boost::out_edges(vertex, m_edgeMesh);

        for (auto edge0 = outEdgesBegin; edge0 != outEdgesEnd; ++edge0) {
            for (auto edge1 = edge0 + 1; edge1 != outEdgesEnd; ++edge1) {
                const auto leftId = boost::target(*edge0, m_edgeMesh);
                const auto rightId = boost::target(*edge1, m_edgeMesh);
                const auto centerId = vertex;

                const auto xLeft = x[leftId];
                const auto xRight = x[rightId];

                const auto du =
                    (Eigen::Map<const Eigen::Vector2d>(xLeft.data()) - Eigen::Map<const Eigen::Vector2d>(xRight.data()))
                        .norm();

                m_elements.push_back(
                    {static_cast<int>(leftId), static_cast<int>(centerId), static_cast<int>(rightId), du});
            }
        }
    }
}

VNCS::Space2D::Real VNCS::Sim2D::BendingForceField::getPotentialEnergy(const sofa::core::MechanicalParams *,
                                                                       const DataVecCoord &xVec) const
{
    VNCS::Real energy = 0;

    const auto &restPositions = getMState()->readRestPositions();
    const auto &displacements = getMState()->readPositions();

    for (const auto &element : m_elements) {
        const auto &xCenter =
            view_eigen(restPositions[element.indexCenter]) + view_eigen(displacements[element.indexCenter]);
        const auto &xLeft = view_eigen(restPositions[element.indexLeft]) + view_eigen(displacements[element.indexLeft]);
        const auto &xRight =
            view_eigen(restPositions[element.indexRight]) + view_eigen(displacements[element.indexRight]);

        const VNCS::Real duLeftRight = element.restLength;

        const Eigen::Vector2d d1 = (xRight - xCenter).normalized();
        const Eigen::Vector2d d2 = (xLeft - xCenter).normalized();

        VNCS::Real proj = -d1.dot(d2);
        if (proj > 1.) {
            proj = 1.;
        }

        const VNCS::Real theta = std::acos(proj);

        const VNCS::Real bendingStiffness = m_bendingStiffness;

        energy += bendingStiffness * theta * theta / duLeftRight;
    }

    return energy;
}

void VNCS::Sim2D::BendingForceField::addForce(const sofa::core::MechanicalParams *mparams,
                                              DataVecDeriv &fVec,
                                              const DataVecCoord &xVec,
                                              const DataVecDeriv &vVec)
{
    const auto msState = getMState();
    const auto &xRestPose = msState->readRestPositions();
    const auto &x = make_read_accessor(xVec);
    auto &f = make_write_accessor(fVec);
    m_hessians.clear();
    for (const auto &element : m_elements) {
        const auto leftIndex = element.indexLeft;
        const auto rightIndex = element.indexRight;
        const auto centerIndex = element.indexCenter;
        const auto duLeftRight = element.restLength;

        const auto xRestPoseLeft = Eigen::Map<const Eigen::Vector2d>(&xRestPose[leftIndex][0]);
        const auto xRestPoseCenter = Eigen::Map<const Eigen::Vector2d>(&xRestPose[centerIndex][0]);
        const auto xRestPoseRight = Eigen::Map<const Eigen::Vector2d>(&xRestPose[rightIndex][0]);

        const Eigen::Vector2d xLeft = Eigen::Map<const Eigen::Vector2d>(&x[leftIndex][0]) + xRestPoseLeft;
        const Eigen::Vector2d xCenter = Eigen::Map<const Eigen::Vector2d>(&x[centerIndex][0]) + xRestPoseCenter;
        const Eigen::Vector2d xRight = Eigen::Map<const Eigen::Vector2d>(&x[rightIndex][0]) + xRestPoseRight;

        Eigen::Vector2d d1 = (xRight - xCenter);
        Eigen::Vector2d d2 = (xLeft - xCenter);

        const VNCS::Real l1 = d1.norm();
        const VNCS::Real l2 = d2.norm();
        d1 /= l1;
        d2 /= l2;

        const Eigen::Matrix2d P1 = Eigen::Matrix2d::Identity() - d1 * d1.transpose();
        const Eigen::Matrix2d P2 = Eigen::Matrix2d::Identity() - d2 * d2.transpose();

        const VNCS::Real proj = std::min(1.0, -d1.dot(d2));

        const VNCS::Real theta = std::acos(proj);

        if (theta != 0.0)  // otherwise error due to division by 0 due to sin(theta), but forces are 0 anyway
        {
            // stuff related with forces
            {
                const VNCS::Real common = -2. * m_bendingStiffness * theta / (duLeftRight * std::sin(theta));
                Eigen::Vector2d fx1 = common / l1 * P1 * d2;
                Eigen::Vector2d fx2 = common / l2 * P2 * d1;

                auto make_deriv = [](const Eigen::Vector2d &v) { return VNCS::Space2D::VecType::Deriv(v[0], v[1]); };

                f[centerIndex] -= make_deriv((fx1 + fx2).head<2>());
                f[leftIndex] += make_deriv((fx2).head<2>());
                f[rightIndex] += make_deriv((fx1).head<2>());
            }

            // stuff related with hessians
            {
                const VNCS::Real sinTheta = std::sin(theta);

                const Eigen::Vector2d P1d2 = P1 * d2;
                const Eigen::Vector2d P2d1 = P2 * d1;

                // lagrangian - lagrangian

                const Eigen::Matrix2d P1d2d1t = P1d2 * d1.transpose();
                const Eigen::Matrix2d P2d1d2t = P2d1 * d2.transpose();

                const VNCS::Real commonLL = 2. * m_bendingStiffness / (duLeftRight * sinTheta);

                Eigen::Matrix2d jacobiansLL[3][3];

                jacobiansLL[2][2] = commonLL / (l1 * l1) *
                                    (theta * (P1d2d1t + P1d2d1t.transpose() - proj * P1) -
                                     1. / sinTheta * (1. - theta / sinTheta * proj) * P1d2 * P1d2.transpose());
                jacobiansLL[0][0] = commonLL / (l2 * l2) *
                                    (theta * (P2d1d2t + P2d1d2t.transpose() - proj * P2) -
                                     1. / sinTheta * (1. - theta / sinTheta * proj) * P2d1 * P2d1.transpose());
                jacobiansLL[2][0] = -commonLL / (l1 * l2) *
                                    (1. / sinTheta * (1. - theta / sinTheta * proj) * P1d2d1t + theta * P1) * P2;
                jacobiansLL[0][2] = jacobiansLL[2][0].transpose();
                jacobiansLL[2][1] = -(jacobiansLL[2][2] + jacobiansLL[2][0]);
                jacobiansLL[0][1] = -(jacobiansLL[0][2] + jacobiansLL[0][0]);
                jacobiansLL[1][2] = -(jacobiansLL[2][2] + jacobiansLL[0][2]);
                jacobiansLL[1][0] = -(jacobiansLL[2][0] + jacobiansLL[0][0]);
                jacobiansLL[1][1] = -(jacobiansLL[2][1] + jacobiansLL[0][1]);

                Eigen::Matrix<VNCS::Real, 6, 6> hessian;

                // Diagonal terms
                hessian.block<2, 2>(0, 0) = jacobiansLL[0][0];
                hessian.block<2, 2>(2, 2) = jacobiansLL[1][1];
                hessian.block<2, 2>(4, 4) = jacobiansLL[2][2];

                // Upper off diagonal
                hessian.block<2, 2>(0, 2) = jacobiansLL[0][1];
                hessian.block<2, 2>(0, 4) = jacobiansLL[0][2];
                hessian.block<2, 2>(2, 4) = jacobiansLL[1][2];

                // Lower off diagonal
                hessian.block<2, 2>(2, 0) = jacobiansLL[1][0];
                hessian.block<2, 2>(4, 0) = jacobiansLL[2][0];
                hessian.block<2, 2>(4, 2) = jacobiansLL[2][1];

                m_hessians.push_back(hessian);
            }
        } else {
            Eigen::Matrix<VNCS::Real, 6, 6> hessian;
            hessian.setZero();
            m_hessians.push_back(hessian);
        }
    }
}

void VNCS::Sim2D::BendingForceField::addDForce(const sofa::core::MechanicalParams *mparams,
                                               DataVecDeriv &fDiffVec,
                                               const DataVecDeriv &xDiffVec)
{
    const auto &dxis = make_read_accessor(xDiffVec);
    auto &dfs = make_write_accessor(fDiffVec);
    for (const auto &[element, hessian] : ranges::views::zip(m_elements, m_hessians)) {
        const auto leftIndex = element.indexLeft;
        const auto rightIndex = element.indexRight;
        const auto centerIndex = element.indexCenter;

        const auto xLeft = Eigen::Map<const Eigen::Vector2d>(&dxis[leftIndex][0]);
        const auto xCenter = Eigen::Map<const Eigen::Vector2d>(&dxis[centerIndex][0]);
        const auto xRight = Eigen::Map<const Eigen::Vector2d>(&dxis[rightIndex][0]);

        Eigen::Matrix<VNCS::Real, 6, 1> dx;
        dx << xLeft, xCenter, xRight;

        Eigen::Matrix<VNCS::Real, 6, 1> df = mparams->kFactor() * hessian * dx;

        dfs[leftIndex][0] += df[0];
        dfs[leftIndex][1] += df[1];

        dfs[centerIndex][0] += df[2];
        dfs[centerIndex][1] += df[3];

        dfs[rightIndex][0] += df[4];
        dfs[rightIndex][1] += df[5];
    }
}

namespace VNCS
{
namespace Sim2D
{
Real BendingForceField::bendingStiffness() const
{
    return m_bendingStiffness;
}

void BendingForceField::setBendingStiffness(Real newBendingStiffness)
{
    m_bendingStiffness = newBendingStiffness;
}

const std::vector<Eigen::Matrix<Real, 6, 6>> &BendingForceField::hessians() const
{
    return m_hessians;
}

void BendingForceField::setEdgeMeshPath(const std::filesystem::path &meshPath)
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
