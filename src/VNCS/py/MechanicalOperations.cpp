#include "MechanicalOperations.h"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>
#include <SofaEigen2Solver/EigenVectorWrapper.h>
#include <Eigen/Core>
#include <VNCS/Types.h>

namespace py = pybind11;

struct MechanicalOperations {
    MechanicalOperations(sofa::core::objectmodel::BaseContext *ctx)
        : m_mop(&m_params, ctx)
        , m_vop(&m_params, ctx)
        , xSofa(std::make_unique<sofa::core::behavior::MultiVecCoord>(&m_vop))
        , vSofa(std::make_unique<sofa::core::behavior::MultiVecDeriv>(&m_vop))
        , fSofa(std::make_unique<sofa::core::behavior::MultiVecDeriv>(&m_vop))
        , dxSofa(std::make_unique<sofa::core::behavior::MultiVecDeriv>(&m_vop))
        , MdxSofa(std::make_unique<sofa::core::behavior::MultiVecDeriv>(&m_vop))
        , KdxSofa(std::make_unique<sofa::core::behavior::MultiVecDeriv>(&m_vop))
    {
        xSofa->realloc(&m_vop, true, true);
        vSofa->realloc(&m_vop, true, true);
        fSofa->realloc(&m_vop, true, true);
        dxSofa->realloc(&m_vop, true, true);
        MdxSofa->realloc(&m_vop, true, true);
        KdxSofa->realloc(&m_vop, true, true);

        m_mop->setX(xSofa->id());
        m_mop->setV(vSofa->id());
        m_mop->setDx(dxSofa->id());
        m_mop->setF(fSofa->id());
    }

    int dimensions()
    {
        sofa::component::linearsolver::DefaultMultiMatrixAccessor accessor;

        m_mop.getMatrixDimension(nullptr, nullptr, &accessor);
        const auto n = static_cast<Eigen::Index>(accessor.getGlobalDimension());
        return n;
    }

    void setX(Eigen::VectorXd x)
    {
        sofa::component::linearsolver::DefaultMultiMatrixAccessor accessor;

        m_mop.getMatrixDimension(nullptr, nullptr, &accessor);
        const auto n = static_cast<Eigen::Index>(accessor.getGlobalDimension());

        if (x.size() != n)
            return;

        sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> xWrapped(x);
        m_mop.baseVector2MultiVector(&xWrapped, xSofa->id(), &accessor);

        m_mop.propagateX(xSofa->id());
    }

    void setV(Eigen::VectorXd v)
    {
        sofa::component::linearsolver::DefaultMultiMatrixAccessor accessor;

        m_mop.getMatrixDimension(nullptr, nullptr, &accessor);
        const auto n = static_cast<Eigen::Index>(accessor.getGlobalDimension());

        if (v.size() != n)
            return;

        sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> vWrapped(v);
        m_mop.baseVector2MultiVector(&vWrapped, vSofa->id(), &accessor);

        m_mop.propagateV(vSofa->id());
    }

    void propagateDx(Eigen::VectorXd dx)
    {
        sofa::component::linearsolver::DefaultMultiMatrixAccessor accessor;

        m_mop.getMatrixDimension(nullptr, nullptr, &accessor);
        const auto n = static_cast<Eigen::Index>(accessor.getGlobalDimension());

        if (dx.size() != n)
            return;

        sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> dxWrapped(dx);
        m_mop.baseVector2MultiVector(&dxWrapped, dxSofa->id(), &accessor);
    }

    std::pair<VNCS::Real, VNCS::Real> computeEnergy()
    {
        // E = 1/2 * v^t * M * v + U
        VNCS::Real kinetic = 0;
        VNCS::Real potential = 0;
        m_mop.computeEnergy(kinetic, potential);
        return std::make_pair(kinetic, potential);
    }

    Eigen::VectorXd computeForce()
    {
        m_mop.computeForce(fSofa->id());

        sofa::component::linearsolver::DefaultMultiMatrixAccessor accessor;
        m_mop.getMatrixDimension(nullptr, nullptr, &accessor);
        const auto n = static_cast<Eigen::Index>(accessor.getGlobalDimension());
        Eigen::VectorXd f;
        f.resize(n);

        sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> fWrapped(f);
        m_mop.multiVector2BaseVector(fSofa->id(), &fWrapped, &accessor);

        return f;
    }

    Eigen::VectorXd Mdx()
    {
        sofa::component::linearsolver::DefaultMultiMatrixAccessor accessor;
        m_mop.addMBKdx(MdxSofa->id(), 1.0, 0.0, 0.0);

        m_mop.getMatrixDimension(nullptr, nullptr, &accessor);
        const auto n = static_cast<Eigen::Index>(accessor.getGlobalDimension());
        Eigen::VectorXd Mdx;
        Mdx.resize(n);

        sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> MdxWrapped(Mdx);
        m_mop.multiVector2BaseVector(MdxSofa->id(), &MdxWrapped, &accessor);

        return Mdx;
    }

    Eigen::VectorXd Mdx(Eigen::VectorXd dx)
    {
        propagateDx(dx);
        return Mdx();
    }

    Eigen::VectorXd Kdx()
    {
        sofa::component::linearsolver::DefaultMultiMatrixAccessor accessor;
        m_mop.addMBKdx(KdxSofa->id(), 0.0, 0.0, 1.0);

        m_mop.getMatrixDimension(nullptr, nullptr, &accessor);
        const auto n = static_cast<Eigen::Index>(accessor.getGlobalDimension());
        Eigen::VectorXd Kdx;
        Kdx.resize(n);

        sofa::component::linearsolver::EigenVectorWrapper<VNCS::Real> KdxWrapped(Kdx);
        m_mop.multiVector2BaseVector(KdxSofa->id(), &KdxWrapped, &accessor);

        return Kdx;
    }

    Eigen::VectorXd Kdx(Eigen::VectorXd dx)
    {
        propagateDx(dx);
        return Kdx();
    }

private:
    sofa::core::MechanicalParams m_params;
    sofa::simulation::common::MechanicalOperations m_mop;
    sofa::simulation::common::VectorOperations m_vop;
    std::unique_ptr<sofa::core::behavior::MultiVecCoord> xSofa;
    std::unique_ptr<sofa::core::behavior::MultiVecDeriv> vSofa;
    std::unique_ptr<sofa::core::behavior::MultiVecDeriv> fSofa;
    std::unique_ptr<sofa::core::behavior::MultiVecDeriv> dxSofa;
    std::unique_ptr<sofa::core::behavior::MultiVecDeriv> MdxSofa;
    std::unique_ptr<sofa::core::behavior::MultiVecDeriv> KdxSofa;
};

PYBIND11_DECLARE_HOLDER_TYPE(T, sofa::core::sptr<T>, true);
void VNCS::py::mechanicalOperations(pybind11::module &m)
{
    ::py::class_<MechanicalOperations>(m, "MechanicalOperations")
        .def(::py::init([](sofa::core::objectmodel::BaseContext *ctx) { return MechanicalOperations(ctx); }))
        .def_property_readonly("dimensions", &MechanicalOperations::dimensions)
        .def("setX", &MechanicalOperations::setX)
        .def("setV", &MechanicalOperations::setV)
        .def_property_readonly("energy", &MechanicalOperations::computeEnergy)
        .def("propagateDx", &MechanicalOperations::propagateDx)
        .def("computeForce", &MechanicalOperations::computeForce)
        .def("setX", &MechanicalOperations::setX)
        .def("Mdx", ::py::overload_cast<Eigen::VectorXd>(&MechanicalOperations::Mdx))
        .def("Mdx", ::py::overload_cast<>(&MechanicalOperations::Mdx))
        .def("Kdx", ::py::overload_cast<Eigen::VectorXd>(&MechanicalOperations::Kdx))
        .def("Kdx", ::py::overload_cast<>(&MechanicalOperations::Kdx));
}
