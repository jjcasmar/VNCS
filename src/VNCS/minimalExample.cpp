#include <PhySim/Physics/DomainDistribution.h>
#include <PhySim/Physics/DoFSet.h>
#include <PhySim/Physics/Simulables/Simulable_FEM_Volume.h>
#include <PhySim/Geometry/Meshes/Mesh_Tetra.h>
#include <PhySim/Geometry/Polytopes/ShapeFunctions.h>
#include <PhySim/Geometry/Polytopes/Quadratures.h>

using namespace PhySim;

int main()
{
    // Create model

    PtrS<Simulable_FEM_Volume> pModel(new Simulable_FEM_Volume());

    // Set mesh

    MatrixXd mV;  // N x 3 matrix with vertex coordinates
    MatrixXi mT;  // M x 4 matrix with tetra connectivity
    PtrS<Mesh_Tetra> pMesh(new Mesh_Tetra(mV, mT, Discretization::Discretization_Tet4));
    pMesh->SetNodesTrait(mV, PhySim::Tag_Position_0);    // Rest = initial config
    pMesh->SetNodesTrait(0 * mV, PhySim::Tag_Velocity);  // Zero initial velocity
    pModel->SetupOptions().m_pMesh = pMesh;

    // Set shape function, quadrature and material type
    // These are constant throughout the whole domain

    // StVK material for 3D polytopes in 3D space
    pModel->SetupOptions().m_pMatModels.reset(
        new DomainDistribution_Constant<PtrS<IMaterialModel>>(MaterialModel_3Din3D_StVK::Instance()));

    // Linear shape function for 4-node tetrahedron
    pModel->SetupOptions().m_pShapeFunctions.reset(
        new DomainDistribution_Constant<PtrS<IShapeFunction>>(ShapeFunction_Tet4::Instance()));

    // Just one quadrature point within a tetrahedron
    pModel->SetupOptions().m_pQuadratures.reset(
        new DomainDistribution_Constant<PtrS<IQuadrature>>(Quadrature_Tet1::Instance()));

    // Set material properties (heterogeneous distribution), in this
    // example, a different one per tetrahedron but it can be set to
    // an indexed library of materials.

    vector<PtrS<ParameterSet>> vmatPars(mT.rows());
    for (int i = 0; i < mT.rows(); ++i) {
        vmatPars[i] = std::make_shared<ParameterSet>();
        vmatPars[i]->InitFromYoungPoisson(1e6, 0.3, 1000);
        vmatPars[i]->AddParameter("ALPHA", 0.42);
    }

    pModel->SetupOptions().m_pMatParams.reset(new DomainDistribution_Variable<PtrS<ParameterSet>>(vmatPars));

    // Setup (initialize) model

    pModel->Setup();

    // Some things you can do

    VectorXd vx;
    VectorXd v0;
    VectorXd vv;
    pModel->GetDOFVector(vx, PhySim::Tag_Position_X);
    pModel->GetDOFVector(v0, PhySim::Tag_Position_0);
    pModel->GetDOFVector(vv, PhySim::Tag_Velocity);

    pModel->SetDOFVector(vx, PhySim::Tag_Position_X);  // Kinematics and mechanics now dirty

    double energy;
    AVectorXd vg;  // Assembled vector (extends Eigen)
    AMatrixSd mH;  // Assembled Matrix (extends Eigen)
    pModel->GetEnergy(energy);
    pModel->GetGradient(vg);
    pModel->GetHessian(mH);

    return 0;
}
