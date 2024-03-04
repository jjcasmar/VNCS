#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <VNCS/Sim3D/TriangleFabricBendingEnergyCommon.h>
#include <array>

#include <sofa/core/objectmodel/Context.h>
#include <sofa/core/VecId.h>
using sofa::core::VecId;
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/VecTypes.h>
using sofa::defaulttype::Vec3Types;
using Coord3 = sofa::defaulttype::Vector3;
using Deriv3 = sofa::defaulttype::Vector3;
using VecCoord3 = sofa::helper::vector<Coord3>;
#include <sofa/helper/system/FileRepository.h>
#include <sofa/simulation/Node.h>

#include <SofaSimulationGraph/DAGSimulation.h>
#include <SofaSimulationGraph/init.h>
#include <sofa/simulation/MechanicalOperations.h>

#include <SofaBaseMechanics/MechanicalObject.h>
#include <SofaBaseMechanics/MechanicalObject.inl>
using MechanicalObject3 = sofa::component::container::MechanicalObject<Vec3Types>;

using sofa::core::objectmodel::Data;
using sofa::core::objectmodel::New;
using sofa::helper::ReadAccessor;
using sofa::helper::WriteAccessor;
using sofa::simulation::Node;

TEST_CASE("TriangleFabricBendingEnergy")
{
    struct Parameter {
        std::array<Coord3, 4> positions;
    };
    MechanicalObject3::SPtr mechanicalObject;
    sofa::core::sptr<sofa::simulation::Node> node;
    sofa::core::sptr<VNCS::Sim3D::TriangleFabricBendingEnergy> bendingEnergy;

    // The graph root node : gravity already exists in a GNode by default
    sofa::simulation::setSimulation(new sofa::simulation::graph::DAGSimulation());
    node = sofa::simulation::getSimulation()->createNewGraph("root");
    node->setGravity(sofa::defaulttype::Vector3(0, 0, 0));

    mechanicalObject = sofa::core::objectmodel::New<MechanicalObject3>();
    WriteAccessor<Data<VecCoord3>> x = *mechanicalObject->write(VecId::restPosition());
    node->addObject(mechanicalObject);
    mechanicalObject->resize(4);
    // get write access to the position vector of mechanical object DOF
    x[0] = Coord3(1, -1, 0);
    x[1] = Coord3(1, 1, 0);
    x[2] = Coord3(-1, -1, 0);
    x[3] = Coord3(-1, 1, 0);

    // Tetrahedron force field
    bendingEnergy = sofa::core::objectmodel::New<VNCS::Sim3D::TriangleFabricBendingEnergy>();
    node->addObject(bendingEnergy);
    bendingEnergy->setBendingStiffness(1);
    bendingEnergy->setMeshFilepath(std::string(ASSETS_DIR) + std::string("/coarse.obj"));

    sofa::simulation::graph::init();
    sofa::simulation::getSimulation()->init(node.get());

    auto param = GENERATE(Parameter{{Coord3{0, 0, 0}, Coord3{0, 0, 0}, Coord3{0, 0, 0}, Coord3{1, -1, std::sqrt(2.0)}}},
                          Parameter{{Coord3{0, 0, 0}, Coord3{0, 0, 0}, Coord3{0, 0, 0}, Coord3{0, 0, 0}}});
    SECTION("Forces")
    {
        WriteAccessor<Data<VecCoord3>> x = *mechanicalObject->write(VecId::position());
        x[0] = param.positions[0];
        x[1] = param.positions[1];
        x[2] = param.positions[2];
        x[3] = param.positions[3];

        sofa::core::MechanicalParams params;
        sofa::simulation::common::MechanicalOperations mop(&params, node.get());

        bendingEnergy->addForce(nullptr,
                                *mechanicalObject->write(VecId::force()),
                                *mechanicalObject->read(VecId::position()),
                                *mechanicalObject->read(VecId::velocity()));

        ReadAccessor<Data<VecCoord3>> force = *mechanicalObject->read(VecId::force());

        std::array<Deriv3, 4> expectedForces = {Deriv3{force[0][0], force[0][1], force[0][2]},
                                                Deriv3{force[1][0], force[1][1], force[1][2]},
                                                Deriv3{force[2][0], force[2][1], force[2][2]},
                                                Deriv3{force[3][0], force[3][1], force[3][2]}};

        // Compute by finite differences
        const auto epsilon = 1e-8;
        for (int nodeId = 0; nodeId < 4; ++nodeId) {
            for (int dirId = 0; dirId < 3; ++dirId) {
                const auto originalValue = x[nodeId][dirId];
                x[nodeId][dirId] -= epsilon;
                const auto backwardEnergy =
                    bendingEnergy->getPotentialEnergy(nullptr, *mechanicalObject->read(VecId::position()));
                x[nodeId][dirId] = originalValue;

                x[nodeId][dirId] += epsilon;
                const auto forwardEnergy =
                    bendingEnergy->getPotentialEnergy(nullptr, *mechanicalObject->read(VecId::position()));
                x[nodeId][dirId] = originalValue;

                const auto diff = -(forwardEnergy - backwardEnergy) / (2.0 * epsilon);
                REQUIRE(expectedForces[nodeId][dirId] - diff == Approx(0.0).margin(1e-5));
            }
        }
    }

    SECTION("Hessian")
    {
        WriteAccessor<Data<VecCoord3>> x = *mechanicalObject->write(VecId::position());
        x[0] = param.positions[0];
        x[1] = param.positions[1];
        x[2] = param.positions[2];
        x[3] = param.positions[3];

        sofa::core::MechanicalParams params;
        sofa::simulation::common::MechanicalOperations mop(&params, node.get());

        bendingEnergy->addForce(nullptr,
                                *mechanicalObject->write(VecId::force()),
                                *mechanicalObject->read(VecId::position()),
                                *mechanicalObject->read(VecId::velocity()));

        ReadAccessor<Data<VecCoord3>> force = *mechanicalObject->read(VecId::force());
        mechanicalObject->resetForce({});

        Eigen::Matrix<VNCS::Real, 12, 12> expectedHessian = bendingEnergy->hessians()[0];
        Eigen::Matrix<VNCS::Real, 12, 12> fdHessian;

        REQUIRE((expectedHessian.transpose().eval() - expectedHessian).isZero());

        // Compute by finite differences
        const auto epsilon = 1e-8;
        for (int nodeId = 0; nodeId < 4; ++nodeId) {
            for (int dirId = 0; dirId < 3; ++dirId) {
                const auto originalValue = x[nodeId][dirId];
                x[nodeId][dirId] -= epsilon;
                bendingEnergy->addForce(nullptr,
                                        *mechanicalObject->write(VecId::force()),
                                        *mechanicalObject->read(VecId::position()),
                                        *mechanicalObject->read(VecId::velocity()));
                ReadAccessor<Data<VecCoord3>> force = *mechanicalObject->read(VecId::force());
                x[nodeId][dirId] = originalValue;

                Eigen::Matrix<VNCS::Real, 12, 1> backwardForce;
                backwardForce << force[0][0],  //
                    force[0][1],               //
                    force[0][2],               //
                    force[1][0],               //
                    force[1][1],               //
                    force[1][2],               //
                    force[2][0],               //
                    force[2][1],               //
                    force[2][2],               //
                    force[3][0],               //
                    force[3][1],               //
                    force[3][2];
                mechanicalObject->resetForce({});

                x[nodeId][dirId] += epsilon;
                bendingEnergy->addForce(nullptr,
                                        *mechanicalObject->write(VecId::force()),
                                        *mechanicalObject->read(VecId::position()),
                                        *mechanicalObject->read(VecId::velocity()));
                force = *mechanicalObject->read(VecId::force());
                x[nodeId][dirId] = originalValue;

                Eigen::Matrix<VNCS::Real, 12, 1> forwardForce;
                forwardForce << force[0][0],  //
                    force[0][1],              //
                    force[0][2],              //
                    force[1][0],              //
                    force[1][1],              //
                    force[1][2],              //
                    force[2][0],              //
                    force[2][1],              //
                    force[2][2],              //
                    force[3][0],              //
                    force[3][1],              //
                    force[3][2];
                mechanicalObject->resetForce({});

                // Fetch the column of the hessian being affected by the movement of the point
                Eigen::Matrix<VNCS::Real, 12, 1> expectedColumn = expectedHessian.block<12, 1>(0, 3 * nodeId + dirId);
                Eigen::Matrix<VNCS::Real, 12, 1> fdColumn = (forwardForce - backwardForce) / (2.0 * epsilon);
                fdHessian.block<12, 1>(0, 3 * nodeId + dirId) = fdColumn;
            }
        }

        REQUIRE((expectedHessian - fdHessian).isZero(1e-5));
    }
}
