#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <VNCS/Sim3D/BendingForceField.h>
#include <VNCS/Sim2D/BendingForceField.h>
#include <array>

#include <sofa/core/objectmodel/Context.h>
#include <sofa/core/VecId.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/system/FileRepository.h>
#include <sofa/simulation/Node.h>

#include <SofaSimulationGraph/DAGSimulation.h>
#include <SofaSimulationGraph/init.h>
#include <sofa/simulation/MechanicalOperations.h>

#include <SofaBaseMechanics/MechanicalObject.h>
#include <SofaBaseMechanics/MechanicalObject.inl>

using sofa::core::VecId;
using sofa::core::objectmodel::Data;
using sofa::core::objectmodel::New;
using sofa::helper::ReadAccessor;
using sofa::helper::WriteAccessor;
using sofa::simulation::Node;

TEST_CASE("Bending")
{
    SECTION("Sim3D")
    {
        using sofa::defaulttype::Vec3Types;
        using Coord3 = sofa::defaulttype::Vector3;
        using Deriv3 = sofa::defaulttype::Vector3;
        using VecCoord3 = sofa::helper::vector<Coord3>;
        using MechanicalObject3 = sofa::component::container::MechanicalObject<Vec3Types>;
        SECTION("must take into account element size")
        {
            VNCS::Real energyBigSize = 0;
            {
                MechanicalObject3::SPtr mechanicalObject;
                sofa::core::sptr<sofa::simulation::Node> node;
                sofa::core::sptr<VNCS::Sim3D::BendingForceField> bendingEnergy;
                // The graph root node : gravity already exists in a GNode by default
                sofa::simulation::setSimulation(new sofa::simulation::graph::DAGSimulation());
                node = sofa::simulation::getSimulation()->createNewGraph("root");
                node->setGravity(sofa::defaulttype::Vector3(0, 0, 0));

                mechanicalObject = sofa::core::objectmodel::New<MechanicalObject3>();
                WriteAccessor<Data<VecCoord3>> x = *mechanicalObject->write(VecId::restPosition());
                WriteAccessor<Data<VecCoord3>> u = *mechanicalObject->write(VecId::position());
                WriteAccessor<Data<VecCoord3>> dx = *mechanicalObject->write(VecId::dx());
                WriteAccessor<Data<VecCoord3>> df = *mechanicalObject->write(VecId::dforce());
                node->addObject(mechanicalObject);
                mechanicalObject->resize(3);
                // get write access to the position vector of mechanical object DOF
                x[0] = Coord3(-1, 0, 0);
                x[1] = Coord3(0, 0, 0);
                x[2] = Coord3(1, 0, 0);

                u[0] = Coord3(-1, 0, 0);
                u[1] = Coord3(0, 0, 0);
                u[2] = Coord3(-1, 1, 0);

                // Tetrahedron force field
                bendingEnergy = sofa::core::objectmodel::New<VNCS::Sim3D::BendingForceField>();
                node->addObject(bendingEnergy);
                bendingEnergy->setBendingStiffness(1);

                sofa::simulation::graph::init();
                sofa::simulation::getSimulation()->init(node.get());

                energyBigSize = bendingEnergy->getPotentialEnergy(nullptr, *mechanicalObject->read(VecId::position()));
            }

            VNCS::Real energySmallSize = 0;
            {
                MechanicalObject3::SPtr mechanicalObject;
                sofa::core::sptr<sofa::simulation::Node> node;
                sofa::core::sptr<VNCS::Sim3D::BendingForceField> bendingEnergy;
                // The graph root node : gravity already exists in a GNode by default
                sofa::simulation::setSimulation(new sofa::simulation::graph::DAGSimulation());
                node = sofa::simulation::getSimulation()->createNewGraph("root");
                node->setGravity(sofa::defaulttype::Vector3(0, 0, 0));

                mechanicalObject = sofa::core::objectmodel::New<MechanicalObject3>();
                WriteAccessor<Data<VecCoord3>> x = *mechanicalObject->write(VecId::restPosition());
                WriteAccessor<Data<VecCoord3>> u = *mechanicalObject->write(VecId::position());
                WriteAccessor<Data<VecCoord3>> dx = *mechanicalObject->write(VecId::dx());
                WriteAccessor<Data<VecCoord3>> df = *mechanicalObject->write(VecId::dforce());
                node->addObject(mechanicalObject);
                mechanicalObject->resize(3);
                // get write access to the position vector of mechanical object DOF
                x[0] = Coord3(-0.5, 0, 0);
                x[1] = Coord3(0, 0, 0);
                x[2] = Coord3(0.5, 0, 0);

                u[0] = Coord3(-0.5, 0, 0);
                u[1] = Coord3(0, 0, 0);
                u[2] = Coord3(-0.5, 0.5, 0);

                // Tetrahedron force field
                bendingEnergy = sofa::core::objectmodel::New<VNCS::Sim3D::BendingForceField>();
                node->addObject(bendingEnergy);
                bendingEnergy->setBendingStiffness(1);

                sofa::simulation::graph::init();
                sofa::simulation::getSimulation()->init(node.get());

                energySmallSize =
                    bendingEnergy->getPotentialEnergy(nullptr, *mechanicalObject->read(VecId::position()));
            }

            REQUIRE(energyBigSize == 2.0 * energySmallSize);
        }

        static constexpr bool IgnoreHessianTest = false;
        struct BendingEnergyFixture {
            struct Parameter {
                std::array<Coord3, 3> positions;
                bool runHessianTest = true;
            };

            BendingEnergyFixture()
            {
                // The graph root node : gravity already exists in a GNode by default
                sofa::simulation::setSimulation(new sofa::simulation::graph::DAGSimulation());
                node = sofa::simulation::getSimulation()->createNewGraph("root");
                node->setGravity(sofa::defaulttype::Vector3(0, 0, 0));

                mechanicalObject = sofa::core::objectmodel::New<MechanicalObject3>();
                WriteAccessor<Data<VecCoord3>> x = *mechanicalObject->write(VecId::restPosition());
                WriteAccessor<Data<VecCoord3>> dx = *mechanicalObject->write(VecId::dx());
                WriteAccessor<Data<VecCoord3>> df = *mechanicalObject->write(VecId::dforce());
                node->addObject(mechanicalObject);
                mechanicalObject->resize(3);
                // get write access to the position vector of mechanical object DOF
                x[0] = Coord3(-1, 0, 0);
                x[1] = Coord3(0, 0, 0);
                x[2] = Coord3(1, 0, 0);

                // Tetrahedron force field
                bendingEnergy = sofa::core::objectmodel::New<VNCS::Sim3D::BendingForceField>();
                node->addObject(bendingEnergy);
                bendingEnergy->setBendingStiffness(1);

                sofa::simulation::graph::init();
                sofa::simulation::getSimulation()->init(node.get());
            }

            ~BendingEnergyFixture() {}

            MechanicalObject3::SPtr mechanicalObject;
            sofa::core::sptr<sofa::simulation::Node> node;
            sofa::core::sptr<VNCS::Sim3D::BendingForceField> bendingEnergy;
        };

        BendingEnergyFixture model;

        auto param = GENERATE(
            BendingEnergyFixture::Parameter{{Coord3{0, 0, 0}, Coord3{0, 0, 0}, Coord3{0, 0, 0}}, IgnoreHessianTest},
            BendingEnergyFixture::Parameter{{Coord3{0, 0, 0}, Coord3{0, 0, 0}, Coord3{0, 0, 1}}});

        SECTION("Forces")
        {
            WriteAccessor<Data<VecCoord3>> x = *model.mechanicalObject->write(VecId::position());
            x[0] = param.positions[0];
            x[1] = param.positions[1];
            x[2] = param.positions[2];

            sofa::core::MechanicalParams params;
            sofa::simulation::common::MechanicalOperations mop(&params, model.node.get());

            model.bendingEnergy->addForce(nullptr,
                                          *model.mechanicalObject->write(VecId::force()),
                                          *model.mechanicalObject->read(VecId::position()),
                                          *model.mechanicalObject->read(VecId::velocity()));

            ReadAccessor<Data<VecCoord3>> force = *model.mechanicalObject->read(VecId::force());

            std::array<Deriv3, 3> expectedForces = {Deriv3{force[0][0], force[0][1], force[0][2]},
                                                    Deriv3{force[1][0], force[1][1], force[1][2]},
                                                    Deriv3{force[2][0], force[2][1], force[2][2]}};

            Eigen::Map<Eigen::Matrix<VNCS::Real, 9, 1>> expectedForcesMap =
                Eigen::Map<Eigen::Matrix<VNCS::Real, 9, 1>>(&expectedForces[0][0]);
            Eigen::Matrix<VNCS::Real, 9, 1> obtainedForcesMap;
            model.mechanicalObject->resetForce(nullptr);

            // Compute by finite differences
            const auto epsilon = 1e-8;
            for (int nodeId = 0; nodeId < 3; ++nodeId) {
                for (int dirId = 0; dirId < 3; ++dirId) {
                    const auto originalValue = x[nodeId][dirId];
                    x[nodeId][dirId] -= epsilon;
                    const auto backwardEnergy = model.bendingEnergy->getPotentialEnergy(
                        nullptr, *model.mechanicalObject->read(VecId::position()));
                    x[nodeId][dirId] = originalValue;

                    x[nodeId][dirId] += epsilon;
                    const auto forwardEnergy = model.bendingEnergy->getPotentialEnergy(
                        nullptr, *model.mechanicalObject->read(VecId::position()));
                    x[nodeId][dirId] = originalValue;

                    const auto diff = (forwardEnergy - backwardEnergy) / (2.0 * epsilon);
                    obtainedForcesMap(3 * nodeId + dirId) = -diff;
                }
            }
            REQUIRE((expectedForcesMap - obtainedForcesMap).isZero(1e-5));
        }

        SECTION("Hessian")
        {
            if (param.runHessianTest) {
                WriteAccessor<Data<VecCoord3>> x = *model.mechanicalObject->write(VecId::position());
                x[0] = param.positions[0];
                x[1] = param.positions[1];
                x[2] = param.positions[2];

                sofa::core::MechanicalParams params;
                sofa::simulation::common::MechanicalOperations mop(&params, model.node.get());

                model.bendingEnergy->addForce(nullptr,
                                              *model.mechanicalObject->write(VecId::force()),
                                              *model.mechanicalObject->read(VecId::position()),
                                              *model.mechanicalObject->read(VecId::velocity()));

                ReadAccessor<Data<VecCoord3>> force = *model.mechanicalObject->read(VecId::force());
                model.mechanicalObject->resetForce({});

                Eigen::Matrix<VNCS::Real, 9, 9> expectedHessian = model.bendingEnergy->hessians()[0];
                REQUIRE((expectedHessian.transpose().eval() - expectedHessian).isZero());

                // Test addDForce works fine
                WriteAccessor<Data<VecCoord3>> dx = *model.mechanicalObject->write(VecId::dx());
                WriteAccessor<Data<VecCoord3>> df = *model.mechanicalObject->write(VecId::dforce());
                Eigen::Map<Eigen::Matrix<VNCS::Real, 9, 1>> dxMap =
                    Eigen::Map<Eigen::Matrix<VNCS::Real, 9, 1>>(&(dx[0][0]));
                Eigen::Map<Eigen::Matrix<VNCS::Real, 9, 1>> dfMap =
                    Eigen::Map<Eigen::Matrix<VNCS::Real, 9, 1>>(&(df[0][0]));
                dxMap = Eigen::Matrix<VNCS::Real, 9, 1>::Random();

                sofa::core::MechanicalParams mparams;
                mparams.setKFactor(2.0);

                model.bendingEnergy->addDForce(&mparams,
                                               *model.mechanicalObject->write(VecId::dforce()),
                                               *model.mechanicalObject->read(VecId::dx()));

                Eigen::Matrix<VNCS::Real, 9, 1> dfExpected = 2.0 * expectedHessian * dxMap;

                REQUIRE((dfMap - dfExpected).isZero(1e-5));

                Eigen::Matrix<VNCS::Real, 9, 9> fdHessian;
                // Compute by finite differences
                const auto epsilon = 1e-8;
                for (int nodeId = 0; nodeId < 3; ++nodeId) {
                    for (int dirId = 0; dirId < 3; ++dirId) {
                        const auto originalValue = x[nodeId][dirId];
                        x[nodeId][dirId] -= epsilon;
                        model.bendingEnergy->addForce(nullptr,
                                                      *model.mechanicalObject->write(VecId::force()),
                                                      *model.mechanicalObject->read(VecId::position()),
                                                      *model.mechanicalObject->read(VecId::velocity()));
                        ReadAccessor<Data<VecCoord3>> force = *model.mechanicalObject->read(VecId::force());
                        x[nodeId][dirId] = originalValue;

                        Eigen::Matrix<VNCS::Real, 9, 1> backwardForce;
                        backwardForce << force[0][0],  //
                            force[0][1],               //
                            force[0][2],               //
                            force[1][0],               //
                            force[1][1],               //
                            force[1][2],               //
                            force[2][0],               //
                            force[2][1],               //
                            force[2][2];               //
                        model.mechanicalObject->resetForce({});

                        x[nodeId][dirId] += epsilon;
                        model.bendingEnergy->addForce(nullptr,
                                                      *model.mechanicalObject->write(VecId::force()),
                                                      *model.mechanicalObject->read(VecId::position()),
                                                      *model.mechanicalObject->read(VecId::velocity()));
                        force = *model.mechanicalObject->read(VecId::force());
                        x[nodeId][dirId] = originalValue;

                        Eigen::Matrix<VNCS::Real, 9, 1> forwardForce;
                        forwardForce << force[0][0],  //
                            force[0][1],              //
                            force[0][2],              //
                            force[1][0],              //
                            force[1][1],              //
                            force[1][2],              //
                            force[2][0],              //
                            force[2][1],              //
                            force[2][2];              //
                        model.mechanicalObject->resetForce({});

                        // Fetch the column of the hessian being affected by the movement of the point
                        Eigen::Matrix<VNCS::Real, 9, 1> fdColumn = (forwardForce - backwardForce) / (2.0 * epsilon);
                        fdHessian.block<9, 1>(0, 3 * nodeId + dirId) = fdColumn;
                    }
                }

                REQUIRE((expectedHessian - fdHessian).isZero(1e-5));
            }
        }
    }

    SECTION("Sim2D")
    {
        using sofa::defaulttype::Vec2Types;
        using Coord2 = sofa::defaulttype::Vector2;
        using Deriv2 = sofa::defaulttype::Vector2;
        using VecCoord = sofa::helper::vector<Coord2>;
        using MechanicalObject2 = sofa::component::container::MechanicalObject<Vec2Types>;
        static constexpr bool IgnoreHessianTest = false;
        struct BendingEnergyFixture {
            struct Parameter {
                std::array<Coord2, 3> positions;
                bool runHessianTest = true;
            };

            BendingEnergyFixture()
            {
                // The graph root node : gravity already exists in a GNode by default
                sofa::simulation::setSimulation(new sofa::simulation::graph::DAGSimulation());
                node = sofa::simulation::getSimulation()->createNewGraph("root");
                node->setGravity(sofa::defaulttype::Vector3(0, 0, 0));

                mechanicalObject = sofa::core::objectmodel::New<MechanicalObject2>();
                WriteAccessor<Data<VecCoord>> x = *mechanicalObject->write(VecId::restPosition());
                WriteAccessor<Data<VecCoord>> dx = *mechanicalObject->write(VecId::dx());
                WriteAccessor<Data<VecCoord>> df = *mechanicalObject->write(VecId::dforce());
                node->addObject(mechanicalObject);
                mechanicalObject->resize(3);
                // get write access to the position vector of mechanical object DOF
                x[0] = Coord2(-1, 0);
                x[1] = Coord2(0, 0);
                x[2] = Coord2(1, 0);

                // Tetrahedron force field
                bendingEnergy = sofa::core::objectmodel::New<VNCS::Sim2D::BendingForceField>();
                node->addObject(bendingEnergy);
                bendingEnergy->setBendingStiffness(1);

                sofa::simulation::graph::init();
                sofa::simulation::getSimulation()->init(node.get());
            }

            ~BendingEnergyFixture() {}

            MechanicalObject2::SPtr mechanicalObject;
            sofa::core::sptr<sofa::simulation::Node> node;
            sofa::core::sptr<VNCS::Sim2D::BendingForceField> bendingEnergy;
        };

        BendingEnergyFixture model;

        auto param =
            GENERATE(BendingEnergyFixture::Parameter{{Coord2{0, 0}, Coord2{0, 0}, Coord2{0, 0}}, IgnoreHessianTest},
                     BendingEnergyFixture::Parameter{{Coord2{0, 0}, Coord2{0, 0}, Coord2{0, 1}}});

        SECTION("Forces")
        {
            WriteAccessor<Data<VecCoord>> x = *model.mechanicalObject->write(VecId::position());
            x[0] = param.positions[0];
            x[1] = param.positions[1];
            x[2] = param.positions[2];

            sofa::core::MechanicalParams params;
            sofa::simulation::common::MechanicalOperations mop(&params, model.node.get());

            model.bendingEnergy->addForce(nullptr,
                                          *model.mechanicalObject->write(VecId::force()),
                                          *model.mechanicalObject->read(VecId::position()),
                                          *model.mechanicalObject->read(VecId::velocity()));

            ReadAccessor<Data<VecCoord>> force = *model.mechanicalObject->read(VecId::force());

            std::array<Deriv2, 3> expectedForces = {
                Deriv2{force[0][0], force[0][1]}, Deriv2{force[1][0], force[1][1]}, Deriv2{force[2][0], force[2][1]}};

            Eigen::Map<Eigen::Matrix<VNCS::Real, 6, 1>> expectedForcesMap =
                Eigen::Map<Eigen::Matrix<VNCS::Real, 6, 1>>(&expectedForces[0][0]);
            Eigen::Matrix<VNCS::Real, 6, 1> obtainedForcesMap;
            model.mechanicalObject->resetForce(nullptr);

            // Compute by finite differences
            const auto epsilon = 1e-8;
            for (int nodeId = 0; nodeId < 3; ++nodeId) {
                for (int dirId = 0; dirId < 2; ++dirId) {
                    const auto originalValue = x[nodeId][dirId];
                    x[nodeId][dirId] -= epsilon;
                    const auto backwardEnergy = model.bendingEnergy->getPotentialEnergy(
                        nullptr, *model.mechanicalObject->read(VecId::position()));
                    x[nodeId][dirId] = originalValue;

                    x[nodeId][dirId] += epsilon;
                    const auto forwardEnergy = model.bendingEnergy->getPotentialEnergy(
                        nullptr, *model.mechanicalObject->read(VecId::position()));
                    x[nodeId][dirId] = originalValue;

                    const auto diff = (forwardEnergy - backwardEnergy) / (2.0 * epsilon);
                    obtainedForcesMap(2 * nodeId + dirId) = -diff;
                }
            }
            REQUIRE((expectedForcesMap - obtainedForcesMap).isZero(1e-5));
        }

        SECTION("Hessian")
        {
            if (param.runHessianTest) {
                WriteAccessor<Data<VecCoord>> x = *model.mechanicalObject->write(VecId::position());
                x[0] = param.positions[0];
                x[1] = param.positions[1];
                x[2] = param.positions[2];

                sofa::core::MechanicalParams params;
                sofa::simulation::common::MechanicalOperations mop(&params, model.node.get());

                model.bendingEnergy->addForce(nullptr,
                                              *model.mechanicalObject->write(VecId::force()),
                                              *model.mechanicalObject->read(VecId::position()),
                                              *model.mechanicalObject->read(VecId::velocity()));

                ReadAccessor<Data<VecCoord>> force = *model.mechanicalObject->read(VecId::force());
                model.mechanicalObject->resetForce({});

                Eigen::Matrix<VNCS::Real, 6, 6> expectedHessian = model.bendingEnergy->hessians()[0];
                REQUIRE((expectedHessian.transpose().eval() - expectedHessian).isZero());

                // Test addDForce works fine
                WriteAccessor<Data<VecCoord>> dx = *model.mechanicalObject->write(VecId::dx());
                WriteAccessor<Data<VecCoord>> df = *model.mechanicalObject->write(VecId::dforce());
                Eigen::Map<Eigen::Matrix<VNCS::Real, 6, 1>> dxMap =
                    Eigen::Map<Eigen::Matrix<VNCS::Real, 6, 1>>(&(dx[0][0]));
                Eigen::Map<Eigen::Matrix<VNCS::Real, 6, 1>> dfMap =
                    Eigen::Map<Eigen::Matrix<VNCS::Real, 6, 1>>(&(df[0][0]));
                dxMap = Eigen::Matrix<VNCS::Real, 6, 1>::Random();

                sofa::core::MechanicalParams mparams;
                mparams.setKFactor(2.0);

                model.bendingEnergy->addDForce(&mparams,
                                               *model.mechanicalObject->write(VecId::dforce()),
                                               *model.mechanicalObject->read(VecId::dx()));

                Eigen::Matrix<VNCS::Real, 6, 1> dfExpected = 2.0 * expectedHessian * dxMap;

                REQUIRE((dfMap - dfExpected).isZero(1e-5));

                Eigen::Matrix<VNCS::Real, 6, 6> fdHessian;
                // Compute by finite differences
                const auto epsilon = 1e-8;
                for (int nodeId = 0; nodeId < 3; ++nodeId) {
                    for (int dirId = 0; dirId < 2; ++dirId) {
                        const auto originalValue = x[nodeId][dirId];
                        x[nodeId][dirId] -= epsilon;
                        model.bendingEnergy->addForce(nullptr,
                                                      *model.mechanicalObject->write(VecId::force()),
                                                      *model.mechanicalObject->read(VecId::position()),
                                                      *model.mechanicalObject->read(VecId::velocity()));
                        ReadAccessor<Data<VecCoord>> force = *model.mechanicalObject->read(VecId::force());
                        x[nodeId][dirId] = originalValue;

                        Eigen::Matrix<VNCS::Real, 6, 1> backwardForce;
                        backwardForce << force[0][0],  //
                            force[0][1],               //
                            force[1][0],               //
                            force[1][1],               //
                            force[2][0],               //
                            force[2][1];               //
                        model.mechanicalObject->resetForce({});

                        x[nodeId][dirId] += epsilon;
                        model.bendingEnergy->addForce(nullptr,
                                                      *model.mechanicalObject->write(VecId::force()),
                                                      *model.mechanicalObject->read(VecId::position()),
                                                      *model.mechanicalObject->read(VecId::velocity()));
                        force = *model.mechanicalObject->read(VecId::force());
                        x[nodeId][dirId] = originalValue;

                        Eigen::Matrix<VNCS::Real, 6, 1> forwardForce;
                        forwardForce << force[0][0],  //
                            force[0][1],              //
                            force[1][0],              //
                            force[1][1],              //
                            force[2][0],              //
                            force[2][1];              //
                        model.mechanicalObject->resetForce({});

                        // Fetch the column of the hessian being affected by the movement of the point
                        Eigen::Matrix<VNCS::Real, 6, 1> fdColumn = (forwardForce - backwardForce) / (2.0 * epsilon);
                        fdHessian.block<6, 1>(0, 2 * nodeId + dirId) = fdColumn;
                    }
                }

                REQUIRE((expectedHessian - fdHessian).isZero(1e-5));
            }
        }
    }
}
