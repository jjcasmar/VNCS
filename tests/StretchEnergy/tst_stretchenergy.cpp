#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <VNCS/Sim3D/StretchForceField.h>
#include <VNCS/Sim2D/StretchForceField.h>
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

TEST_CASE("StretchEnergy")
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
                sofa::core::sptr<VNCS::Sim3D::StretchForceField> bendingEnergy;
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
                u[2] = Coord3(1, 0, 0);

                // Tetrahedron force field
                bendingEnergy = sofa::core::objectmodel::New<VNCS::Sim3D::StretchForceField>();
                node->addObject(bendingEnergy);
                bendingEnergy->setStretchStiffness(1);

                sofa::simulation::graph::init();
                sofa::simulation::getSimulation()->init(node.get());

                energyBigSize = bendingEnergy->getPotentialEnergy(nullptr, *mechanicalObject->read(VecId::position()));
            }

            VNCS::Real energySmallSize = 0;
            {
                MechanicalObject3::SPtr mechanicalObject;
                sofa::core::sptr<sofa::simulation::Node> node;
                sofa::core::sptr<VNCS::Sim3D::StretchForceField> bendingEnergy;
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
                u[2] = Coord3(0.5, 0, 0);

                // Tetrahedron force field
                bendingEnergy = sofa::core::objectmodel::New<VNCS::Sim3D::StretchForceField>();
                node->addObject(bendingEnergy);
                bendingEnergy->setStretchStiffness(1);

                sofa::simulation::graph::init();
                sofa::simulation::getSimulation()->init(node.get());

                energySmallSize =
                    bendingEnergy->getPotentialEnergy(nullptr, *mechanicalObject->read(VecId::position()));
            }

            REQUIRE(energyBigSize == 2.0 * energySmallSize);
        }

        static constexpr bool IgnoreHessianTest = false;
        struct Parameter {
            std::array<Coord3, 2> positions;
        };

        MechanicalObject3::SPtr mechanicalObject;
        sofa::core::sptr<sofa::simulation::Node> node;
        sofa::core::sptr<VNCS::Sim3D::StretchForceField> energy;

        // The graph root node : gravity already exists in a GNode by default
        sofa::simulation::setSimulation(new sofa::simulation::graph::DAGSimulation());
        node = sofa::simulation::getSimulation()->createNewGraph("root");
        node->setGravity(sofa::defaulttype::Vector3(0, 0, 0));

        mechanicalObject = sofa::core::objectmodel::New<MechanicalObject3>();
        WriteAccessor<Data<VecCoord3>> x = *mechanicalObject->write(VecId::restPosition());
        WriteAccessor<Data<VecCoord3>> dx = *mechanicalObject->write(VecId::dx());
        WriteAccessor<Data<VecCoord3>> df = *mechanicalObject->write(VecId::dforce());
        node->addObject(mechanicalObject);
        mechanicalObject->resize(2);
        // get write access to the position vector of mechanical object DOF
        x[0] = Coord3(-1, 0, 0);
        x[1] = Coord3(1, 0, 0);

        // Tetrahedron force field
        energy = sofa::core::objectmodel::New<VNCS::Sim3D::StretchForceField>();
        node->addObject(energy);
        energy->setStretchStiffness(0.9);

        sofa::simulation::graph::init();
        sofa::simulation::getSimulation()->init(node.get());

        auto param = GENERATE(Parameter{{Coord3{0, 0, 0}, Coord3{0, 0, 0}}},
                              Parameter{{Coord3{0, 0, 0}, Coord3{1, 0, 0}}},
                              Parameter{{Coord3{0, 0, 0}, Coord3{0, 1, 0}}},
                              Parameter{{Coord3{0, 0, 0}, Coord3{0, 0, 1}}},
                              Parameter{{Coord3{1, 0, 0}, Coord3{0, 0, 0}}},
                              Parameter{{Coord3{0, 1, 0}, Coord3{0, 0, 0}}},
                              Parameter{{Coord3{0, 0, 1}, Coord3{0, 0, 0}}},
                              Parameter{{Coord3{1, 0, 0}, Coord3{1, 0, 0}}},
                              Parameter{{Coord3{0, 1, 0}, Coord3{0, 1, 0}}},
                              Parameter{{Coord3{0, 0, 1}, Coord3{0, 0, 1}}});

        SECTION("Forces")
        {
            WriteAccessor<Data<VecCoord3>> x = *mechanicalObject->write(VecId::position());
            x[0] = param.positions[0];
            x[1] = param.positions[1];

            sofa::core::MechanicalParams params;
            sofa::simulation::common::MechanicalOperations mop(&params, node.get());

            energy->addForce(nullptr,
                             *mechanicalObject->write(VecId::force()),
                             *mechanicalObject->read(VecId::position()),
                             *mechanicalObject->read(VecId::velocity()));

            ReadAccessor<Data<VecCoord3>> force = *mechanicalObject->read(VecId::force());

            std::array<Deriv3, 2> expectedForces = {Deriv3{force[0][0], force[0][1], force[0][2]},
                                                    Deriv3{force[1][0], force[1][1], force[1][2]}};

            Eigen::Map<Eigen::Matrix<VNCS::Real, 6, 1>> expectedForcesMap =
                Eigen::Map<Eigen::Matrix<VNCS::Real, 6, 1>>(&expectedForces[0][0]);
            Eigen::Matrix<VNCS::Real, 6, 1> obtainedForcesMap;
            mechanicalObject->resetForce(nullptr);

            // Compute by finite differences
            const auto epsilon = 1e-8;
            for (int nodeId = 0; nodeId < 2; ++nodeId) {
                for (int dirId = 0; dirId < 3; ++dirId) {
                    const auto originalValue = x[nodeId][dirId];
                    x[nodeId][dirId] -= epsilon;
                    const auto backwardEnergy =
                        energy->getPotentialEnergy(nullptr, *mechanicalObject->read(VecId::position()));
                    x[nodeId][dirId] = originalValue;

                    x[nodeId][dirId] += epsilon;
                    const auto forwardEnergy =
                        energy->getPotentialEnergy(nullptr, *mechanicalObject->read(VecId::position()));
                    x[nodeId][dirId] = originalValue;

                    const auto diff = (forwardEnergy - backwardEnergy) / (2.0 * epsilon);
                    obtainedForcesMap(3 * nodeId + dirId) = -diff;
                }
            }
            REQUIRE((expectedForcesMap - obtainedForcesMap).isZero(1e-5));
        }

        SECTION("Hessian")
        {
            WriteAccessor<Data<VecCoord3>> x = *mechanicalObject->write(VecId::position());
            x[0] = param.positions[0];
            x[1] = param.positions[1];

            sofa::core::MechanicalParams params;
            sofa::simulation::common::MechanicalOperations mop(&params, node.get());

            energy->addForce(nullptr,
                             *mechanicalObject->write(VecId::force()),
                             *mechanicalObject->read(VecId::position()),
                             *mechanicalObject->read(VecId::velocity()));

            ReadAccessor<Data<VecCoord3>> force = *mechanicalObject->read(VecId::force());
            mechanicalObject->resetForce({});

            Eigen::Matrix<VNCS::Real, 6, 6> expectedHessian = energy->hessians()[0];
            REQUIRE((expectedHessian.transpose().eval() - expectedHessian).isZero());

            // Test addDForce works fine
            WriteAccessor<Data<VecCoord3>> dx = *mechanicalObject->write(VecId::dx());
            WriteAccessor<Data<VecCoord3>> df = *mechanicalObject->write(VecId::dforce());
            Eigen::Map<Eigen::Matrix<VNCS::Real, 6, 1>> dxMap =
                Eigen::Map<Eigen::Matrix<VNCS::Real, 6, 1>>(&(dx[0][0]));
            Eigen::Map<Eigen::Matrix<VNCS::Real, 6, 1>> dfMap =
                Eigen::Map<Eigen::Matrix<VNCS::Real, 6, 1>>(&(df[0][0]));
            dxMap = Eigen::Matrix<VNCS::Real, 6, 1>::Random();

            sofa::core::MechanicalParams mparams;
            mparams.setKFactor(2.0);

            energy->addDForce(
                &mparams, *mechanicalObject->write(VecId::dforce()), *mechanicalObject->read(VecId::dx()));

            Eigen::Matrix<VNCS::Real, 6, 1> dfExpected = 2.0 * expectedHessian * dxMap;

            REQUIRE((dfMap - dfExpected).isZero(1e-5));

            Eigen::Matrix<VNCS::Real, 6, 6> fdHessian;
            // Compute by finite differences
            const auto epsilon = 1e-8;
            for (int nodeId = 0; nodeId < 2; ++nodeId) {
                for (int dirId = 0; dirId < 3; ++dirId) {
                    const auto originalValue = x[nodeId][dirId];
                    x[nodeId][dirId] -= epsilon;
                    energy->addForce(nullptr,
                                     *mechanicalObject->write(VecId::force()),
                                     *mechanicalObject->read(VecId::position()),
                                     *mechanicalObject->read(VecId::velocity()));
                    ReadAccessor<Data<VecCoord3>> force = *mechanicalObject->read(VecId::force());
                    x[nodeId][dirId] = originalValue;

                    Eigen::Matrix<VNCS::Real, 6, 1> backwardForce;
                    backwardForce << force[0][0],  //
                        force[0][1],               //
                        force[0][2],               //
                        force[1][0],               //
                        force[1][1],               //
                        force[1][2];               //
                    mechanicalObject->resetForce({});

                    x[nodeId][dirId] += epsilon;
                    energy->addForce(nullptr,
                                     *mechanicalObject->write(VecId::force()),
                                     *mechanicalObject->read(VecId::position()),
                                     *mechanicalObject->read(VecId::velocity()));
                    force = *mechanicalObject->read(VecId::force());
                    x[nodeId][dirId] = originalValue;

                    Eigen::Matrix<VNCS::Real, 6, 1> forwardForce;
                    forwardForce << force[0][0],  //
                        force[0][1],              //
                        force[0][2],              //
                        force[1][0],              //
                        force[1][1],              //
                        force[1][2];              //
                    mechanicalObject->resetForce({});

                    // Fetch the column of the hessian being affected by the movement of the point
                    Eigen::Matrix<VNCS::Real, 6, 1> fdColumn = (forwardForce - backwardForce) / (2.0 * epsilon);
                    fdHessian.block<6, 1>(0, 3 * nodeId + dirId) = fdColumn;
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

        SECTION("must take into account element size")
        {
            VNCS::Real energyBigSize = 0;
            {
                MechanicalObject2::SPtr mechanicalObject;
                sofa::core::sptr<sofa::simulation::Node> node;
                sofa::core::sptr<VNCS::Sim2D::StretchForceField> bendingEnergy;
                // The graph root node : gravity already exists in a GNode by default
                sofa::simulation::setSimulation(new sofa::simulation::graph::DAGSimulation());
                node = sofa::simulation::getSimulation()->createNewGraph("root");
                node->setGravity(sofa::defaulttype::Vector3(0, 0, 0));

                mechanicalObject = sofa::core::objectmodel::New<MechanicalObject2>();
                WriteAccessor<Data<VecCoord>> x = *mechanicalObject->write(VecId::restPosition());
                WriteAccessor<Data<VecCoord>> u = *mechanicalObject->write(VecId::position());
                WriteAccessor<Data<VecCoord>> dx = *mechanicalObject->write(VecId::dx());
                WriteAccessor<Data<VecCoord>> df = *mechanicalObject->write(VecId::dforce());
                node->addObject(mechanicalObject);
                mechanicalObject->resize(3);
                // get write access to the position vector of mechanical object DOF
                x[0] = Coord2(-1, 0);
                x[1] = Coord2(0, 0);
                x[2] = Coord2(1, 0);

                u[0] = Coord2(-1, 0);
                u[1] = Coord2(0, 0);
                u[2] = Coord2(1, 0);

                // Tetrahedron force field
                bendingEnergy = sofa::core::objectmodel::New<VNCS::Sim2D::StretchForceField>();
                node->addObject(bendingEnergy);
                bendingEnergy->setStretchStiffness(1);

                sofa::simulation::graph::init();
                sofa::simulation::getSimulation()->init(node.get());

                energyBigSize = bendingEnergy->getPotentialEnergy(nullptr, *mechanicalObject->read(VecId::position()));
            }

            VNCS::Real energySmallSize = 0;
            {
                MechanicalObject2::SPtr mechanicalObject;
                sofa::core::sptr<sofa::simulation::Node> node;
                sofa::core::sptr<VNCS::Sim2D::StretchForceField> bendingEnergy;
                // The graph root node : gravity already exists in a GNode by default
                sofa::simulation::setSimulation(new sofa::simulation::graph::DAGSimulation());
                node = sofa::simulation::getSimulation()->createNewGraph("root");
                node->setGravity(sofa::defaulttype::Vector3(0, 0, 0));

                mechanicalObject = sofa::core::objectmodel::New<MechanicalObject2>();
                WriteAccessor<Data<VecCoord>> x = *mechanicalObject->write(VecId::restPosition());
                WriteAccessor<Data<VecCoord>> u = *mechanicalObject->write(VecId::position());
                WriteAccessor<Data<VecCoord>> dx = *mechanicalObject->write(VecId::dx());
                WriteAccessor<Data<VecCoord>> df = *mechanicalObject->write(VecId::dforce());
                node->addObject(mechanicalObject);
                mechanicalObject->resize(3);
                // get write access to the position vector of mechanical object DOF
                x[0] = Coord2(-0.5, 0);
                x[1] = Coord2(0, 0);
                x[2] = Coord2(0.5, 0);

                u[0] = Coord2(-0.5, 0);
                u[1] = Coord2(0, 0);
                u[2] = Coord2(0.5, 0);

                // Tetrahedron force field
                bendingEnergy = sofa::core::objectmodel::New<VNCS::Sim2D::StretchForceField>();
                node->addObject(bendingEnergy);
                bendingEnergy->setStretchStiffness(1);

                sofa::simulation::graph::init();
                sofa::simulation::getSimulation()->init(node.get());

                energySmallSize =
                    bendingEnergy->getPotentialEnergy(nullptr, *mechanicalObject->read(VecId::position()));
            }

            REQUIRE(energyBigSize == 2.0 * energySmallSize);
        }

        struct Parameter {
            std::array<Coord2, 2> positions;
        };

        MechanicalObject2::SPtr mechanicalObject;
        sofa::core::sptr<sofa::simulation::Node> node;
        sofa::core::sptr<VNCS::Sim2D::StretchForceField> energy;

        // The graph root node : gravity already exists in a GNode by default
        sofa::simulation::setSimulation(new sofa::simulation::graph::DAGSimulation());
        node = sofa::simulation::getSimulation()->createNewGraph("root");
        node->setGravity(sofa::defaulttype::Vector3(0, 0, 0));

        mechanicalObject = sofa::core::objectmodel::New<MechanicalObject2>();
        WriteAccessor<Data<VecCoord>> x = *mechanicalObject->write(VecId::restPosition());
        WriteAccessor<Data<VecCoord>> dx = *mechanicalObject->write(VecId::dx());
        WriteAccessor<Data<VecCoord>> df = *mechanicalObject->write(VecId::dforce());
        node->addObject(mechanicalObject);
        mechanicalObject->resize(2);
        // get write access to the position vector of mechanical object DOF
        x[0] = Coord2(-1, 0);
        x[1] = Coord2(1, 0);

        // Tetrahedron force field
        energy = sofa::core::objectmodel::New<VNCS::Sim2D::StretchForceField>();
        node->addObject(energy);
        energy->setStretchStiffness(0.9);

        sofa::simulation::graph::init();
        sofa::simulation::getSimulation()->init(node.get());

        auto param = GENERATE(Parameter{{Coord2{0, 0}, Coord2{0, 0}}},
                              Parameter{{Coord2{0, 0}, Coord2{1, 0}}},
                              Parameter{{Coord2{0, 0}, Coord2{0, 1}}},
                              Parameter{{Coord2{0, 0}, Coord2{0, 0}}},
                              Parameter{{Coord2{1, 0}, Coord2{0, 0}}},
                              Parameter{{Coord2{0, 1}, Coord2{0, 0}}},
                              Parameter{{Coord2{0, 0}, Coord2{0, 0}}},
                              Parameter{{Coord2{1, 0}, Coord2{1, 0}}},
                              Parameter{{Coord2{0, 1}, Coord2{0, 1}}},
                              Parameter{{Coord2{0, 0}, Coord2{0, 0}}});

        SECTION("Forces")
        {
            WriteAccessor<Data<VecCoord>> x = *mechanicalObject->write(VecId::position());
            x[0] = param.positions[0];
            x[1] = param.positions[1];

            sofa::core::MechanicalParams params;
            sofa::simulation::common::MechanicalOperations mop(&params, node.get());

            energy->addForce(nullptr,
                             *mechanicalObject->write(VecId::force()),
                             *mechanicalObject->read(VecId::position()),
                             *mechanicalObject->read(VecId::velocity()));

            ReadAccessor<Data<VecCoord>> force = *mechanicalObject->read(VecId::force());

            std::array<Deriv2, 2> expectedForces = {Deriv2{force[0][0], force[0][1]}, Deriv2{force[1][0], force[1][1]}};

            Eigen::Map<Eigen::Matrix<VNCS::Real, 4, 1>> expectedForcesMap =
                Eigen::Map<Eigen::Matrix<VNCS::Real, 4, 1>>(&expectedForces[0][0]);
            Eigen::Matrix<VNCS::Real, 4, 1> obtainedForcesMap;
            mechanicalObject->resetForce(nullptr);

            // Compute by finite differences
            const auto epsilon = 1e-8;
            for (int nodeId = 0; nodeId < 2; ++nodeId) {
                for (int dirId = 0; dirId < 2; ++dirId) {
                    const auto originalValue = x[nodeId][dirId];
                    x[nodeId][dirId] -= epsilon;
                    const auto backwardEnergy =
                        energy->getPotentialEnergy(nullptr, *mechanicalObject->read(VecId::position()));
                    x[nodeId][dirId] = originalValue;

                    x[nodeId][dirId] += epsilon;
                    const auto forwardEnergy =
                        energy->getPotentialEnergy(nullptr, *mechanicalObject->read(VecId::position()));
                    x[nodeId][dirId] = originalValue;

                    const auto diff = (forwardEnergy - backwardEnergy) / (2.0 * epsilon);
                    obtainedForcesMap(2 * nodeId + dirId) = -diff;
                }
            }
            REQUIRE((expectedForcesMap - obtainedForcesMap).isZero(1e-5));
        }

        SECTION("Hessian")
        {
            WriteAccessor<Data<VecCoord>> x = *mechanicalObject->write(VecId::position());
            x[0] = param.positions[0];
            x[1] = param.positions[1];

            sofa::core::MechanicalParams params;
            sofa::simulation::common::MechanicalOperations mop(&params, node.get());

            energy->addForce(nullptr,
                             *mechanicalObject->write(VecId::force()),
                             *mechanicalObject->read(VecId::position()),
                             *mechanicalObject->read(VecId::velocity()));

            ReadAccessor<Data<VecCoord>> force = *mechanicalObject->read(VecId::force());
            mechanicalObject->resetForce({});

            Eigen::Matrix<VNCS::Real, 4, 4> expectedHessian = energy->hessians()[0];
            REQUIRE((expectedHessian.transpose().eval() - expectedHessian).isZero());

            // Test addDForce works fine
            WriteAccessor<Data<VecCoord>> dx = *mechanicalObject->write(VecId::dx());
            WriteAccessor<Data<VecCoord>> df = *mechanicalObject->write(VecId::dforce());
            Eigen::Map<Eigen::Matrix<VNCS::Real, 4, 1>> dxMap =
                Eigen::Map<Eigen::Matrix<VNCS::Real, 4, 1>>(&(dx[0][0]));
            Eigen::Map<Eigen::Matrix<VNCS::Real, 4, 1>> dfMap =
                Eigen::Map<Eigen::Matrix<VNCS::Real, 4, 1>>(&(df[0][0]));
            dxMap = Eigen::Matrix<VNCS::Real, 4, 1>::Random();

            sofa::core::MechanicalParams mparams;
            mparams.setKFactor(2.0);

            energy->addDForce(
                &mparams, *mechanicalObject->write(VecId::dforce()), *mechanicalObject->read(VecId::dx()));

            Eigen::Matrix<VNCS::Real, 4, 1> dfExpected = 2.0 * expectedHessian * dxMap;

            REQUIRE((dfMap - dfExpected).isZero(1e-5));

            Eigen::Matrix<VNCS::Real, 4, 4> fdHessian;
            // Compute by finite differences
            const auto epsilon = 1e-8;
            for (int nodeId = 0; nodeId < 2; ++nodeId) {
                for (int dirId = 0; dirId < 2; ++dirId) {
                    const auto originalValue = x[nodeId][dirId];
                    x[nodeId][dirId] -= epsilon;
                    energy->addForce(nullptr,
                                     *mechanicalObject->write(VecId::force()),
                                     *mechanicalObject->read(VecId::position()),
                                     *mechanicalObject->read(VecId::velocity()));
                    ReadAccessor<Data<VecCoord>> force = *mechanicalObject->read(VecId::force());
                    x[nodeId][dirId] = originalValue;

                    Eigen::Matrix<VNCS::Real, 4, 1> backwardForce;
                    backwardForce << force[0][0],  //
                        force[0][1],               //
                        force[1][0],               //
                        force[1][1];               //
                    mechanicalObject->resetForce({});

                    x[nodeId][dirId] += epsilon;
                    energy->addForce(nullptr,
                                     *mechanicalObject->write(VecId::force()),
                                     *mechanicalObject->read(VecId::position()),
                                     *mechanicalObject->read(VecId::velocity()));
                    force = *mechanicalObject->read(VecId::force());
                    x[nodeId][dirId] = originalValue;

                    Eigen::Matrix<VNCS::Real, 4, 1> forwardForce;
                    forwardForce << force[0][0],  //
                        force[0][1],              //
                        force[1][0],              //
                        force[1][1];              //
                    mechanicalObject->resetForce({});

                    // Fetch the column of the hessian being affected by the movement of the point
                    Eigen::Matrix<VNCS::Real, 4, 1> fdColumn = (forwardForce - backwardForce) / (2.0 * epsilon);
                    fdHessian.block<4, 1>(0, 2 * nodeId + dirId) = fdColumn;
                }

                REQUIRE((expectedHessian - fdHessian).isZero(1e-5));
            }
        }
    }
}
