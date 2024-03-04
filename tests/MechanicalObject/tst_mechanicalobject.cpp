#include <gtest/gtest.h>
#include <VNCS/2D/MechanicalObject.h>

TEST(MechanicalObject, 2D)
{
    // Create a Clustering matrix
    auto mo = sofa::core::objectmodel::New<VNCS::Sim2D::MO>();
    mo->resize(10);
    auto restPositions = mo->readRestPositions();
    auto positions = mo->readPositions();
    auto velocities = mo->readVelocities();
    auto forces = mo->readForces();

    EXPECT_EQ(restPositions.size(), 10);
    EXPECT_EQ(positions.size(), 10);
    EXPECT_EQ(velocities.size(), 10);
    EXPECT_EQ(forces.size(), 10);

    // Check we can actually write in the vectors
    auto writerRestPositions = mo->writeRestPositions();
    EXPECT_EQ(writerRestPositions.size(), 10);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
