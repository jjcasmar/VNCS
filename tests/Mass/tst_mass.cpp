#include <gtest/gtest.h>
#include <VNCS/2D/Mass.h>

#include <VNCS/2D/SimCreator.h>

TEST(Mass, 2D)
{
    VNCS::Sim2D::Mass mass;
}

TEST(Mass, 3D)
{
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
