#include <gtest/gtest.h>
#include <VNCS/2D/DisplacementMap.h>

TEST(DisplacementMap, 2D)
{
    // Create a Clustering matrix
    VNCS::Sim2D::DisplacementMap uMap;
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
