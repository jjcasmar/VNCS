#include <gtest/gtest.h>
#include <VNCS/2D/ConcatMap.h>

TEST(ConcatMap, 2D)
{
    // Create a Clustering matrix
    VNCS::Sim2D::ConcatMap concatMap;
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
