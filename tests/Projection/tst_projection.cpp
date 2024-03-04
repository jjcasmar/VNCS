#include <gtest/gtest.h>
#include <VNCS/2D/Projection.h>
#include <VNCS/2D/Types.h>

TEST(Projection, 2D)
{
    // Create a Clustering matrix
    VNCS::Sim2D::Projection projection;
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
