#include <gtest/gtest.h>
#include <VNCS/2D/Types.h>
#include <VNCS/2D/IntegrationScheme.h>

TEST(IntegrationScheme2D, Gauss1)
{
    // Test we correctly integrate a constant function with gauss1 scheme
    const auto fConstant = [](Eigen::Vector2d p) { return 5.0; };

    {
        VNCS::Sim2D::Kernel::Real result = 0.0;
        for (const auto &p : VNCS::Sim2D::Tri_3::Gauss1::Scheme)
            result += p.weight * fConstant(p.point);

        EXPECT_DOUBLE_EQ(5 * 1.0 / 2.0, result);
    }

    // Test we correctly integrate a linear function with gauss1 scheme
    {
        const auto fLinear = [](Eigen::Vector2d p) { return 5.0 * p[0] + 4.0 * p[1] + 1.0; };

        VNCS::Sim2D::Kernel::Real result = 0.0;
        for (const auto &p : VNCS::Sim2D::Tri_3::Gauss1::Scheme)
            result += p.weight * fLinear(p.point);

        EXPECT_DOUBLE_EQ(2.0, result);
    }

    // Test we fail with a quadratic function
    {
        const auto fLinear = [](Eigen::Vector2d p) { return 5.0 * p[0] * p[0] + 4.0 * p[1] * p[1] + 1.0; };

        VNCS::Sim2D::Kernel::Real result = 0.0;
        for (const auto &p : VNCS::Sim2D::Tri_3::Gauss1::Scheme)
            result += p.weight * fLinear(p.point);

        EXPECT_TRUE(1.25 != result);
    }
}

TEST(IntegrationScheme2D, Gauss3)
{
    // Test we correctly integrate a constant function with gauss1 scheme
    const auto fConstant = [](Eigen::Vector2d p) { return 5.0; };

    {
        VNCS::Sim2D::Kernel::Real result = 0.0;
        for (const auto &p : VNCS::Sim2D::Tri_3::Gauss3::Scheme)
            result += p.weight * fConstant(p.point);

        EXPECT_DOUBLE_EQ(5 * 1.0 / 2.0, result);
    }

    // Test we correctly integrate a linear function with gauss1 scheme
    {
        const auto fLinear = [](Eigen::Vector2d p) { return 5.0 * p[0] + 4.0 * p[1] + 1.0; };

        VNCS::Sim2D::Kernel::Real result = 0.0;
        for (const auto &p : VNCS::Sim2D::Tri_3::Gauss3::Scheme)
            result += p.weight * fLinear(p.point);

        EXPECT_DOUBLE_EQ(2.0, result);
    }

    // Test we fail with a quadratic function
    {
        const auto fQuadratic = [](Eigen::Vector2d p) { return 5.0 * p[0] * p[0] + 4.0 * p[1] * p[1] + 1.0; };

        VNCS::Sim2D::Kernel::Real result = 0.0;
        for (const auto &p : VNCS::Sim2D::Tri_3::Gauss3::Scheme)
            result += p.weight * fQuadratic(p.point);

        EXPECT_DOUBLE_EQ(1.25, result);
    }

    // Test we fail with a quadratic function
    {
        const auto fCubid = [](Eigen::Vector2d p) { return 5.0 * p[0] * p[0] * p[1] + 4.0 * p[1] * p[1] + 1.0; };

        VNCS::Sim2D::Kernel::Real result = 0.0;
        for (const auto &p : VNCS::Sim2D::Tri_3::Gauss3::Scheme)
            result += p.weight * fCubid(p.point);

        EXPECT_TRUE(0.9167 != result);
    }
}

TEST(IntegrationScheme2D, Gauss4)
{
    // Test we correctly integrate a constant function with gauss1 scheme
    const auto fConstant = [](Eigen::Vector2d p) { return 5.0; };

    {
        VNCS::Sim2D::Kernel::Real result = 0.0;
        for (const auto &p : VNCS::Sim2D::Tri_3::Gauss4::Scheme)
            result += p.weight * fConstant(p.point);

        EXPECT_DOUBLE_EQ(5 * 1.0 / 2.0, result);
    }

    // Test we correctly integrate a linear function with gauss1 scheme
    {
        const auto fLinear = [](Eigen::Vector2d p) { return 5.0 * p[0] + 4.0 * p[1] + 1.0; };

        VNCS::Sim2D::Kernel::Real result = 0.0;
        for (const auto &p : VNCS::Sim2D::Tri_3::Gauss4::Scheme)
            result += p.weight * fLinear(p.point);

        EXPECT_DOUBLE_EQ(2.0, result);
    }

    // Test we fail with a quadratic function
    {
        const auto fQuadratic = [](Eigen::Vector2d p) { return 5.0 * p[0] * p[0] + 4.0 * p[1] * p[1] + 1.0; };

        VNCS::Sim2D::Kernel::Real result = 0.0;
        for (const auto &p : VNCS::Sim2D::Tri_3::Gauss4::Scheme)
            result += p.weight * fQuadratic(p.point);

        EXPECT_DOUBLE_EQ(1.25, result);
    }

    // Test we fail with a quadratic function
    {
        const auto fCubid = [](Eigen::Vector2d p) { return 5.0 * p[0] * p[0] * p[1] + 4.0 * p[1] * p[1] + 1.0; };

        VNCS::Sim2D::Kernel::Real result = 0.0;
        for (const auto &p : VNCS::Sim2D::Tri_3::Gauss4::Scheme)
            result += p.weight * fCubid(p.point);

        EXPECT_DOUBLE_EQ(0.91666666666666666, result);
    }
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
