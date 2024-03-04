
add_executable(tst_projection tst_projection.cpp)
target_link_libraries(tst_projection VNCS2D GTest::gtest)
add_test(NAME Projection COMMAND tst_projection)
