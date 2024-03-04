
add_executable(tst_mass tst_mass.cpp)
target_link_libraries(tst_mass VNCS2D GTest::gtest)
add_test(NAME Mass COMMAND tst_mass)
