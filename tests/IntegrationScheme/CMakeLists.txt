
add_executable(tst_integrationscheme tst_integrationscheme.cpp)
target_link_libraries(tst_integrationscheme VNCS2D GTest::gtest)
add_test(NAME UfMap COMMAND tst_integrationscheme)
