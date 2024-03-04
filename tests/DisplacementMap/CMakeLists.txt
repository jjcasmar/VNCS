
add_executable(tst_displacementmap tst_displacementmap.cpp)
target_link_libraries(tst_displacementmap VNCS2D GTest::gtest)
add_test(NAME DisplacementMap COMMAND tst_displacementmap)
