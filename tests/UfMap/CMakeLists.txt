
add_executable(tst_ufmap tst_ufmap.cpp)
target_link_libraries(tst_ufmap VNCS2D GTest::gtest)
add_test(NAME UfMap COMMAND tst_ufmap)
