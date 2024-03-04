
add_executable(tst_concatmap tst_concatmap.cpp)
target_link_libraries(tst_concatmap VNCS2D GTest::gtest)
add_test(NAME ConcatMap COMMAND tst_concatmap)
