
add_executable(tst_mechanicalobject tst_mechanicalobject.cpp)
target_link_libraries(tst_mechanicalobject VNCS2D GTest::gtest)
add_test(NAME MechanicalObject COMMAND tst_mechanicalobject)
