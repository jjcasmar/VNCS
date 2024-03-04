
add_executable(tst_edgemesh tst_edgemesh.cpp)
target_link_libraries(tst_edgemesh VNCS Catch2::Catch2)
target_compile_definitions(tst_edgemesh PRIVATE ASSETS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")
catch_discover_tests(tst_edgemesh)
