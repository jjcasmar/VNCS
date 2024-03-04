
add_executable(tst_bendingenergy tst_bendingenergy.cpp)
target_link_libraries(tst_bendingenergy VNCS Catch2::Catch2 SofaSimulationGraph)
target_compile_definitions(tst_bendingenergy PRIVATE ASSETS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")
catch_discover_tests(tst_bendingenergy)
