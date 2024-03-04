
add_executable(tst_stvkenergy tst_stvkenergy.cpp)
target_link_libraries(tst_stvkenergy VNCS Catch2::Catch2 SofaSimulationGraph)
target_compile_definitions(tst_stvkenergy PRIVATE ASSETS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")
catch_discover_tests(tst_stvkenergy)
