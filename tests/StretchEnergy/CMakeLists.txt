
add_executable(tst_stretchenergy tst_stretchenergy.cpp)
target_link_libraries(tst_stretchenergy VNCS Catch2::Catch2 SofaSimulationGraph)
target_compile_definitions(tst_stretchenergy PRIVATE ASSETS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")
catch_discover_tests(tst_stretchenergy)
