cmake_minimum_required(VERSION 3.12)
project(VNCS)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-z,defs")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Eigen3 REQUIRED)
find_package(fmt REQUIRED)
find_package(spdlog REQUIRED)
find_package(spdlog_setup REQUIRED)
find_package(SofaFramework REQUIRED)
find_package(SofaBase REQUIRED)
find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)
find_package(range-v3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(pybind11 REQUIRED)

include(GNUInstallDirs)

add_subdirectory(src)

option(VNCS_ENABLE_TESTS ON "Enable testing")
if (VNCS_ENABLE_TESTS)
    find_package(Catch2 REQUIRED)
    include(CTest)
    include(Catch)
    enable_testing()
    add_subdirectory(tests)
endif()
