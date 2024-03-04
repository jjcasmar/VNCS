set(SOURCES
        )

set(HEADERS
        )

add_library(PyVNCS3D STATIC ${SOURCES} ${HEADERS})
target_link_libraries(PyVNCS3D PUBLIC VNCS3D pybind11::module pybind11::headers)

target_include_directories(VNCS3D PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/../..
        )
