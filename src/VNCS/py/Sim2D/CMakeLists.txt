set(SOURCES
        )

set(HEADERS
        )

add_library(PyVNCS2D STATIC ${SOURCES} ${HEADERS})
target_link_libraries(PyVNCS2D PUBLIC VNCS2D pybind11::module pybind11::headers)

target_include_directories(VNCS2D PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/../..
        )
