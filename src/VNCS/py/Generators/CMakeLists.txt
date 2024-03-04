set(SOURCES
        )

set(HEADERS
        )

add_library(PyGenerators STATIC ${SOURCES} ${HEADERS})
target_link_libraries(PyGenerators PUBLIC Generators pybind11::module pybind11::headers)

target_include_directories(PyGenerators PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/../..
        )
