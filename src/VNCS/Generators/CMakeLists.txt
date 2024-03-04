set(SOURCES
    Gen202/Gen202.cpp
    Gen202/SamplingPointsExporter.cpp
    Gen202/GenerateSampler.cpp

    Gen212/Gen212.cpp
    Gen212/SamplingPointsExporter.cpp
    Gen212/GenerateSampler.cpp

    Gen222/Gen222.cpp
    Gen222/VisualSamplersExporter.cpp
    Gen222/SamplingPointsExporter.cpp
    Gen222/GenerateSampler.cpp

    Gen333/Gen333.cpp
    Gen333/VisualSamplersExporter.cpp
    Gen333/SamplingPointsExporter.cpp

    Gen323/Gen323.cpp
    Gen323/SamplingPointsExporter.cpp
    Gen323/GenerateSampler.cpp

    Gen223/Gen223.cpp
    Gen223/SamplingPointsExporter.cpp
    Gen223/GenerateSampler.cpp

    Gen223Barycentric/Gen223.cpp
    Gen223Barycentric/SamplingPointsExporter.cpp
    Gen223Barycentric/GenerateSampler.cpp

    CDT.cpp
    TetraMeshGenerator.cpp
    PMPTriangleIntersection.cpp
    TetraMeshIntersection.cpp
    TetraAABBTree.cpp
    )

set(PUBLIC_HEADERS
    BlendingField.h

    Gen202/Gen202.h
    Gen202/SamplingPointsExporter.h
    Gen202/GenerateSampler.h

    Gen212/Gen212.h
    Gen212/SamplingPointsExporter.h
    Gen212/GenerateSampler.h

    Gen222/Gen222.h
    Gen222/VisualSamplersExporter.h
    Gen222/SamplingPointsExporter.h
    Gen222/GenerateSampler.h

    Gen333/Gen333.h
    Gen333/VisualSamplersExporter.h
    Gen333/SamplingPointsExporter.h

    Gen323/Gen323.h
    Gen323/SamplingPointsExporter.h
    Gen323/GenerateSampler.h

    Gen223/Gen223.h
    Gen223/SamplingPointsExporter.h
    Gen223/GenerateSampler.h

    Gen223Barycentric/Gen223.h
    Gen223Barycentric/SamplingPointsExporter.h
    Gen223Barycentric/GenerateSampler.h

    CDT.h
    TetraMeshGenerator.h
    PMPTriangleIntersection.h
    AdaptiveMeshCriteria.h
    KMeans.h
    IntegrationScheme.h
    TetraMeshIntersection.h
    TetraAABBTree.h
    )

add_library(Generators SHARED ${SOURCES} ${PUBLIC_HEADERS})
set_target_properties(Generators PROPERTIES PUBLIC_HEADER "${PUBLIC_HEADERS}")
target_link_libraries(Generators PUBLIC VNCS)

target_include_directories(Generators PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/../..
        )

#set_source_files_properties(Gen323/SamplingPointsExporter.cpp PROPERTIES COMPILE_FLAGS "-g -O0")

install(TARGETS Generators
        EXPORT VNCSGeneratorsTarget
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Generators/)
