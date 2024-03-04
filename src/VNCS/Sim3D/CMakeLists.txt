set(SOURCES
    StVKForceField.cpp
    StableNeoHookean.cpp
        )

set(HEADERS
        )

set(PUBLIC_HEADERS
    StVKForceField.h
    FinalPositionMap.h
    StableNeoHookean.h
        )


add_library(VNCS3D SHARED ${SOURCES} ${HEADERS} ${PUBLIC_HEADERS})
set_target_properties(VNCS3D PROPERTIES PUBLIC_HEADER "${PUBLIC_HEADERS}")
target_link_libraries(VNCS3D
    PUBLIC
        VNCS
    PRIVATE
        OpenMP::OpenMP_CXX
        )

install(TARGETS VNCS3D
        EXPORT VNCS2DTarget
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/VNCS/3D)
