set(HEADERS
    Types.h
    Mass.h
    SamplingPoints.h
    Projection.h
    UfMap.h
    UMap.h
    ConcatMap.h
    FinalPositionMap.h
    DeformationGradientTypes.h
    FMap.h
    ConcatMap.h
    Spaces.h
    MechanicalObject.h
    Preconditioner.h
    ConjugateGradientSolver.h
    NBESolver.h
    CoarseStateSetter.h
    Logger.h
    autodiff.h
    DataExtensions.h
    EdgeMesh.h

    Sim2D/Vec2D23D.h
    Sim2D/StVKForceField.h
    Sim2D/FinalPositionMap.h
    Sim2D/BendingForceField.h
    Sim2D/StretchForceField.h

    Sim3D/FinalPositionMap.h
    Sim3D/StableNeoHookean.h
    Sim3D/StVKForceField.h
    Sim3D/TriangleFabricBendingEnergyCommon.h
    Sim3D/StretchForceField.h
    Sim3D/BendingForceField.h
    )

set(SOURCES
    ConjugateGradientSolver.cpp
    NBESolver.cpp

    Sim2D/StVKForceField.cpp
    Sim2D/Vec2D23D.cpp
    Sim2D/BendingForceField.cpp
    Sim2D/StretchForceField.cpp

    Sim3D/StableNeoHookean.cpp
    Sim3D/StVKForceField.cpp
    Sim3D/TriangleFabricBendingEnergyCommon.cpp
    Sim3D/StretchForceField.cpp
    Sim3D/BendingForceField.cpp
    )

add_library(VNCS ${SOURCES} ${HEADERS})
target_include_directories(VNCS PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/../
        )
target_link_libraries(VNCS PUBLIC
        SofaCore
        SofaBaseMechanics
        SofaBaseTopology
        Eigen3::Eigen
        CGAL::CGAL
        range-v3::range-v3
        spdlog::spdlog
        Eigen3::Eigen
        nlohmann_json::nlohmann_json
        spdlog_setup::spdlog_setup
        OpenMP::OpenMP_CXX
        Boost::boost
        )

add_subdirectory(Generators)
add_subdirectory(py)
