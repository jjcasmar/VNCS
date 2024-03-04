set (HEADERS
    SamplingPoints.h
    ConjugateGradient.h
    MechanicalOperations.h

    Sim2D/MechanicalObject.h
    Sim2D/DeformationGradientMO.h
    Sim2D/Projection.h
    Sim2D/UfMap.h
    Sim2D/ConcatMap.h
    Sim2D/DisplacementMap.h
    Sim2D/FinalPositionMap.h
    Sim2D/Vec2D23D.h
    Sim2D/Mass.h
    Sim2D/DeformationGradientMap.h
    Sim2D/EnhaceRelationMap.h
    Sim2D/StretchEnergy.h
    Sim2D/BendingEnergy.h

    Sim3D/MechanicalObject.h
    Sim3D/DeformationGradientMO.h
    Sim3D/Projection.h
    Sim3D/UfMap.h
    Sim3D/ConcatMap.h
    Sim3D/DisplacementMap.h
    Sim3D/FinalPositionMap.h
    Sim3D/Mass.h
    Sim3D/DeformationGradientMap.h
    Sim3D/StVK.h
    Sim3D/StretchEnergy.h
    Sim3D/BendingEnergy.h

    Generators/Gen222.h
    Generators/Gen212.h
    Generators/Gen202.h
    Generators/Gen333.h
    Generators/Gen323.h
    Generators/Gen223.h
    Generators/Gen223Barycentric.h
    Generators/AdaptiveMeshCriteriaFunctor.h
    Generators/BlendingField.h
    )

set (SOURCES
    module.cpp
    SamplingPoints.cpp
    ConjugateGradient.cpp
    MechanicalOperations.cpp

    Sim2D/MechanicalObject.cpp
    Sim2D/DeformationGradientMO.cpp
    Sim2D/Projection.cpp
    Sim2D/UfMap.cpp
    Sim2D/ConcatMap.cpp
    Sim2D/DisplacementMap.cpp
    Sim2D/FinalPositionMap.cpp
    Sim2D/Vec2D23D.cpp
    Sim2D/Mass.cpp
    Sim2D/DeformationGradientMap.cpp
    Sim2D/StVK.cpp
    Sim2D/EnhaceRelationMap.cpp
    Sim2D/StretchEnergy.cpp
    Sim2D/BendingEnergy.cpp

    Sim3D/MechanicalObject.cpp
    Sim3D/DeformationGradientMO.cpp
    Sim3D/Projection.cpp
    Sim3D/UfMap.cpp
    Sim3D/ConcatMap.cpp
    Sim3D/DisplacementMap.cpp
    Sim3D/FinalPositionMap.cpp
    Sim3D/Mass.cpp
    Sim3D/DeformationGradientMap.cpp
    Sim3D/StVK.cpp
    Sim3D/StretchEnergy.cpp
    Sim3D/BendingEnergy.cpp

    Generators/Gen222.cpp
    Generators/Gen212.cpp
    Generators/Gen202.cpp
    Generators/Gen333.cpp
    Generators/Gen323.cpp
    Generators/Gen223.cpp
    Generators/Gen223Barycentric.cpp
    Generators/AdaptiveMeshCriteriaFunctor.cpp
    Generators/BlendingField.cpp
    )

pybind11_add_module(PyVNCS ${SOURCES} ${HEADERS})
target_link_libraries(PyVNCS PUBLIC VNCS Generators pybind11::module pybind11::headers)
