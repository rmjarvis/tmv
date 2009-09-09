SET(EX1 Vector.cpp)
SET(EX2 Matrix.cpp)
SET(EX3 Division.cpp)
SET(EX4 BandMatrix.cpp)
SET(EX5 SymMatrix.cpp)

INCLUDE_DIRECTORIES(../include)

ADD_EXECUTABLE(tmvvector ${EX1})
TARGET_LINK_LIBRARIES(tmvvector tmv)

ADD_EXECUTABLE(tmvmatrix ${EX2})
TARGET_LINK_LIBRARIES(tmvmatrix tmv)

ADD_EXECUTABLE(tmvdivision ${EX3})
TARGET_LINK_LIBRARIES(tmvdivision tmv)

ADD_EXECUTABLE(tmvbandmatrix ${EX4})
TARGET_LINK_LIBRARIES(tmvbandmatrix tmv tmv_symband)

ADD_EXECUTABLE(tmvsymmatrix ${EX5})
TARGET_LINK_LIBRARIES(tmvsymmatrix tmv tmv_symband)


