SET(EX1 Vector.cpp)
SET(EX2 Matrix.cpp)
SET(EX3 Division.cpp)
SET(EX4 BandMatrix.cpp)
SET(EX5 SymMatrix.cpp)

INCLUDE_DIRECTORIES(../include)

ADD_EXECUTABLE(vector ${EX1})
TARGET_LINK_LIBRARIES(vector tmv)

ADD_EXECUTABLE(matrix ${EX2})
TARGET_LINK_LIBRARIES(matrix tmv)

ADD_EXECUTABLE(division ${EX3})
TARGET_LINK_LIBRARIES(division tmv)

ADD_EXECUTABLE(bandmatrix ${EX4})
TARGET_LINK_LIBRARIES(bandmatrix tmv tmv_symband)

ADD_EXECUTABLE(symmatrix ${EX5})
TARGET_LINK_LIBRARIES(symmatrix tmv tmv_symband)


