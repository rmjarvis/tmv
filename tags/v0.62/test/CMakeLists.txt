SET(TEST1_VM TMV_TestVector.cpp TMV_TestMatrix.cpp TMV_TestMatrixArith.cpp TMV_TestMatrixDiv.cpp  TMV_TestMatrixDecomp.cpp)
SET(TEST1_DIAG TMV_TestDiag.cpp TMV_TestDiagArith_A.cpp TMV_TestDiagArith_B1.cpp TMV_TestDiagArith_B2.cpp TMV_TestDiagDiv_A.cpp TMV_TestDiagDiv_B1.cpp TMV_TestDiagDiv_B2.cpp )
SET(TEST1_TRI TMV_TestTri.cpp TMV_TestTriArith_A1.cpp TMV_TestTriArith_A2.cpp TMV_TestTriArith_B1.cpp TMV_TestTriArith_B2.cpp TMV_TestTriArith_C1.cpp TMV_TestTriArith_C2.cpp TMV_TestTriDiv_A1.cpp TMV_TestTriDiv_A2.cpp TMV_TestTriDiv_B1.cpp TMV_TestTriDiv_B2.cpp TMV_TestTriDiv_C1.cpp TMV_TestTriDiv_C2.cpp )
SET(TEST1 TMV_Test1.cpp ${TEST1_VM} ${TEST1_DIAG} ${TEST1_TRI})
SET(TEST1a TMV_Test1a.cpp ${TEST1_VM})
SET(TEST1b TMV_Test1b.cpp ${TEST1_DIAG})
SET(TEST1c TMV_Test1c.cpp ${TEST1_TRI})
SET(TEST2_BAND  TMV_TestBand.cpp TMV_TestBandArith_A.cpp TMV_TestBandArith_B1.cpp TMV_TestBandArith_B2.cpp TMV_TestBandArith_C1.cpp TMV_TestBandArith_C2.cpp TMV_TestBandArith_D1.cpp TMV_TestBandArith_D2.cpp TMV_TestBandDiv.cpp TMV_TestBandDiv_A.cpp TMV_TestBandDiv_B1.cpp TMV_TestBandDiv_B2.cpp TMV_TestBandDiv_C1.cpp TMV_TestBandDiv_C2.cpp TMV_TestBandDiv_D1.cpp TMV_TestBandDiv_D2.cpp TMV_TestBandDecomp.cpp)
SET(TEST2_SYM TMV_TestSym.cpp TMV_TestSymArith_A.cpp TMV_TestSymArith_B1.cpp TMV_TestSymArith_B2.cpp TMV_TestSymArith_C1.cpp TMV_TestSymArith_C2.cpp TMV_TestSymArith_D1.cpp TMV_TestSymArith_D2.cpp TMV_TestSymArith_E1.cpp TMV_TestSymArith_E2.cpp TMV_TestSymDiv.cpp TMV_TestSymDiv_A.cpp TMV_TestSymDiv_B1.cpp TMV_TestSymDiv_B2.cpp TMV_TestSymDiv_C1.cpp TMV_TestSymDiv_C2.cpp TMV_TestSymDiv_D1.cpp TMV_TestSymDiv_D2.cpp TMV_TestSymDiv_E1.cpp TMV_TestSymDiv_E2.cpp TMV_TestSymDecomp.cpp)
SET(TEST2_SYMBAND TMV_TestSymBand.cpp TMV_TestSymBandArith_A.cpp TMV_TestSymBandArith_B1.cpp TMV_TestSymBandArith_B2.cpp TMV_TestSymBandArith_C1.cpp TMV_TestSymBandArith_C2.cpp TMV_TestSymBandArith_D1.cpp TMV_TestSymBandArith_D2.cpp TMV_TestSymBandArith_E1.cpp TMV_TestSymBandArith_E2.cpp TMV_TestSymBandArith_F1.cpp TMV_TestSymBandArith_F2.cpp TMV_TestSymBandDiv.cpp TMV_TestSymBandDiv_A.cpp TMV_TestSymBandDiv_B1.cpp TMV_TestSymBandDiv_B2.cpp TMV_TestSymBandDiv_C1.cpp TMV_TestSymBandDiv_C2.cpp TMV_TestSymBandDiv_D1.cpp TMV_TestSymBandDiv_D2.cpp TMV_TestSymBandDiv_E1.cpp TMV_TestSymBandDiv_E2.cpp TMV_TestSymBandDiv_F1.cpp TMV_TestSymBandDiv_F2.cpp TMV_TestSymBandDecomp.cpp)
SET(TEST2 TMV_Test2.cpp ${TEST2_BAND} ${TEST2_SYM} ${TEST2_SYMBAND})
SET(TEST2a TMV_Test2a.cpp ${TEST2_BAND})
SET(TEST2b TMV_Test2b.cpp ${TEST2_SYM})
SET(TEST2c TMV_Test2c.cpp ${TEST2_SYMBAND})
SET(TEST3_A TMV_TestSmallVector.cpp TMV_TestSmallMatrix.cpp TMV_TestSmallMatrix_Sub.cpp)
SET(TEST3a TMV_Test3a.cpp ${TEST3_A})
SET(TEST3_B TMV_TestSmallMatrixA.cpp TMV_TestSmallMatrixArith_A1.cpp TMV_TestSmallMatrixArith_A2a.cpp TMV_TestSmallMatrixArith_A2b.cpp TMV_TestSmallMatrixArith_A3a.cpp TMV_TestSmallMatrixArith_A3b.cpp TMV_TestSmallMatrixArith_A4a.cpp TMV_TestSmallMatrixArith_A4b.cpp TMV_TestSmallMatrixArith_A5a.cpp TMV_TestSmallMatrixArith_A5b.cpp TMV_TestSmallMatrixArith_A6a.cpp TMV_TestSmallMatrixArith_A6b.cpp TMV_TestSmallMatrixArith_A7.cpp)
SET(TEST3b TMV_Test3b.cpp ${TEST3_B})
SET(TEST3_C TMV_TestSmallMatrixB.cpp TMV_TestSmallMatrixArith_B1.cpp TMV_TestSmallMatrixArith_B2a.cpp TMV_TestSmallMatrixArith_B2b.cpp TMV_TestSmallMatrixArith_B3a.cpp TMV_TestSmallMatrixArith_B3b.cpp TMV_TestSmallMatrixArith_B4a.cpp TMV_TestSmallMatrixArith_B4b.cpp TMV_TestSmallMatrixArith_B4c.cpp TMV_TestSmallMatrixArith_B4d.cpp TMV_TestSmallMatrixArith_B5a.cpp TMV_TestSmallMatrixArith_B5b.cpp TMV_TestSmallMatrixArith_B5c.cpp TMV_TestSmallMatrixArith_B5d.cpp TMV_TestSmallMatrixArith_B6a.cpp TMV_TestSmallMatrixArith_B6b.cpp TMV_TestSmallMatrixArith_B6c.cpp TMV_TestSmallMatrixArith_B6d.cpp TMV_TestSmallMatrixArith_B7.cpp)
SET(TEST3c TMV_Test3c.cpp ${TEST3_C})
SET(TEST3_D TMV_TestSmallMatrixDiv.cpp TMV_TestSmallMatrixDivA.cpp TMV_TestSmallMatrixDiv_A1.cpp TMV_TestSmallMatrixDiv_A2.cpp TMV_TestSmallMatrixDiv_A3.cpp TMV_TestSmallMatrixDiv_A4.cpp TMV_TestSmallMatrixDiv_A5.cpp)
SET(TEST3d TMV_Test3d.cpp ${TEST3_D})
SET(TEST3_E TMV_TestSmallMatrixDivB.cpp TMV_TestSmallMatrixDiv_B1a.cpp TMV_TestSmallMatrixDiv_B1b.cpp TMV_TestSmallMatrixDiv_B1c.cpp TMV_TestSmallMatrixDiv_B2a.cpp TMV_TestSmallMatrixDiv_B2b.cpp TMV_TestSmallMatrixDiv_B3a.cpp TMV_TestSmallMatrixDiv_B3b.cpp TMV_TestSmallMatrixDiv_B4a.cpp TMV_TestSmallMatrixDiv_B4b.cpp TMV_TestSmallMatrixDiv_B5a.cpp TMV_TestSmallMatrixDiv_B5b.cpp)
SET(TEST3e TMV_Test3e.cpp ${TEST3_E})
SET(TEST3 TMV_Test3.cpp ${TEST3_A} ${TEST3_B} ${TEST3_C} ${TEST3_D} ${TEST3_E})


INCLUDE_DIRECTORIES(../include .)

# Test 1
ADD_EXECUTABLE(tmvtest1 ${TEST1})
TARGET_LINK_LIBRARIES(tmvtest1 tmv)

ADD_EXECUTABLE(tmvtest1a ${TEST1a})
TARGET_LINK_LIBRARIES(tmvtest1a tmv)

ADD_EXECUTABLE(tmvtest1b ${TEST1b})
TARGET_LINK_LIBRARIES(tmvtest1b tmv)

ADD_EXECUTABLE(tmvtest1c ${TEST1c})
TARGET_LINK_LIBRARIES(tmvtest1c tmv)

# Test 2
ADD_EXECUTABLE(tmvtest2 ${TEST2})
TARGET_LINK_LIBRARIES(tmvtest2 tmv tmv_symband)

ADD_EXECUTABLE(tmvtest2a ${TEST2a})
TARGET_LINK_LIBRARIES(tmvtest2a tmv tmv_symband)

ADD_EXECUTABLE(tmvtest2b ${TEST2b})
TARGET_LINK_LIBRARIES(tmvtest2b tmv tmv_symband)

ADD_EXECUTABLE(tmvtest2c ${TEST2c})
TARGET_LINK_LIBRARIES(tmvtest2c tmv tmv_symband)

# Test 3
ADD_EXECUTABLE(tmvtest3 ${TEST3})
TARGET_LINK_LIBRARIES(tmvtest3 tmv)

ADD_EXECUTABLE(tmvtest3a ${TEST3a})
TARGET_LINK_LIBRARIES(tmvtest3a tmv)

ADD_EXECUTABLE(tmvtest3b ${TEST3b})
TARGET_LINK_LIBRARIES(tmvtest3b tmv)

ADD_EXECUTABLE(tmvtest3c ${TEST3c})
TARGET_LINK_LIBRARIES(tmvtest3c tmv)

ADD_EXECUTABLE(tmvtest3d ${TEST3d})
TARGET_LINK_LIBRARIES(tmvtest3d tmv)

ADD_EXECUTABLE(tmvtest3e ${TEST3e})
TARGET_LINK_LIBRARIES(tmvtest3e tmv)


