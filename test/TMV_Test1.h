// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#define ALIASOK

#include "tmv/TMV_Base.h"

template <class T> void TestAllVector();
template <class T> void TestAllMatrix();
template <class T> void TestMatrixArith_1();
template <class T> void TestMatrixArith_2();
template <class T> void TestMatrixArith_3();
template <class T> void TestMatrixArith_4();
template <class T> void TestMatrixArith_5();
template <class T> void TestMatrixArith_6();
template <class T> void TestMatrixArith_7();
template <class T> void TestMatrixArith_8();
template <class T> void TestAllMatrixDiv();
template <class T, tmv::StorageType stor> void TestMatrixDecomp();

template <class T> void TestDiagMatrix();
template <class T> void TestDiagMatrixArith_A1();
template <class T> void TestDiagMatrixArith_A2();
template <class T> void TestDiagMatrixArith_A3();
template <class T> void TestDiagMatrixArith_A4();
template <class T> void TestDiagMatrixArith_A5();
template <class T> void TestDiagMatrixArith_A6();
template <class T> void TestDiagMatrixArith_B4a();
template <class T> void TestDiagMatrixArith_B4b();
template <class T> void TestDiagMatrixArith_B5a();
template <class T> void TestDiagMatrixArith_B5b();
template <class T> void TestDiagMatrixArith_B6a();
template <class T> void TestDiagMatrixArith_B6b();
template <class T> void TestDiagDiv();
template <class T> void TestDiagDiv_A();
template <class T> void TestDiagDiv_B1();
template <class T> void TestDiagDiv_B2();

template <class T> void TestTriMatrix();
template <class T> void TestTriMatrixArith_A1();
template <class T> void TestTriMatrixArith_A2();
template <class T> void TestTriMatrixArith_B1();
template <class T> void TestTriMatrixArith_B2();
template <class T> void TestTriMatrixArith_C1();
template <class T> void TestTriMatrixArith_C2();
template <class T> void TestAllTriDiv();
template <class T> void TestTriDiv_A1();
template <class T> void TestTriDiv_A2();
template <class T> void TestTriDiv_B1();
template <class T> void TestTriDiv_B2();
template <class T> void TestTriDiv_C1();
template <class T> void TestTriDiv_C2();
template <class T, tmv::DiagType D> void TestTriDiv();

