#define ALIASOK

template <class T> void TestAllVector();
template <class T> void TestAllMatrix();
template <class T> void TestAllMatrixDiv();
template <class T> void TestAllMatrixArith();
template <class T, tmv::StorageType stor> void TestMatrixDecomp();

template <class T> void TestDiagMatrix();
template <class T> void TestDiagMatrixArith_A();
template <class T> void TestDiagMatrixArith_B1();
template <class T> void TestDiagMatrixArith_B2();
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


