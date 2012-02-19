#define ALIASOK

#include "tmv/TMV_Base.h"

template <class T> void TestBandMatrix();
template <class T> void TestBandMatrixArith_A();
template <class T> void TestBandMatrixArith_B1();
template <class T> void TestBandMatrixArith_B2();
template <class T> void TestBandMatrixArith_C1();
template <class T> void TestBandMatrixArith_C2();
template <class T> void TestBandMatrixArith_D1();
template <class T> void TestBandMatrixArith_D2();
template <class T> void TestAllBandDiv();
template <class T> void TestBandDiv(tmv::DivType dt);
template <class T, tmv::StorageType stor> void TestBandDecomp();
template <class T> void TestBandDiv_A(tmv::DivType dt);
template <class T> void TestBandDiv_B1(tmv::DivType dt);
template <class T> void TestBandDiv_B2(tmv::DivType dt);
template <class T> void TestBandDiv_C1(tmv::DivType dt);
template <class T> void TestBandDiv_C2(tmv::DivType dt);
template <class T> void TestBandDiv_D1(tmv::DivType dt);
template <class T> void TestBandDiv_D2(tmv::DivType dt);

enum PosDefCode { PosDef, InDef, Sing };
inline std::string PDLabel(PosDefCode pdc)
{
    if (pdc == PosDef) return "Positive Definite";
    else if (pdc == InDef) return "Indefinite";
    else return "Singular";
}

template <class T> void TestSymMatrix();
template <class T> void TestSymMatrixArith_A();
template <class T> void TestSymMatrixArith_B1();
template <class T> void TestSymMatrixArith_B2();
template <class T> void TestSymMatrixArith_C1();
template <class T> void TestSymMatrixArith_C2();
template <class T> void TestSymMatrixArith_D1();
template <class T> void TestSymMatrixArith_D2();
template <class T> void TestSymMatrixArith_E1();
template <class T> void TestSymMatrixArith_E2();
template <class T> void TestAllSymDiv();
template <class T> void TestSymDiv(tmv::DivType dt, PosDefCode pc);
template <class T, tmv::UpLoType uplo, tmv::StorageType stor>
void TestHermDecomp();
template <class T, tmv::UpLoType uplo, tmv::StorageType stor>
void TestSymDecomp();
template <class T, tmv::StorageType stor> void TestPolar();
template <class T> void TestSymDiv_A(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymDiv_B1(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymDiv_B2(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymDiv_C1(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymDiv_C2(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymDiv_D1(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymDiv_D2(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymDiv_E1(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymDiv_E2(tmv::DivType dt, PosDefCode pc);

template <class T> void TestSymBandMatrix();
template <class T> void TestSymBandMatrixArith_A();
template <class T> void TestSymBandMatrixArith_B1();
template <class T> void TestSymBandMatrixArith_B2();
template <class T> void TestSymBandMatrixArith_C1();
template <class T> void TestSymBandMatrixArith_C2();
template <class T> void TestSymBandMatrixArith_D1();
template <class T> void TestSymBandMatrixArith_D2();
template <class T> void TestSymBandMatrixArith_E1();
template <class T> void TestSymBandMatrixArith_E2();
template <class T> void TestSymBandMatrixArith_F1();
template <class T> void TestSymBandMatrixArith_F2();
template <class T> void TestAllSymBandDiv();
template <class T, tmv::UpLoType uplo, tmv::StorageType stor>
void TestHermBandDecomp();
template <class T, tmv::UpLoType uplo, tmv::StorageType stor>
void TestSymBandDecomp();
template <class T> void TestSymBandDiv(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_A(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_B1(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_B2(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_C1(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_C2(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_D1(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_D2(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_E1(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_E2(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_F1(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_F2(tmv::DivType dt, PosDefCode pc);

