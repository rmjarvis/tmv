
#define TESTDIV
#define TESTDIAG
#define TESTTRI
#define TESTBAND
#define TESTSYM

#include "TMV.h"
//#include "TMV_Test.h"
using tmv::Matrix;

#ifdef TESTDIAG
#include "TMV_Diag.h"
#endif

#ifdef TESTTRI
#include "TMV_Tri.h"
using tmv::UpperTriMatrix;
#endif

#ifdef TESTBAND
#include "TMV_Band.h"
using tmv::BandMatrix;
#endif

#ifdef TESTSYM
#include "TMV_Sym.h"
using tmv::SymMatrix;
#endif

#include "MemDebug.h"

#ifdef MEMDEBUG
AllocList* allocList=0;
#endif

template <class T> extern void TestAllVector();
template <class T> extern void TestAllMatrix();
template <class T> extern void TestAllMatrixDiv();
extern void TestPermutation();
template <class T> extern void TestDiagMatrix();
template <class T> extern void TestDiagDiv();
template <class T> extern void TestTriMatrix();
template <class T> extern void TestAllTriDiv();
template <class T> extern void TestBandMatrix();
template <class T> extern void TestAllBandDiv();
template <class T> extern void TestSymMatrix();
template <class T> extern void TestAllSymDiv();

int main() {

#ifdef MEMDEBUG
  atexit(&DumpUnfreed);
#endif

//  TestAllMatrixDiv<double>();
//#define SKIPREST

#ifndef SKIPREST
  TestAllVector<double>();
  TestAllMatrix<double>();
#ifdef TESTDIAG
  TestDiagMatrix<double>();
  TestDiagDiv<double>();
#endif
#ifdef TESTTRI
  TestTriMatrix<double>();
  TestAllTriDiv<double>();
#endif
#ifdef TESTDIV
  TestAllMatrixDiv<double>();
#endif
#ifdef TESTBAND
  TestBandMatrix<double>();
  TestAllBandDiv<double>();
#endif
#ifdef TESTSYM
  TestSymMatrix<double>();
  TestAllSymDiv<double>();
#endif

#ifndef NOFLOAT
  TestAllVector<float>();
  TestAllMatrix<float>();
#ifdef TESTDIV
  TestAllMatrixDiv<float>();
#endif
#ifdef TESTDIAG
  TestDiagMatrix<float>();
  TestDiagDiv<float>();
#endif
#ifdef TESTTRI
  TestTriMatrix<float>();
  TestAllTriDiv<float>();
#endif
#ifdef TESTBAND
  TestBandMatrix<float>();
  TestAllBandDiv<float>();
#endif
#ifdef TESTSYM
  TestSymMatrix<float>();
  TestAllSymDiv<float>();
#endif
#endif

#ifdef LONGDOUBLE
  TestAllVector<long double>();
  TestAllMatrix<long double>();
#ifdef TESTDIV
  TestAllMatrixDiv<long double>();
#endif
#ifdef TESTDIAG
  TestDiagMatrix<long double>();
  TestDiagDiv<long double>();
#endif
#ifdef TESTTRI
  TestTriMatrix<long double>();
  TestAllTriDiv<long double>();
#endif
#ifdef TESTBAND
  TestBandMatrix<long double>();
  TestAllBandDiv<long double>();
#endif
#ifdef TESTSYM
  TestSymMatrix<long double>();
  TestAllSymDiv<long double>();
#endif
#endif 

#endif // SKIPREST

  return 0;
}
