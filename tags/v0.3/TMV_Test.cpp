
//#define SHOWCHECK
//#define SHOWACC
//#define SHOWTESTS
#define TESTDIV
#define TESTDIAG
#define TESTTRI
#define TESTBAND

#include "TMV.h"
#include "TMV_Test.h"
using tmv::Matrix;

#ifdef TESTDIAG
#include "TMV_Diag.h"
#endif

#ifdef TESTTRI
#include "TMV_Tri.h"
using tmv::UpperTriMatrix;
#endif

#include "MemDebug.h"
#ifdef TESTBAND
#include "TMV_Band.h"
using tmv::BandMatrix;
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

int main() {

#ifdef MEMDEBUG
  atexit(&DumpUnfreed);
#endif

  TestAllVector<double>();
  TestAllMatrix<double>();
  TestPermutation();
#ifdef TESTTRI
  TestTriMatrix<double>();
  TestAllTriDiv<double>();
#endif
#ifdef TESTDIAG
  TestDiagMatrix<double>();
  TestDiagDiv<double>();
#endif
#ifdef TESTDIV
  TestAllMatrixDiv<double>();
#endif
#ifdef TESTBAND
  TestBandMatrix<double>();
  TestAllBandDiv<double>();
#endif

#ifndef NOFLOAT
  TestAllVector<float>();
  TestAllMatrix<float>();
#ifdef TESTDIAG
  TestDiagMatrix<float>();
  TestDiagDiv<float>();
#endif
#ifdef TESTTRI
  TestTriMatrix<float>();
  TestAllTriDiv<float>();
#endif
#ifdef TESTDIV
  TestAllMatrixDiv<float>();
#endif
#ifdef TESTBAND
  TestBandMatrix<float>();
  TestAllBandDiv<float>();
#endif
#endif

#ifdef LONGDOUBLE
  TestAllVector<long double>();
  TestAllMatrix<long double>();
#ifdef TESTDIAG
  TestDiagMatrix<long double>();
  TestDiagDiv<long double>();
#endif
#ifdef TESTTRI
  TestTriMatrix<long double>();
  TestAllTriDiv<long double>();
#endif
#ifdef TESTDIV
  TestAllMatrixDiv<long double>();
#endif
#ifdef TESTBAND
  TestBandMatrix<long double>();
  TestAllBandDiv<long double>();
#endif
#endif 

  return 0;
}
