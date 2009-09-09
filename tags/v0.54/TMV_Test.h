#ifndef TMV_TEST_H
#define TMV_TEST_H

//#define LONGDOUBLE

#define EPS (10*tmv::Epsilon<T>())

#include "TMV.h"

extern bool showtests;
extern bool showacc;
extern bool showdiv;
extern bool doallarith;
extern bool donorm2;
extern bool showstartdone;

inline void Assert(bool x,std::string s)
{
  if (!(x)) throw tmv::Error(std::string("Error in test: ") + s);
  else if (showtests) std::cout<<"Passed: "<<s<<std::endl;
}

template <class T> extern void TestAllVector();
template <class T> extern void TestAllMatrix();
template <class T> extern void TestAllMatrixDiv();
template <class T> extern void TestDiagMatrix();
template <class T> extern void TestDiagDiv();
template <class T> extern void TestTriMatrix();
template <class T> extern void TestAllTriDiv();
template <class T> extern void TestBandMatrix();
template <class T> extern void TestAllBandDiv();
template <class T> extern void TestSymMatrix();
template <class T> extern void TestAllSymDiv();

template <class T> extern void TestTriMatrixArith_A();
template <class T> extern void TestTriMatrixArith_B();
template <class T> extern void TestBandMatrixArith_A();
template <class T> extern void TestBandMatrixArith_B();
template <class T> extern void TestSymMatrixArith_A();
template <class T> extern void TestSymMatrixArith_B();


#endif
