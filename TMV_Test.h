#ifndef TMV_TEST_H
#define TMV_TEST_H

//#define LONGDOUBLE

#include <iostream>
#include <complex>
#include <vector>
#include <string>
#include "TMV.h"

using std::cerr;
using std::cout;
using std::endl;
using std::complex;
using std::abs;
using std::vector;
using std::ostream;
using std::string;

#define EPS (10*tmv::Epsilon<T>())

extern bool showtests;
extern bool showacc;
extern bool showdiv;
extern bool doallarith;
extern bool donorm2;
extern bool showstartdone;

inline void Assert(bool x,string s)
{
  if (!(x)) tmv::tmv_error("Error in test: ",s.c_str());
  else if (showtests) cout<<"Passed: "<<s<<endl;
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
