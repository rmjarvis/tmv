//#define SHOWTESTS
#define TESTDIV

#include "TMV.h"
#include "TMV_Diag.h"
#include "TMV_Test.h"
#include "TMV_TestMatrixDiv.h"

using tmv::Matrix;
using tmv::DiagMatrix;
using tmv::RowMajor;

template <class T> void TestDiagDiv()
{
  const int N = 10;

  DiagMatrix<T> a(N);
  DiagMatrix<T> b(N);
  for (int i=0; i<N; ++i) {
    a(i,i) = 3.+5*i;
    b(i,i) = 5.-2*i;
  }

  DiagMatrix<complex<T> > ca = a*complex<T>(1,2);
  DiagMatrix<complex<T> > cb = b*complex<T>(-5,-1);

  Matrix<T> p(N,N,5);
  p.diag().AddToAll(N*10);
  Matrix<complex<T> > cp = p*complex<T>(2,3);
  cp += ca;

  Matrix<T> q(2*N,N,-2);
  q.Rows(0,N) += p;
  q.Rows(N,2*N) -= p;
  q.Rows(N/2,3*N/2) += T(4)*p;
  Matrix<complex<T> > cq = q*complex<T>(-1,4);
  cq.Rows(0,N) -= cp;
  cq.Rows(N,2*N) -= cp;
  cq.Rows(N/2,3*N/2) -= cp;

  Matrix<T> r(N,0,1);
  Matrix<complex<T> > cr(N,0,1);

  TestMatrixDivArith<T>(tmv::LU,a,b,ca,cb,"Diag/Diag");
  TestMatrixDivArith<T>(tmv::LU,a,p,ca,cp,"SquareM/Diag");
  TestMatrixDivArith<T>(tmv::LU,a,q,ca,cq,"NonSqaureM/Diag");
  TestMatrixDivArith<T>(tmv::LU,a,r,ca,cr,"DegenM/Diag");
  TestMatrixDivArith<T>(tmv::LU,p,b,cp,cb,"Diag/SquareM");

  cout<<"DiagMatrix<"<<tmv::Type(T())<<"> Division passed all tests\n";
}

template void TestDiagDiv<double>();
#ifndef NOFLOAT
template void TestDiagDiv<float>();
#endif
#ifdef LONGDOUBLE
template void TestDiagDiv<long double>();
#endif
