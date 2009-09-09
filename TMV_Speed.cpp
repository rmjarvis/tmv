#include "TMV.h"
#include "TMV_Diag.h"

#include "MemDebug.h"

#ifdef MEMDEBUG
AllocList* allocList=0;
#endif

#include <sys/time.h>

#define MATRIXDIV
//#define SING_MAT
#define DOINVERSE

const int M = 1500;
const int N = 1500;
//const int M = 6;
//const int N = 6;
#define S tmv::ColMajor

template <class T> void inline SpeedTest(const tmv::Matrix<T,S>& A,
    const tmv::Matrix<T,S>& A2)
{
  static tmv::Vector<T> v1(N,3.);
  static tmv::Vector<T> v2(N,3.);
  static tmv::Vector<T> v3(M);
  static tmv::Vector<T> v4(N);
  static tmv::Vector<T> v5(M);
  static tmv::Vector<T> v6(N);

#ifdef MATRIXDIV
  static tmv::Matrix<T,tmv::ColMajor> m1(N,N,3.);
  static tmv::Matrix<T,tmv::ColMajor> m2(N,N,3.);
  static tmv::Matrix<T,tmv::ColMajor> m3(M,N);
  static tmv::Matrix<T,tmv::ColMajor> m4(N,N);
  static tmv::Matrix<T,tmv::ColMajor> m5(N,M);
  static tmv::Matrix<T,tmv::ColMajor> m6(N,N);
#endif
#ifdef DOINVERSE
  static tmv::Matrix<T,tmv::ColMajor> Ainv(N,M);
#endif

  //std::cerr<<"A = "<<A<<std::endl;
  //std::cerr<<"A2 = "<<A2<<std::endl;
  timeval tp;

  gettimeofday(&tp,0);
  double t1 = tp.tv_sec + tp.tv_usec/1.e6;
  A.SetDiv();
  //A.CheckDecomp(A2,&std::cout);
  gettimeofday(&tp,0);
  double t2 = tp.tv_sec + tp.tv_usec/1.e6;

#ifdef SING_MAT
  v3 = A2*v1;
  v4 = v3/A;
  v1 = v4;
  v5 = v2%A;
  v6 = v5*A2;
  v2 = v6;
  m3 = A2*m1;
  m4 = m3/A;
  m1 = m4;
  m5 = m2%A;
  m6 = m5*A2;
  m2 = m6;
#endif
  gettimeofday(&tp,0);
  double t2b = tp.tv_sec + tp.tv_usec/1.e6;
  //std::cout<<"v1 = "<<v1<<std::endl;
  v3 = A2*v1;
  //std::cout<<"v3 = A*v1 "<<v3<<std::endl;
  v4 = v3/A;
  //std::cout<<"v4 = v3/A "<<v4<<std::endl;
  //tmv::Vector<T> v4mv1 = v4-v1;
  //v4mv1.Clip(0.0001);
  //std::cout<<"v4-v1 = "<<v4mv1<<std::endl;
  //std::cout<<"v2 = "<<v2<<std::endl;
  v5 = v2%A;
  //std::cout<<"v5 = v2%A "<<v5<<std::endl;
  v6 = v5*A2;
  //std::cout<<"v6 = v5*A "<<v6<<std::endl;
  //tmv::Vector<T> v6mv2 = v6-v2;
  //v6mv2.Clip(0.0001);
  //std::cout<<"v6-v2 = "<<v6mv2<<std::endl;
  gettimeofday(&tp,0);
  double t3 = tp.tv_sec + tp.tv_usec/1.e6;

  std::cout<<"Norm(v4-v1) = "<<Norm(v4-v1)<<std::endl;
  std::cout<<"Norm(v6-v2) = "<<Norm(v6-v2)<<std::endl;

#ifdef MATRIXDIV
  gettimeofday(&tp,0);
  double t3b = tp.tv_sec + tp.tv_usec/1.e6;
  m3 = A2*m1;
  gettimeofday(&tp,0);
  double t4 = tp.tv_sec + tp.tv_usec/1.e6;
  m4 = m3/A;
  m5 = m2%A;
  gettimeofday(&tp,0);
  double t5 = tp.tv_sec + tp.tv_usec/1.e6;
  m6 = m5*A2;
  gettimeofday(&tp,0);
  double t6 = tp.tv_sec + tp.tv_usec/1.e6;

  std::cout<<"Norm(m4-m1) = "<<Norm(m4-m1)<<std::endl;
  std::cout<<"Norm(m6-m2) = "<<Norm(m6-m2)<<std::endl;
#endif
#ifdef DOINVERSE
  gettimeofday(&tp,0);
  double t6b = tp.tv_sec + tp.tv_usec/1.e6;
  Ainv = A.Inverse();
  gettimeofday(&tp,0);
  double t7 = tp.tv_sec + tp.tv_usec/1.e6;
  std::cout<<"Norm(Ainv*A-1) = "<<Norm(Ainv*A2-T(1))<<std::endl;
#endif

  std::cout<<"Time for Decomp = "<<t2-t1<<std::endl;
  std::cout<<"Time for Vector Division/Multiplication = "<<t3-t2b<<std::endl;
#ifdef MATRIXDIV
  std::cout<<"Time for Matrix Multiplication = "<<t4-t3b+t6-t5<<std::endl;
  std::cout<<"Time for Matrix Division = "<<t5-t4<<std::endl;
#endif
#ifdef DOINVERSE
  std::cout<<"Time for Inverse = "<<t7-t6b<<std::endl;
#endif
}

int main() {

#ifdef MEMDEBUG
  atexit(&DumpUnfreed);
#endif

  tmv::Matrix<double,S> A(M,N,0.8);

  //for(int i=0;i<M;i++) for(int j=0;j<N;j++) A(i,j) = 2*i-j;
  A.diag().AddToAll(100.*N);
  A.diag(-N/5).AddToAll(200.*N);
  A.row(3*N/5).AddToAll(-300.*N);
  A.SubVector(4*N/5,N/3,-1,2,N/3).AddToAll(400.*N);
#ifdef SING_MAT
  A.col(0).Zero();
  A.col(N/5) *= 1.e-8;
  A.col(2*N/5) *= 1.e-18;
  A.col(3*N/5) *= 1.e-15;
  A.col(3*N/5) += A.col(4*N/5);
#endif

  std::cout<<"EPS*Norm(A) = "<<tmv::Epsilon<double>()*Norm(A)<<std::endl;

#ifndef SING_MAT
  if (M==N) {
    std::cout<<"\nLU Division: \n";
    tmv::Matrix<double,S> A2 = A;
    A2.DivideUsing(tmv::LU);
    A2.DivideInPlace();
    A2.SaveDiv();
    SpeedTest(A2,A);
  }

  {
    std::cout<<"\nQR Division: \n";
    tmv::Matrix<double,S> A2 = A;
    A2.DivideUsing(tmv::QR);
    A2.DivideInPlace();
    A2.SaveDiv();
    SpeedTest(A2,A);
  }
#endif

  {
    std::cout<<"\nQRP Division (Loose): \n";
    tmv::Matrix<double,S> A2 = A;
    A2.DivideUsing(tmv::QRP);
    A2.DivideInPlace();
    A2.SaveDiv();
    tmv::StrictQRP = false;
    SpeedTest(A2,A);
  }

  {
    std::cout<<"\nQRP Division (Strict): \n";
    tmv::Matrix<double,S> A2 = A;
    A2.DivideUsing(tmv::QRP);
    A2.DivideInPlace();
    A2.SaveDiv();
    tmv::StrictQRP = true;
    SpeedTest(A2,A);
  }

  {
    std::cout<<"\nSV Division: \n";
    tmv::Matrix<double,S> A2 = A;
    A2.DivideUsing(tmv::SV);
    A2.DivideInPlace();
    A2.SaveDiv();
    SpeedTest(A2,A);
  }

  return 0;
}
