#include "TMV.h"
#include "TMV_Sym.h"

#include "MemDebug.h"

#ifdef MEMDEBUG
AllocList* allocList=0;
#endif

#include <mkl.h>

#include <sys/time.h>

#define MATRIXDIV
//#define SING_MAT

const int N = 2000;
const int NRHS = 20;
//const int N = 7;
#define U tmv::Lower
#define S tmv::ColMajor

template <class T> inline void SpeedTest(const tmv::SymMatrix<T,U,S>& A,
    const tmv::SymMatrix<T,U,S>& A2)
{
  static tmv::Vector<T> v1(N,3.);
  static tmv::Vector<T> v2(N,3.);
  static tmv::Vector<T> v3(N);
  static tmv::Vector<T> v4(N);
  static tmv::Vector<T> v5(N);
  static tmv::Vector<T> v6(N);

#ifdef MATRIXDIV
  static tmv::Matrix<T,tmv::ColMajor> m1(N,NRHS,3.);
  static tmv::Matrix<T,tmv::ColMajor> m2(N,NRHS,3.);
  static tmv::Matrix<T,tmv::ColMajor> m3(N,NRHS);
  static tmv::Matrix<T,tmv::ColMajor> m4(N,NRHS);
  static tmv::Matrix<T,tmv::ColMajor> m5(N,NRHS);
  static tmv::Matrix<T,tmv::ColMajor> m6(N,NRHS);
#endif

  //std::cout<<"A = "<<A<<std::endl;

  tmv::auto_array<T> lapA(new T[A.size()*A.size()]);
  memmove(lapA.get(),A.cptr(),A.size()*A.size()*sizeof(T));
  tmv::auto_array<T> lapAinv(new T[A.size()*A.size()]);
  tmv::SymMatrix<T,U,S> Ainv(A.size());

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
#endif

  //std::cout<<"v1 = "<<v1<<std::endl;
  v3 = A2*v1;
  gettimeofday(&tp,0);
  double t2b = tp.tv_sec + tp.tv_usec/1.e6;
  //std::cout<<"v3 = A*v1 "<<v3<<std::endl;
  v4 = v3/A;
  //std::cout<<"v4 = v3/A "<<v4<<std::endl;
  //tmv::Vector<T> v4mv1 = v4-v1;
  //v4mv1.Clip(0.0001);
  //std::cout<<"v4-v1 = "<<v4mv1<<std::endl;
  //std::cout<<"v2 = "<<v2<<std::endl;
  //v5 = v2%A;
  //std::cout<<"v5 = v2%A "<<v5<<std::endl;
  //v6 = v5*A2;
  //std::cout<<"v6 = v5*A "<<v6<<std::endl;
  //tmv::Vector<T> v6mv2 = v6-v2;
  //v6mv2.Clip(0.0001);
  //std::cout<<"v6-v2 = "<<v6mv2<<std::endl;
  gettimeofday(&tp,0);
  double t3 = tp.tv_sec + tp.tv_usec/1.e6;

  std::cout<<"Norm(v4-v1) = "<<Norm(v4-v1)<<std::endl;
  //std::cout<<"Norm(v6-v2) = "<<Norm(v6-v2)<<std::endl;

#ifdef MATRIXDIV
  m3 = A2*m1;
  gettimeofday(&tp,0);
  double t4 = tp.tv_sec + tp.tv_usec/1.e6;
  m4 = m3/A;
  //m5 = m2%A;
  gettimeofday(&tp,0);
  double t5 = tp.tv_sec + tp.tv_usec/1.e6;
  //m6 = m5*A2;
  gettimeofday(&tp,0);
  double t6 = tp.tv_sec + tp.tv_usec/1.e6;

  std::cout<<"Norm(m4-m1) = "<<Norm(m4-m1)<<std::endl;
  //std::cout<<"Norm(m6-m2) = "<<Norm(m6-m2)<<std::endl;

  gettimeofday(&tp,0);
  double t6b = tp.tv_sec + tp.tv_usec/1.e6;
  Ainv = A.Inverse();
  gettimeofday(&tp,0);
  double t7 = tp.tv_sec + tp.tv_usec/1.e6;

  //std::cerr<<"A = "<<A2<<std::endl;
  //std::cerr<<"Ainv = "<<Ainv<<std::endl;
  //std::cerr<<"Ainv*A = "<<Ainv*A2<<std::endl;
  //std::cerr<<"Ainv*A-1 = "<<(Ainv*A2-T(1))<<std::endl;
  tmv::Matrix<T> minv = A.Inverse();
  //std::cerr<<"minv*A = "<<minv*A2<<std::endl;
  //std::cerr<<"minv*A-1 = "<<(minv*A2-T(1))<<std::endl;

  std::cout<<"Norm(A^-1 A - 1) = "<<Norm(Ainv*A2-T(1))<<std::endl;;
  std::cout<<"Norm(A^-1 A - 1) = "<<Norm(minv*A2-T(1))<<std::endl;;
#endif

  std::cout<<"Time for Decomp = "<<t2-t1<<std::endl;
  std::cout<<"Time for Vector Division = "<<t3-t2b<<std::endl;
#ifdef MATRIXDIV
  std::cout<<"Time for Matrix Multiplication = "<<t4-t3+t6-t5<<std::endl;
  std::cout<<"Time for Matrix Division = "<<t5-t4<<std::endl;
  std::cout<<"Time for Matrix Inverse = "<<t7-t6b<<std::endl;
#endif

  std::cout<<"Lap version:\n";
  gettimeofday(&tp,0);
  t1 = tp.tv_sec + tp.tv_usec/1.e6;
  int n = A.size();
  int lda = A.stepj();
  int lap_p[A.size()];
#ifndef LAPNOWORK
  int lwork = n*64;
  double work[lwork];
#endif
  LAPNAME(dsytrf) (LAPCM LAPCH_LO,LAPV(n),LAPP(lapA.get()),LAPV(lda),
      LAPP(lap_p) LAPWK(work) LAPVWK(lwork) LAPINFO LAP1);
  if (tmv::Lap_info) std::cout<<"dsytrf info = "<<tmv::Lap_info<<std::endl;
  gettimeofday(&tp,0);
  t2 = tp.tv_sec + tp.tv_usec/1.e6;

  //std::cout<<"v1 = "<<v1<<std::endl;
  v3 = A2*v1;
  gettimeofday(&tp,0);
  t2b = tp.tv_sec + tp.tv_usec/1.e6;
  //std::cout<<"v3 = A*v1 "<<v3<<std::endl;
  v4 = v3;
  int nrhs = 1;
  int ldb = n;
  LAPNAME(dsytrs) (LAPCM LAPCH_LO,LAPV(n),LAPV(nrhs),LAPP(lapA.get()),LAPV(lda),
      LAPP(lap_p),LAPP(v4.ptr()),LAPV(ldb) LAPINFO LAP1);
  if (tmv::Lap_info) std::cout<<"dsytrs info = "<<tmv::Lap_info<<std::endl;
  //std::cout<<"v4 = v3/A "<<v4<<std::endl;
  //std::cout<<"v2 = "<<v2<<std::endl;
  //v5 = v2;
  //LAPNAME(dsytrs) (LAPCM LAPCH_LO,LAPV(n),LAPV(nrhs),LAPP(lapA.get()),LAPV(lda),
  //    LAPP(lap_p),LAPP(v5.ptr()),LAPV(ldb) LAPINFO LAP1);
  //if (tmv::Lap_info) std::cout<<"dsytrs info = "<<tmv::Lap_info<<std::endl;
  //std::cout<<"v5 = v2%A "<<v5<<std::endl;
  //v6 = v5*A2;
  //std::cout<<"v6 = v5*A "<<v6<<std::endl;
  gettimeofday(&tp,0);
  t3 = tp.tv_sec + tp.tv_usec/1.e6;

  std::cout<<"Norm(v4-v1) = "<<Norm(v4-v1)<<std::endl;
  //std::cout<<"Norm(v6-v2) = "<<Norm(v6-v2)<<std::endl;

#ifdef MATRIXDIV
  m3 = A2*m1;
  gettimeofday(&tp,0);
  t4 = tp.tv_sec + tp.tv_usec/1.e6;
  m4 = m3;
  nrhs = m4.rowsize();
  ldb = m4.stepj();
  LAPNAME(dsytrs) (LAPCM LAPCH_LO,LAPV(n),LAPV(nrhs),LAPP(lapA.get()),LAPV(lda),
      LAPP(lap_p),LAPP(m4.ptr()),LAPV(ldb) LAPINFO LAP1);
  if (tmv::Lap_info) std::cout<<"dsytrs info = "<<tmv::Lap_info<<std::endl;
  //m5.Transpose() = m2;
  //ldb = m5.stepj();
  //dsytrs(&u,&n,&nrhs,lapA.get(),&lda,lap_p,m5.ptr(),&ldb,&info);
  gettimeofday(&tp,0);
  t5 = tp.tv_sec + tp.tv_usec/1.e6;
  //m6 = m5.Transpose()*A2;
  gettimeofday(&tp,0);
  t6 = tp.tv_sec + tp.tv_usec/1.e6;

  std::cout<<"Norm(m4-m1) = "<<Norm(m4-m1)<<std::endl;
  //std::cout<<"Norm(m6-m2) = "<<Norm(m6-m2)<<std::endl;

  gettimeofday(&tp,0);
  t6b = tp.tv_sec + tp.tv_usec/1.e6;
  memmove(lapAinv.get(),lapA.get(),A.size()*A.size()*sizeof(T));
  LAPNAME(dsytri) (LAPCM LAPCH_LO,LAPV(n),LAPP(lapAinv.get()),LAPV(lda),
      LAPP(lap_p) LAPWK(work) LAPINFO LAP1);
  if (tmv::Lap_info) std::cout<<"dsytri info = "<<tmv::Lap_info<<std::endl;
  gettimeofday(&tp,0);
  t7 = tp.tv_sec + tp.tv_usec/1.e6;
  
  std::cout<<"Norm(A^-1 A - 1) = "<<Norm(
      tmv::SymMatrixViewOf(lapAinv.get(),A.size(),U,S)*A2 - T(1))<<std::endl;;
#endif

  std::cout<<"Time for Decomp = "<<t2-t1<<std::endl;
  std::cout<<"Time for Vector Division = "<<t3-t2b<<std::endl;
#ifdef MATRIXDIV
  std::cout<<"Time for Matrix Multiplication = "<<t4-t3+t6-t5<<std::endl;
  std::cout<<"Time for Matrix Division = "<<t5-t4<<std::endl;
  std::cout<<"Time for Matrix Inverse = "<<t7-t6b<<std::endl;
#endif

}

int main() {

#ifdef MEMDEBUG
  atexit(&DumpUnfreed);
#endif

  tmv::SymMatrix<double,U,S> A(N,0.8);

  A.diag().AddToAll(100.*N);
  A.diag(-N/5).AddToAll(200.*N);
  A.row(3*N/5,0,3*N/5).AddToAll(-300.*N);
  A.SubVector(4*N/5,N/8,-1,2,N/6).AddToAll(400.*N);
#ifdef SING_MAT
  A.col(0,0,N).Zero();
  A.col(N/5,N/5,N) *= 1.e-8;
  A.row(N/5,0,N/5) *= 1.e-8;
  A.col(2*N/5,2*N/5,N) *= 1.e-18;
  A.row(2*N/5,0,2*N/5) *= 1.e-18;
  A.col(3*N/5,3*N/5,N) *= 1.e-15;
  A.row(3*N/5,0,3*N/5) *= 1.e-15;
  A.col(3*N/5,4*N/5,N) += A.col(4*N/5,4*N/5,N);
  A.col(3*N/5,3*N/5,4*N/5) += A.row(4*N/5,3*N/5,4*N/5);
  A.row(3*N/5,0,3*N/5) += A.row(4*N/5,0,3*N/5);
#endif

  std::cout<<"EPS*Norm(A) = "<<tmv::Epsilon<double>()*Norm(A)<<std::endl;

  {
    std::cout<<"\nLU Division: \n";
    tmv::SymMatrix<double,U,S> A2 = A;
    A2.DivideUsing(tmv::LU);
    A2.DivideInPlace();
    A2.SaveDiv();
    SpeedTest(A2,A);
  }

  return 0;
}
