#include "TMV.h"
using tmv::Matrix;
using tmv::Vector;

#include "MemDebug.h"

#ifdef MEMDEBUG
AllocList* allocList=0;
#endif

#include <sys/time.h>

const int M = 500;
const int N = 500;
//const size_t M = 50;
//const size_t N = 50;
//const size_t M = 10;
//const size_t N = 2;

template <class T> void SpeedTest(const Matrix<T>& A)
{
  static Vector<T> v1(N,3.);
  static Vector<T> v2(N,3.);
  static Vector<T> v3(M);
  static Vector<T> v4(N);
  static Vector<T> v5(M);
  static Vector<T> v6(N);

  static Matrix<T> m1(N,M,3.,tmv::RowMajor);
  static Matrix<T> m2(M,N,3.,tmv::RowMajor);
  static Matrix<T> m3(M,M);
  static Matrix<T> m4(N,M);
  static Matrix<T> m5(M,M);
  static Matrix<T> m6(M,N);

  double t1,t2,t3,t4,t5,t6;
  timeval tp;

  gettimeofday(&tp,0);
  t1 = tp.tv_sec + tp.tv_usec/1.e6;
  A.SetDiv();
  gettimeofday(&tp,0);
  t2 = tp.tv_sec + tp.tv_usec/1.e6;
    
  v3 = A*v1;
  v4 = v3/A;
  v5 = v2%A;
  v6 = v5*A;
  gettimeofday(&tp,0);
  t3 = tp.tv_sec + tp.tv_usec/1.e6;

  m3 = A*m1;
  gettimeofday(&tp,0);
  t4 = tp.tv_sec + tp.tv_usec/1.e6;
  m4 = m3/A;
  m5 = m2%A;
  gettimeofday(&tp,0);
  t5 = tp.tv_sec + tp.tv_usec/1.e6;
  m6 = m5*A;
  gettimeofday(&tp,0);
  t6 = tp.tv_sec + tp.tv_usec/1.e6;

  std::cout<<"Norm(v4-v1) = "<<Norm(v4-v1)<<std::endl;
  std::cout<<"Norm(v6-v2) = "<<Norm(v6-v2)<<std::endl;
  std::cout<<"Norm(m4-m1) = "<<Norm(m4-m1)<<std::endl;
  std::cout<<"Norm(m6-m2) = "<<Norm(m6-m2)<<std::endl;
  std::cout<<"Time for Decomp = "<<t2-t1<<std::endl;
  std::cout<<"Time for Vector Division/Multiplication = "<<t3-t2<<std::endl;
  std::cout<<"Time for Matrix Multiplication = "<<t4-t3+t6-t5<<std::endl;
  std::cout<<"Time for Matrix Division = "<<t5-t4<<std::endl;
}

int main() {

#ifdef MEMDEBUG
  atexit(&DumpUnfreed);
#endif

  Matrix<double> A(M,N,tmv::ColMajor);

  double t1,t2;
  timeval tp;

  gettimeofday(&tp,0);
  t1 = tp.tv_sec + tp.tv_usec/1.e6;
  for(int i=0;i<M;i++) for(int j=0;j<N;j++) A(i,j) = 2*i-j;
  A.diag().AddToAll(100.*(M+N));
  gettimeofday(&tp,0);
  t2 = tp.tv_sec + tp.tv_usec/1.e6;

  std::cout<<"EPS*Norm(A) = "<<tmv::Epsilon<double>()*Norm(A)<<std::endl;

  std::cout<<"Time to make A = "<<t2-t1<<std::endl;

  std::cout<<"\nLU Division: \n";
  A.DivideUsing(tmv::LU);
  SpeedTest(A);

  std::cout<<"\nQR Division: \n";
  A.DivideUsing(tmv::QR);
  SpeedTest(A);

  std::cout<<"\nQRP Division: \n";
  A.DivideUsing(tmv::QRP);
  SpeedTest(A);

  std::cout<<"\nSV Division: \n";
  A.DivideUsing(tmv::SV);
  SpeedTest(A);

  return 0;
}
