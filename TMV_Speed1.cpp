#include "TMV.h"
using tmv::Matrix;
using tmv::Vector;
using tmv::VectorView;
using tmv::cerr;
using tmv::cout;
using tmv::endl;

#include "MemDebug.h"

#ifdef MEMDEBUG
AllocList* allocList=0;
#endif

#include <sys/time.h>

// Use N >= M
const int M = 50;
const int N = 500000;
#define S tmv::RowMajor

template <class T> T Prod(const tmv::GenVector<T>& v1,
    const tmv::GenVector<T>& v2)
{
  if (v1.size() == 1) return v1(0) * v2(0);
  else {
    size_t N = v1.size();
    size_t N1 = N/2;
    return Prod(v1.SubVector(0,N1),v2.SubVector(0,N1)) +
      Prod(v1.SubVector(N1,N),v2.SubVector(N1,N));
  }
}

template <class T> void SpeedTest(const Matrix<T,S>& A)
{
  static Vector<T> v(N,3.);
  static Vector<T> w(M);
  static Vector<T> w2(M);
  static Vector<T> u(M);
  static Matrix<T,tmv::ColMajor> m1(N,2,3.);
  VectorView<T> v1 = m1.col(0);
  VectorView<T> u1 = m1.col(1,0,M);
  static Matrix<T,tmv::RowMajor> m2(N,2,3.);
  VectorView<T> v2 = m2.col(0);
  VectorView<T> u2 = m2.col(1,0,M);
  static Matrix<T,tmv::RowMajor> m16(N,16,3.);
  VectorView<T> v16 = m16.col(0);
  VectorView<T> u16 = m16.col(1,0,M);

  timeval tp;

  for(size_t i=0;i<M;i++) {
    double temp = 0.;
    for(size_t j=0;j<N;j++) temp += A(i,j)*v(j);
    w2(i) = temp;
    w(i) = Prod(A.row(i),v);
  }
  cerr<<"Norm(w-w2) = "<<Norm(w-w2)<<endl;

  gettimeofday(&tp,0);
  double t1 = tp.tv_sec + tp.tv_usec/1.e6;

  u = A * v;
  u1 = A * v1;

  gettimeofday(&tp,0);
  double t2 = tp.tv_sec + tp.tv_usec/1.e6;

  u2 = A * v2;
  u16 = A * v16;

  gettimeofday(&tp,0);
  double t3 = tp.tv_sec + tp.tv_usec/1.e6;

  cout<<"Norm(v) = "<<Norm(v)<<endl;
  cout<<"Norm(A) = "<<Norm(A)<<endl;
  cout<<"Norm(w) = "<<Norm(w)<<endl;
  cout<<"Norm(u-w) = "<<Norm(u-w)<<endl;
  cout<<"Norm(u1-w) = "<<Norm(u1-w)<<endl;
  cout<<"Norm(u2-w) = "<<Norm(u2-w)<<endl;
  cout<<"Norm(u16-w) = "<<Norm(u16-w)<<endl;

  cout<<"Norm(u1-w)/Norm(v)/Norm(A) = "<<Norm(u1-w)/Norm(v)/Norm(A)/tmv::Epsilon<double>()<<" * epsilon\n";
  cout<<"Norm(u16-w)/Norm(v)/Norm(A) = "<<Norm(u16-w)/Norm(v)/Norm(A)/tmv::Epsilon<double>()<<" * epsilon\n";

  cout<<"Time for Multiplication UU = "<<t2-t1<<endl;
  cout<<"Time for Multiplication SS = "<<t3-t2<<endl;
}

int main() {

#ifdef MEMDEBUG
  atexit(&DumpUnfreed);
#endif

  Matrix<double,S> A(M,N,0.8);

  //for(int i=0;i<M;i++) for(int j=0;j<N;j++) A(i,j) = 2*i-j;
  A.diag().AddToAll(100.*M);
  A.diag(-M/5).AddToAll(200.*M);
  A.row(3*M/5).AddToAll(-300.*M);
  A.SubVector(4*M/5,M/3,-1,2,M/3).AddToAll(400.*M);

  cout<<"EPS*Norm(A) = "<<tmv::Epsilon<double>()*Norm(A)<<endl;

  SpeedTest(A);

  return 0;
}
