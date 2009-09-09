#include "TMV.h"
using tmv::Matrix;
using tmv::Vector;
using tmv::cerr;
using tmv::cout;
using tmv::endl;
#include "TMV_Diag.h"

#include "MemDebug.h"

#ifdef MEMDEBUG
AllocList* allocList=0;
#endif

#include <sys/time.h>

#define MATRIXDIV
//#define SING_MAT

const int M = 1500;
const int N = 1500;
//const int M = 6;
//const int N = 6;
#define S tmv::ColMajor

template <class T> void SpeedTest(const Matrix<T,S>& A,
    const Matrix<T,S>& A2)
{
  static Vector<T> v1(N,3.);
  static Vector<T> v2(N,3.);
  static Vector<T> v3(M);
  static Vector<T> v4(N);
  static Vector<T> v5(M);
  static Vector<T> v6(N);

#ifdef MATRIXDIV
  static Matrix<T,tmv::RowMajor> m1(N,N,3.);
  static Matrix<T,tmv::RowMajor> m2(N,N,3.);
  static Matrix<T> m3(M,N);
  static Matrix<T> m4(N,N);
  static Matrix<T> m5(N,M);
  static Matrix<T> m6(N,N);
#endif

  //cerr<<"A = "<<A<<endl;
  //cerr<<"A2 = "<<A2<<endl;
  timeval tp;

  gettimeofday(&tp,0);
  double t1 = tp.tv_sec + tp.tv_usec/1.e6;
  A.SetDiv();
  //cerr<<"After set\n";
  //A.CheckDecomp(A2,&cout);
  //cerr<<"After checkdecomp\n";
  gettimeofday(&tp,0);
  double t2 = tp.tv_sec + tp.tv_usec/1.e6;

  //cerr<<"A = "<<A2<<endl;
  //Matrix<T> U = A.SVFD().GetU();
  //cerr<<"U = "<<U<<endl;
  //tmv::DiagMatrix<T> D = DiagMatrixViewOf(A.SVFD().GetS());
  //cerr<<"D = "<<D<<endl;
  //Matrix<T> V = A.SVFD().GetV();
  //cerr<<"V = "<<V<<endl;
  //Matrix<T> UU = A2*V.Adjoint();
  //cerr<<"AVt = "<<UU<<endl;
  //Matrix<T> DD = D;
  //QR_Decompose(UU.View(),DD.View());
  //cerr<<"UU = "<<UU<<endl;
  //cerr<<"DD = "<<DD<<endl;
  //abort();
 
#ifdef SING_MAT
  v3 = A2*v1;
  v4 = v3/A;
  v1 = v4;
#endif
  //cout<<"v1 = "<<v1<<endl;
  v3 = A2*v1;
  //cout<<"v3 = A*v1 "<<v3<<endl;
  v4 = v3/A;
  //cout<<"v4 = v3/A "<<v4<<endl;
  //Vector<T> v4mv1 = v4-v1;
  //v4mv1.Clip(0.0001);
  //cout<<"v4-v1 = "<<v4mv1<<endl;
#ifdef SING_MAT
  v5 = v2%A;
  v6 = v5*A2;
  v2 = v6;
#endif
  //cout<<"v2 = "<<v2<<endl;
  v5 = v2%A;
  //cout<<"v5 = v2%A "<<v5<<endl;
  v6 = v5*A2;
  //cout<<"v6 = v5*A "<<v6<<endl;
  //Vector<T> v6mv2 = v6-v2;
  //v6mv2.Clip(0.0001);
  //cout<<"v6-v2 = "<<v6mv2<<endl;
  gettimeofday(&tp,0);
  double t3 = tp.tv_sec + tp.tv_usec/1.e6;

  cout<<"Norm(v4-v1) = "<<Norm(v4-v1)<<endl;
  cout<<"Norm(v6-v2) = "<<Norm(v6-v2)<<endl;

#ifdef MATRIXDIV
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

  cout<<"Norm(m4-m1) = "<<Norm(m4-m1)<<endl;
  cout<<"Norm(m6-m2) = "<<Norm(m6-m2)<<endl;
#endif

  cout<<"Time for Decomp = "<<t2-t1<<endl;
  cout<<"Time for Vector Division/Multiplication = "<<t3-t2<<endl;
#ifdef MATRIXDIV
  cout<<"Time for Matrix Multiplication = "<<t4-t3+t6-t5<<endl;
  cout<<"Time for Matrix Division = "<<t5-t4<<endl;
#endif
}

int main() {

#ifdef MEMDEBUG
  atexit(&DumpUnfreed);
#endif

  Matrix<double,S> A(M,N,0.8);

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

  cout<<"EPS*Norm(A) = "<<tmv::Epsilon<double>()*Norm(A)<<endl;

  /*
  if (M==N) {
    cout<<"\nLU Division: (Recursive) \n";
    Matrix<double,S> A2 = A;
    A2.DivideUsing(tmv::LU);
    A2.DivideInPlace();
    tmv::RecursiveLU = true;
    SpeedTest(A2,A);
  }

  if (M==N) {
    cout<<"\nLU Division: (Non-Recursive) \n";
    Matrix<double,S> A2 = A;
    A2.DivideUsing(tmv::LU);
    A2.DivideInPlace();
    tmv::RecursiveLU = false;
    SpeedTest(A2,A);
  }

  {
    cout<<"\nQR Division: (Recursive) \n";
    Matrix<double,S> A2 = A;
    A2.DivideUsing(tmv::QR);
    A2.DivideInPlace();
    tmv::RecursiveQR = true;
    SpeedTest(A2,A);
  }

  {
    cout<<"\nQR Division: (Non-Recursive) \n";
    Matrix<double,S> A2 = A;
    A2.DivideUsing(tmv::QR);
    A2.DivideInPlace();
    tmv::RecursiveQR = false;
    SpeedTest(A2,A);
  }

  {
    cout<<"\nQRP Division (Loose): \n";
    Matrix<double,S> A2 = A;
    A2.DivideUsing(tmv::QRP);
    A2.DivideInPlace();
    tmv::StrictQRP = false;
    SpeedTest(A2,A);
  }

  {
    cout<<"\nQRP Division (Strict): \n";
    Matrix<double,S> A2 = A;
    A2.DivideUsing(tmv::QRP);
    A2.DivideInPlace();
    tmv::StrictQRP = true;
    SpeedTest(A2,A);
  }
  */

  {
    cout<<"\nSVF Division: \n";
    Matrix<double,S> A2 = A;
    A2.DivideUsing(tmv::SVF);
    A2.DivideInPlace();
    SpeedTest(A2,A);
  }

  /*
  {
    cout<<"\nSV Division: \n";
    Matrix<double,S> A2 = A;
    A2.DivideUsing(tmv::SV);
    A2.DivideInPlace();
    SpeedTest(A2,A);
  }*/

  return 0;
}
