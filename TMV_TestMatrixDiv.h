#include "TMV.h"
using tmv::Matrix;
using tmv::Vector;
using tmv::StorageType;
using tmv::RowMajor;
using tmv::ColMajor;
using tmv::Type;
using tmv::DivType;

template <class SM1, class SM2> bool CanLDiv(const SM1& a, const SM2& b)
{ return a.colsize() == b.colsize(); }

template <class T, class SM2> bool CanLDiv(const Vector<T>& a, const SM2& b)
{ return a.size() == b.colsize(); }

template <class SM1, class SM2> bool CanLDivEq(const SM1& a, const SM2& b)
{ return CanLDiv(a,b) && b.IsSquare(); }

template <class SM1, class SM2> bool CanRDiv(const SM1& a, const SM2& b)
{ return a.rowsize() == b.rowsize(); }

template <class T, class SM2> bool CanRDiv(const Vector<T>& a, const SM2& b)
{ return a.size() == b.rowsize(); }

template <class SM1, class SM2> bool CanRDivEq(const SM1& a, const SM2& b)
{ return CanRDiv(a,b) && b.IsSquare(); }

template <class SM1, class SM2, class T, StorageType S, class T2, StorageType S2> void DoTestMatrixDivArithMM1b(
    const SM1& a, const SM2& b, const Matrix<T,S>& m1,
    const Matrix<T2,S2>& m2, string label)
{
  if (showstartdone) {
    cout<<"Start MM1b: "<<label<<endl;
    cout<<"b.dt = "<<b.GetDivType()<<endl;
    cout<<"m2.dt = "<<m2.GetDivType()<<endl;
  }
  b.SaveDiv();
  m2.SaveDiv();
  b.SetDiv();
  m2.SetDiv();

//#define XDEBUG

#ifdef XDEBUG
  cout <<"b.dt = "<<Text(b.GetDivType())<<" = "<<b.GetDiv()->Type()<<endl;
  cout <<"m.dt = "<<Text(m2.GetDivType())<<" = "<<m2.GetDiv()->Type()<<endl;
  cout<<"a = "<<Type(a)<<"  "<<a<<endl;
  cout<<"b = "<<Type(b)<<"  "<<b<<endl;
  cout<<"m1 = "<<Type(m1)<<"  "<<m1<<endl;
  cout<<"m2 = "<<Type(m2)<<"  "<<m2<<endl;
#endif
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a!=m1");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b!=m2");

  Matrix<T2,ColMajor> m2inv = Inverse(m2);
  double eps = EPS*Norm(m1)*std::max(m2.colsize(),m2.rowsize())*Norm(m2)*Norm(m2inv);

  if (CanLDiv(a,b)) {
#ifdef XDEBUG
    cout<<"LDiv:\n";
    cout<<"a = "<<Type(a)<<a<<endl;
    cout<<"m1 = "<<Type(m1)<<m1<<endl;
    cout<<"b = "<<Type(b)<<b<<endl;
    cout<<"m2 = "<<Type(m2)<<m2<<endl;
    cout<<"m1/m2 = "<<m1/m2<<endl;
    cout<<"a/m2 = "<<a/m2<<endl;
    cout<<"m1/b = "<<m1/b<<endl;
    cout<<"a/b = "<<a/b<<endl;
    cout<<"Norm(diff) = "<<Norm(a/b-m1/m2)<<" <?  "<<eps<<endl;
    cout<<"Norm(diff) = "<<Norm(a/m2-m1/m2)<<" <?  "<<eps<<endl;
    cout<<"Norm(diff) = "<<Norm(m1/b-m1/m2)<<" <?  "<<eps<<endl;
#endif
    Assert(Norm((a/b)-(m1/m2)) <= eps,label+" a/b");
    Assert(Norm((a/m2)-(m1/m2)) <= eps,label+" a/m");
    Assert(Norm((m1/b)-(m1/m2)) <= eps,label+" m/b");
    Assert(Norm((b.Inverse()*a)-(m1/m2)) <= eps,label+" b^-1*a");
    Assert(Norm((m2.Inverse()*a)-(m1/m2)) <= eps,label+" m^-1*a");
    Assert(Norm((b.Inverse()*m1)-(m1/m2)) <= eps,label+" b^-1*m");
  }
  if (CanRDiv(a,b)) {
#ifdef XDEBUG
    cout<<"RDiv:\n";
    cout<<"a = "<<Type(a)<<a<<endl;
    cout<<"m1 = "<<Type(m1)<<m1<<endl;
    cout<<"b = "<<Type(b)<<b<<endl;
    cout<<"m2 = "<<Type(m2)<<m2<<endl;
    cout<<"m1%m2 = "<<m1%m2<<endl;
    cout<<"a%m2 = "<<a%m2<<endl;
    cout<<"m1%b = "<<m1%b<<endl;
    cout<<"a%b = "<<a%b<<endl;
    cout<<"Norm(diff) = "<<Norm(a%b-m1%m2)<<" <?  "<<eps<<endl;
    cout<<"Norm(diff) = "<<Norm(a%m2-m1%m2)<<" <?  "<<eps<<endl;
    cout<<"Norm(diff) = "<<Norm(m1%b-m1%m2)<<" <?  "<<eps<<endl;
#endif
    Assert(Norm((a%b)-(m1%m2)) <= eps,label+" a%b");
    Assert(Norm((a%m2)-(m1%m2)) <= eps,label+" a%m");
    Assert(Norm((m1%b)-(m1%m2)) <= eps,label+" m%b");
    Assert(Norm((a*b.Inverse())-(m1%m2)) <= eps,label+" a*b^-1");
    Assert(Norm((a*m2.Inverse())-(m1%m2)) <= eps,label+" a*m^-1");
    Assert(Norm((m1*b.Inverse())-(m1%m2)) <= eps,label+" m*b^-1");
  }

#ifdef XDEBUG
#undef XDEBUG
#endif
  if (showstartdone) 
    cout<<"Done MM1b"<<endl;
}

template <class SM1, class SM2, class T, StorageType S, class T2, StorageType S2> void DoTestMatrixDivArithMM1a(
    DivType dt, const SM1& a, const SM2& b, const Matrix<T,S>& m1,
    const Matrix<T2,S2>& m2, string label)
{
  if (showstartdone) {
    cout<<"Start MM1a: "<<label<<endl;
    cout<<"b.dt = "<<b.GetDivType()<<endl;
    cout<<"m2.dt = "<<m2.GetDivType()<<endl;
  }
  b.DivideUsing(dt);
  DoTestMatrixDivArithMM1b(a,b,m1,m2,label);
  if (showstartdone) 
    cout<<"Done MM1a"<<endl;
}

template <class SM1, class SM2, class T, class T2> void DoTestMatrixDivArithMM1(
    DivType dt, const SM1& a, const SM2& b, const Matrix<T,RowMajor>& m1,
    const Matrix<T2,RowMajor>& m2, string label)
{
  if (showstartdone) {
    cout<<"Start MM1: "<<label<<endl;
    cout<<"b.dt = "<<b.GetDivType()<<endl;
    cout<<"m2.dt = "<<m2.GetDivType()<<endl;
  }
  b.DivideUsing(dt);
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a!=m1");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b!=m2");

  Matrix<T,ColMajor> m1t = m1.Transpose();
  Matrix<T,RowMajor> m1c = m1.Conjugate();
  Matrix<T,ColMajor> m1a = m1.Adjoint();
  Matrix<T2,ColMajor> m2t = m2.Transpose();
  Matrix<T2,RowMajor> m2c = m2.Conjugate();
  Matrix<T2,ColMajor> m2a = m2.Adjoint();

  DoTestMatrixDivArithMM1a(dt,a.View(),b.View(),m1,m2,label);
  DoTestMatrixDivArithMM1a(dt,Transpose(a),b.View(),m1t,m2,label+" TransA");
  DoTestMatrixDivArithMM1a(dt,Conjugate(a),b.View(),m1c,m2,label+" ConjA");
  DoTestMatrixDivArithMM1a(dt,Adjoint(a),b.View(),m1a,m2,label+" AdjA");

  DoTestMatrixDivArithMM1a(dt,a.View(),Transpose(b),m1,m2t,label+" TransB");
  DoTestMatrixDivArithMM1a(dt,Transpose(a),Transpose(b),m1t,m2t,label+" TransA TransB");
  DoTestMatrixDivArithMM1a(dt,Conjugate(a),Transpose(b),m1c,m2t,label+" ConjA TransB");
  DoTestMatrixDivArithMM1a(dt,Adjoint(a),Transpose(b),m1a,m2t,label+" AdjA TransB");

  DoTestMatrixDivArithMM1a(dt,a.View(),Conjugate(b),m1,m2c,label+" ConjB");
  DoTestMatrixDivArithMM1a(dt,Transpose(a),Conjugate(b),m1t,m2c,label+" TransA ConjB");
  DoTestMatrixDivArithMM1a(dt,Conjugate(a),Conjugate(b),m1c,m2c,label+" ConjA ConjB");
  DoTestMatrixDivArithMM1a(dt,Adjoint(a),Conjugate(b),m1a,m2c,label+" AdjA ConjB");

  DoTestMatrixDivArithMM1a(dt,a.View(),Adjoint(b),m1,m2a,label+" AdjB");
  DoTestMatrixDivArithMM1a(dt,Transpose(a),Adjoint(b),m1t,m2a,label+" TransA AdjB");
  DoTestMatrixDivArithMM1a(dt,Conjugate(a),Adjoint(b),m1c,m2a,label+" ConjA AdjB");
  DoTestMatrixDivArithMM1a(dt,Adjoint(a),Adjoint(b),m1a,m2a,label+" AdjA AdjB");
  if (showstartdone) 
    cout<<"Done MM1"<<endl;
}

template <class SM1, class SM2, class T, StorageType S, class T2, StorageType S2> void DoTestMatrixDivArithMM2b(
    const SM1& a, const SM2& b, const Matrix<T,S>& m1,
    const Matrix<T2,S2>& m2, string label)
{
  if (showstartdone) {
    cout<<"Start MM2b: "<<label<<endl;
    cout<<"b.dt = "<<b.GetDivType()<<endl;
    cout<<"m2.dt = "<<m2.GetDivType()<<endl;
  }

//#define XDEBUG

#ifdef XDEBUG
  cout<<"a = "<<Type(a)<<"  "<<a<<endl;
  cout<<"b = "<<Type(b)<<"  "<<b<<endl;
  cout<<"m1 = "<<Type(m1)<<"  "<<m1<<endl;
  cout<<"m2 = "<<Type(m2)<<"  "<<m2<<endl;
#endif

  m2.SaveDiv();
  b.SaveDiv();

  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a!=m1");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b!=m2");

  Matrix<T2,ColMajor> m2inv = Inverse(m2);
  double eps = EPS*Norm(m1)*std::max(m2.colsize(),m2.rowsize())*Norm(m2)*Norm(m2inv);

  if (CanLDivEq(m1,b)) {
    Matrix<T,S> m3 = m1;
    Matrix<T,S> m4 = m1;
    m3 /= b;
    m4 = Matrix<T,S>(m4/m2);
    Assert(Norm(m3-m4) <= eps,label+" m/=b");
  }
  if (CanRDivEq(m1,b)) {
    Matrix<T,S> m3 = m1;
    Matrix<T,S> m4 = m1;
    m3 %= b;
    m4 = Matrix<T,S>(m4%m2);
#ifdef XDEBUG
    cout<<"m3 = "<<m3<<endl;
    cout<<"m4 = "<<m4<<endl;
    cout<<"m3*b = "<<m3*b<<endl;
    cout<<"m4*m2 = "<<m4*m2<<endl;
    cout<<"Norm(m3-m4) = "<<Norm(m3-m4)<<endl;
    cout<<"EPS*Norm(m3)*Norm(m2) = "<<EPS*Norm(m3)*Norm(m2)<<endl;
#endif
    Assert(Norm(m3-m4) <= eps,label+" m%=b");
    m3 = m1;
    m4 = m1;
    m3 *= b.Inverse();
    m4 = Matrix<T,S>(m4%m2);
    Assert(Norm(m3-m4) <= eps,label+" m*=b^-1");
  }

#ifdef XDEBUG
#undef XDEBUG
#endif
  if (showstartdone) 
    cout<<"Done MM2b"<<endl;
}

template <class SM1, class SM2, class T, StorageType S, class T2, StorageType S2> void DoTestMatrixDivArithMM2a(
    DivType dt, const SM1& a, const SM2& b, const Matrix<T,S>& m1,
    const Matrix<T2,S2>& m2, string label)
{
  if (showstartdone) {
    cout<<"Start MM2a: "<<label<<endl;
    cout<<"b.dt = "<<b.GetDivType()<<endl;
    cout<<"m2.dt = "<<m2.GetDivType()<<endl;
  }
  b.DivideUsing(dt);
  DoTestMatrixDivArithMM2b(a,b,m1,m2,label);
  if (showstartdone) 
    cout<<"Done MM2a"<<endl;
}

template <class SM1, class SM2, class T, class T2> void DoTestMatrixDivArithMM2(
    DivType dt, const SM1& a, const SM2& b, const Matrix<T,RowMajor>& m1,
    const Matrix<T2,RowMajor>& m2, string label)
{
  b.DivideUsing(dt);
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a!=m1");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b!=m2");

  Matrix<T,ColMajor> m1t = m1.Transpose();
  Matrix<T,RowMajor> m1c = m1.Conjugate();
  Matrix<T,ColMajor> m1a = m1.Adjoint();
  Matrix<T2,ColMajor> m2t = m2.Transpose();
  Matrix<T2,RowMajor> m2c = m2.Conjugate();
  Matrix<T2,ColMajor> m2a = m2.Adjoint();

  DoTestMatrixDivArithMM2a(dt,a.View(),b.View(),m1,m2,label);
  DoTestMatrixDivArithMM2a(dt,Transpose(a),b.View(),m1t,m2,label+" TransA");
  DoTestMatrixDivArithMM2a(dt,Conjugate(a),b.View(),m1c,m2,label+" ConjA");
  DoTestMatrixDivArithMM2a(dt,Adjoint(a),b.View(),m1a,m2,label+" AdjA");

  DoTestMatrixDivArithMM2a(dt,a.View(),Transpose(b),m1,m2t,label+" TransB");
  DoTestMatrixDivArithMM2a(dt,Transpose(a),Transpose(b),m1t,m2t,label+" TransA TransB");
  DoTestMatrixDivArithMM2a(dt,Conjugate(a),Transpose(b),m1c,m2t,label+" ConjA TransB");
  DoTestMatrixDivArithMM2a(dt,Adjoint(a),Transpose(b),m1a,m2t,label+" AdjA TransB");

  DoTestMatrixDivArithMM2a(dt,a.View(),Conjugate(b),m1,m2c,label+" ConjB");
  DoTestMatrixDivArithMM2a(dt,Transpose(a),Conjugate(b),m1t,m2c,label+" TransA ConjB");
  DoTestMatrixDivArithMM2a(dt,Conjugate(a),Conjugate(b),m1c,m2c,label+" ConjA ConjB");
  DoTestMatrixDivArithMM2a(dt,Adjoint(a),Conjugate(b),m1a,m2c,label+" AdjA ConjB");

  DoTestMatrixDivArithMM2a(dt,a.View(),Adjoint(b),m1,m2a,label+" AdjB");
  DoTestMatrixDivArithMM2a(dt,Transpose(a),Adjoint(b),m1t,m2a,label+" TransA AdjB");
  DoTestMatrixDivArithMM2a(dt,Conjugate(a),Adjoint(b),m1c,m2a,label+" ConjA AdjB");
  DoTestMatrixDivArithMM2a(dt,Adjoint(a),Adjoint(b),m1a,m2a,label+" AdjA AdjB");
}

template <class SM, class T, class T2, StorageType S2> void DoTestMatrixDivArithMV1b(
    const Vector<T>& v, const SM& b, const Matrix<T2,S2>& m, string label)
{
  if (showstartdone) {
    cout<<"Start MV1b: "<<label<<endl;
    cout<<"b.dt = "<<b.GetDivType()<<endl;
    cout<<"m.dt = "<<m.GetDivType()<<endl;
  }

  m.SaveDiv();
  b.SaveDiv();

  Assert(Norm(b-m) <= EPS*Norm(m),label+" b!=m");

  Matrix<T2,ColMajor> minv = Inverse(m);
  double eps = EPS*Norm(v)*std::max(m.colsize(),m.rowsize())*Norm(m)*Norm(minv);

  if (CanLDiv(v,b)) {
//#define XDEBUG
#ifdef XDEBUG
    cerr<<"v = "<<Type(v)<<"  step "<<v.step()<<"  "<<v<<endl;
    cerr<<"b = "<<Type(b)<<"  "<<b<<endl;
    cerr<<"m = "<<Type(m)<<"  "<<m<<endl;
    cerr<<"v/b = "<<v/b<<endl;
    cerr<<"v/m = "<<v/m<<endl;
    cerr<<"Norm(diff) = "<<Norm((v/b)-(v/m))<<endl;
    cerr<<"eps = "<<EPS<<"*"<<Norm(v)<<"*"<<std::max(m.colsize(),m.rowsize())<<"*"<<Norm(m)<<"*"<<Norm(minv)<<"="<<eps<<endl;
#undef XDEBUG
#endif
    Assert(Norm((v/b)-(v/m)) <= eps,label+" v/m");
    Assert(Norm((b.Inverse()*v)-(v/m)) <= eps,label+" m^-1*v");
  }
  if (CanRDiv(v,b)) {
//#define XDEBUG
#ifdef XDEBUG
    cerr<<"v = "<<Type(v)<<"  step "<<v.step()<<"  "<<v<<endl;
    cerr<<"b = "<<Type(b)<<"  "<<b<<endl;
    cerr<<"m = "<<Type(m)<<"  "<<m<<endl;
    cerr<<"v%b = "<<v%b<<endl;
    cerr<<"v%m = "<<v%m<<endl;
    cerr<<"(v%b)*m = "<<(v%b)*m<<endl;
    cerr<<"(v%b)*b = "<<(v%b)*b<<endl;
    cerr<<"(v%m)*m = "<<(v%m)*m<<endl;
    cerr<<"(v%m)*b = "<<(v%m)*b<<endl;
    cerr<<"diff = "<<(v%b)-(v%m)<<endl;
    cerr<<"Norm(diff) = "<<Norm((v%b)-(v%m))<<endl;
    cerr<<"eps = "<<EPS<<"*"<<Norm(v)<<"*"<<std::max(m.colsize(),m.rowsize())<<"*"<<Norm(m)<<"*"<<Norm(minv)<<"="<<eps<<endl;
#undef XDEBUG
#endif
    Assert(Norm((v%b)-(v%m)) <= eps,label+" v%m");
    Assert(Norm((v*b.Inverse())-(v%m)) <= eps,label+" v*m^-1");
  }

  if (showstartdone) 
    cout<<"Done MV1b"<<endl;
}

template <class SM, class T, class T2, StorageType S2> void DoTestMatrixDivArithMV1a(
    DivType dt, const Vector<T>& v, const SM& b, 
    const Matrix<T2,S2>& m, string label)
{
  if (showstartdone) {
    cout<<"Start MV1a: "<<label<<endl;
    cout<<"b.dt = "<<b.GetDivType()<<endl;
    cout<<"m.dt = "<<m.GetDivType()<<endl;
  }
  b.DivideUsing(dt);
  DoTestMatrixDivArithMV1b(v,b,m,label);
  if (showstartdone) 
    cout<<"Done MV1a"<<endl;
}

template <class SM, class T, class T2> void DoTestMatrixDivArithMV1(
    DivType dt, const Vector<T>& v, const SM& b,
    const Matrix<T2,RowMajor>& m, string label)
{
  b.DivideUsing(dt);
  Assert(Norm(b-m) <= EPS*Norm(m),label+" b!=m");

  Matrix<T2,ColMajor> mt = m.Transpose();
  Matrix<T2,RowMajor> mc = m.Conjugate();
  Matrix<T2,ColMajor> ma = m.Adjoint();

  DoTestMatrixDivArithMV1a(dt,v,b.View(),m,label);
  DoTestMatrixDivArithMV1a(dt,v,Transpose(b),mt,label+" Trans");
  DoTestMatrixDivArithMV1a(dt,v,Conjugate(b),mc,label+" Conj");
  DoTestMatrixDivArithMV1a(dt,v,Adjoint(b),ma,label+" Adj");
}

template <class SM, class T, class T2, StorageType S2> void DoTestMatrixDivArithMV2b(
    const Vector<T>& v, const SM& b, const Matrix<T2,S2>& m, string label)
{
  if (showstartdone) {
    cout<<"Start MV2b: "<<label<<endl;
    cout<<"b.dt = "<<b.GetDivType()<<endl;
    cout<<"m.dt = "<<m.GetDivType()<<endl;
  }
  Assert(Norm(b-m) <= EPS*Norm(m),label+" b!=m");

  m.SaveDiv();
  b.SaveDiv();

  Matrix<T2,ColMajor> minv = Inverse(m);
  double eps = EPS*Norm(v)*std::max(m.colsize(),m.rowsize())*Norm(m)*Norm(minv);

  if (CanLDivEq(v,m)) {
    Vector<T> v1 = v;
    Vector<T> v2 = v;
    v1 /= b;
    v2 = Vector<T>(v2/m);
    Assert(Norm(v1-v2) <= eps,label+" v/=m");
  }
  if (CanRDivEq(v,m)) {
    Vector<T> v1 = v;
    Vector<T> v2 = v;
    v1 %= b;
    v2 = Vector<T>(v2%m);
    Assert(Norm(v1-v2) <= eps,label+" v%=m");
    v1 = v;
    v2 = v;
    v1 *= b.Inverse();
    v2 = Vector<T>(v2%m);
    Assert(Norm(v1-v2) <= eps,label+" v*=m^-1");
  }
  if (showstartdone) 
    cout<<"Done MV2b"<<endl;
}

template <class SM, class T, class T2, StorageType S2> void DoTestMatrixDivArithMV2a(
    DivType dt, const Vector<T>& v, const SM& b,
    const Matrix<T2,S2>& m, string label)
{
  if (showstartdone) {
    cout<<"Start MV2a: "<<label<<endl;
    cout<<"b.dt = "<<b.GetDivType()<<endl;
    cout<<"m.dt = "<<m.GetDivType()<<endl;
  }
  b.DivideUsing(dt);
  DoTestMatrixDivArithMV2b(v,b,m,label);
  if (showstartdone) 
    cout<<"Done MV2a"<<endl;
}

template <class SM, class T, class T2> void DoTestMatrixDivArithMV2(
    DivType dt, const Vector<T>& v, const SM& b,
    const Matrix<T2,RowMajor>& m, string label)
{
  b.DivideUsing(dt);
  Assert(Norm(b-m) <= EPS*Norm(m),label+" b!=m");

  Matrix<T2,ColMajor> mt = m.Transpose();
  Matrix<T2,RowMajor> mc = m.Conjugate();
  Matrix<T2,ColMajor> ma = m.Adjoint();

  DoTestMatrixDivArithMV2a(dt,v,b.View(),m,label);
  DoTestMatrixDivArithMV2a(dt,v,Transpose(b),mt,label+" TransA");
  DoTestMatrixDivArithMV2a(dt,v,Conjugate(b),mc,label+" ConjA");
  DoTestMatrixDivArithMV2a(dt,v,Adjoint(b),ma,label+" AdjA");
}

template <class SM1, class T, StorageType S, class T2> void DoTestMatrixDivArithMXb(
    const SM1& a, const Matrix<T,S>& m, T2 x, string label)
{
  if (showstartdone) {
    cout<<"Start MXb: "<<label<<endl;
    cout<<"a.dt = "<<a.GetDivType()<<endl;
    cout<<"m.dt = "<<m.GetDivType()<<endl;
  }

  m.SaveDiv();
  a.SaveDiv();

  Matrix<T,ColMajor> minv = m.Inverse();
  double normm = Norm(m);
  double eps = EPS*std::max(m.colsize(),m.rowsize())*Norm(m)*Norm(minv);

//#define XDEBUG

#ifdef XDEBUG
  cout<<"eps = "<<eps<<endl;
  cout<<"x = "<<x<<endl;
  cout<<"a = "<<Type(a)<<"  "<<a<<endl;
  cout<<"m = "<<Type(m)<<"  "<<m<<endl;
  cout<<"x/a = "<<x/a<<endl;
  cout<<"x/m = "<<x/m<<endl;
  cout<<"a*(x/a) = "<<a*(x/a)<<endl;
  cout<<"m*(x/m) = "<<m*(x/m)<<endl;
  cout<<"Norm(diff) = "<<Norm((x/a)-(x/m))<<endl;
  cout<<"eps*x = "<<eps*abs(x)<<endl;
#endif
  Assert(Norm((x/a)-(x/m)) <= eps*abs(x),label+" x/a");
  Assert(Norm((a.Inverse()*x)-(x/m)) <= eps*abs(x),label+" a^-1*x");
#ifdef XDEBUG
  cout<<"x%a = "<<x%a<<endl;
  cout<<"x%m = "<<x%m<<endl;
  cout<<"x%a*a = "<<(x%a)*a<<endl;
  cout<<"x%m*m = "<<(x%m)*m<<endl;
  cout<<"Norm(diff) = "<<Norm((x%a)-(x%m))<<endl;
#endif
  Assert(Norm((x%a)-(x%m)) <= eps*abs(x),label+" x%a");
  Assert(Norm((x*a.Inverse())-(x%m)) <= eps*abs(x),label+" x*a^-1");
#ifdef XDEBUG
  cout<<"a = "<<Type(a)<<"  "<<a<<endl;
  cout<<"1/a = "<<T2(1)/a<<endl;
  cout<<"a*(1/a) = "<<a*(T2(1)/a)<<endl;
  cout<<"(a*(1/a))*a = "<<(a*(T2(1)/a))*a<<endl;
  cout<<"(a*(1/a))*a-a = "<<(a*(T2(1)/a))*a-a<<endl;
  cout<<Norm(a*((T2(1)/a)*a)-a)<<"  "<<eps*normm<<endl;
  cout<<Norm((a*(T2(1)/a))*a-a)<<"  "<<eps*normm<<endl;
  cout<<Norm(a*((T2(1)%a)*a)-a)<<"  "<<eps*normm<<endl;
  cout<<Norm((a*(T2(1)%a))*a-a)<<"  "<<eps*normm<<endl;
  cout<<Norm((T2(1)/a)*a-Transpose((T2(1)/a)*a))<<"  "<<eps<<endl;
  cout<<Norm((T2(1)%a)*a-Transpose((T2(1)%a)*a))<<"  "<<eps<<endl;
  cout<<Norm(a*(T2(1)/a)-Transpose(a*(T2(1)/a)))<<"  "<<eps<<endl;
  cout<<Norm(a*(T2(1)%a)-Transpose(a*(T2(1)%a)))<<"  "<<eps<<endl;
#endif
  Assert(Norm(a*((T2(1)/a)*a)-a) <= eps*normm,label+" a*(1/a)*a");
  Assert(Norm((a*(T2(1)/a))*a-a) <= eps*normm,label+" a*(1/a)*a");
  Assert(Norm(a*((T2(1)%a)*a)-a) <= eps*normm,label+" a*(1%a)*a");
  Assert(Norm((a*(T2(1)%a))*a-a) <= eps*normm,label+" a*(1%a)*a");
  Assert(Norm((T2(1)/a)*a-Transpose((T2(1)/a)*a)) <= 2*eps,label+" (1/a)*a-((1/a)*a)T");
  Assert(Norm((T2(1)%a)*a-Transpose((T2(1)%a)*a)) <= 2*eps,label+" (1/a)*a-((1/a)*a)T");
  Assert(Norm(a*(T2(1)/a)-Transpose(a*(T2(1)/a))) <= 2*eps,label+" a*(1/a)-(a*(1/a))T");
  Assert(Norm(a*(T2(1)%a)-Transpose(a*(T2(1)%a))) <= 2*eps,label+" a*(1/a)-(a*(1/a))T");

#ifdef XDEBUG
#undef XDEBUG
#endif
  if (showstartdone) 
    cout<<"Done MXb"<<endl;
}

template <class SM1, class T, StorageType S, class T2> void DoTestMatrixDivArithMXa(
    DivType dt, const SM1& a, const Matrix<T,S>& m, T2 x, string label)
{
  if (showstartdone) {
    cout<<"Start MXa: "<<label<<endl;
    cout<<"a.dt = "<<a.GetDivType()<<endl;
    cout<<"m.dt = "<<m.GetDivType()<<endl;
  }
  a.DivideUsing(dt);
  DoTestMatrixDivArithMXb(a,m,x,label);
  if (showstartdone) 
    cout<<"Done MXa"<<endl;
}

template <class SM1, class T, class T2> void DoTestMatrixDivArithMX(
    DivType dt, const SM1& a, const Matrix<T,RowMajor>& m, T2 x, string label)
{
  if (showstartdone) {
    cout<<"Start MX: "<<label<<endl;
    cout<<"a.dt = "<<a.GetDivType()<<endl;
    cout<<"m.dt = "<<m.GetDivType()<<endl;
  }
  Matrix<T,ColMajor> mt = m.Transpose();
  Matrix<T,RowMajor> mc = m.Conjugate();
  Matrix<T,ColMajor> ma = m.Adjoint();
  DoTestMatrixDivArithMXa(dt,a.View(),m,x,label);
  DoTestMatrixDivArithMXa(dt,Transpose(a),mt,x,label+" Trans");
  DoTestMatrixDivArithMXa(dt,Conjugate(a),mc,x,label+" Conj");
  DoTestMatrixDivArithMXa(dt,Adjoint(a),ma,x,label+" Adj");
  if (showstartdone) 
    cout<<"Done MX"<<endl;
}

template <class T, class SM1, class SM2, class CSM1, class CSM2> void TestMatrixDivArith(
    DivType dt, const SM1& a, const SM2& b, const CSM1& ca, const CSM2& cb, 
    string label)
{
  if (showstartdone) {
    cout<<"Start Test Div: "<<label<<endl;
    cout<<"a = "<<Type(a)<<"  "<<a<<endl;
    cout<<"b = "<<Type(b)<<"  "<<b<<endl;
    cout<<"ca = "<<Type(ca)<<"  "<<ca<<endl;
    cout<<"cb = "<<Type(cb)<<"  "<<cb<<endl;
    cout<<"a.dt = "<<a.GetDivType()<<endl;
    cout<<"b.dt = "<<b.GetDivType()<<endl;
    cout<<"ca.dt = "<<ca.GetDivType()<<endl;
    cout<<"cb.dt = "<<cb.GetDivType()<<endl;
  }
  
  T x = 12;
  complex<T> z1(9,-2);

  a.SaveDiv();
  b.SaveDiv();
  ca.SaveDiv();
  cb.SaveDiv();

  Matrix<T,RowMajor> m1(a);
  Matrix<T,RowMajor> m2(b);
  Matrix<complex<T>,RowMajor> cm1(ca);
  Matrix<complex<T>,RowMajor> cm2(cb);

  m1.SaveDiv();
  m2.SaveDiv();
  cm1.SaveDiv();
  cm2.SaveDiv();

  Vector<T> v(b.colsize());
  Vector<complex<T> > cv(b.colsize());
  if (b.rowsize() > 0) {
    v = m2.col(0);
    cv = cm2.col(0);
  }
  else {
    v.SetAllTo(1);
    cv.SetAllTo(complex<T>(1,2));
  }

  Vector<T> w(b.rowsize());
  Vector<complex<T> > cw(b.rowsize());
  if (b.colsize() > 0) {
    w = m2.row(0);
    cw = cm2.row(0);
  }
  else {
    w.SetAllTo(1);
    cw.SetAllTo(complex<T>(1,2));
  }

  DoTestMatrixDivArithMX(dt,a.View(),m1,x,label+" R,R");
  DoTestMatrixDivArithMX(dt,a.View(),m1,z1,label+" R,C");
  DoTestMatrixDivArithMX(dt,ca.View(),cm1,x,label+" C,R");
  DoTestMatrixDivArithMX(dt,ca.View(),cm1,z1,label+" C,C");

  DoTestMatrixDivArithMM1(dt,b.View(),a.View(),m2,m1,label+" R,R");
  DoTestMatrixDivArithMM2(dt,b.View(),a.View(),m2,m1,label+" R,R");
  DoTestMatrixDivArithMM1(dt,b.View(),ca.View(),m2,cm1,label+" R,C");
  DoTestMatrixDivArithMM1(dt,cb.View(),a.View(),cm2,m1,label+" C,R");
  DoTestMatrixDivArithMM2(dt,cb.View(),a.View(),cm2,m1,label+" C,R");
  DoTestMatrixDivArithMM1(dt,cb.View(),ca.View(),cm2,cm1,label+" C,C");
  DoTestMatrixDivArithMM2(dt,cb.View(),ca.View(),cm2,cm1,label+" C,C");

#ifdef XTEST
  DoTestMatrixDivArithMV1(dt,v,a.View(),m1,label+" R,R");
  DoTestMatrixDivArithMV1(dt,v,ca.View(),cm1,label+" R,C");
  DoTestMatrixDivArithMV1(dt,cv,a.View(),m1,label+" C,R");
  DoTestMatrixDivArithMV1(dt,cv,ca.View(),cm1,label+" C,C");
  DoTestMatrixDivArithMV2(dt,v,a.View(),m1,label+" R,R");
  DoTestMatrixDivArithMV2(dt,cv,a.View(),m1,label+" C,R");
  DoTestMatrixDivArithMV2(dt,cv,ca.View(),cm1,label+" C,C");

  DoTestMatrixDivArithMV1(dt,w,a.View(),m1,label+" R,R");
  DoTestMatrixDivArithMV1(dt,w,ca.View(),cm1,label+" R,C");
  DoTestMatrixDivArithMV1(dt,cw,a.View(),m1,label+" C,R");
  DoTestMatrixDivArithMV1(dt,cw,ca.View(),cm1,label+" C,C");
  DoTestMatrixDivArithMV2(dt,w,a.View(),m1,label+" R,R");
  DoTestMatrixDivArithMV2(dt,cw,a.View(),m1,label+" C,R");
  DoTestMatrixDivArithMV2(dt,cw,ca.View(),cm1,label+" C,C");
#endif
}

