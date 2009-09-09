//#define SHOWSTARTDONE
#define DONORM2

#include "TMV.h"
using tmv::Matrix;
using tmv::Vector;
using tmv::StorageType;
using tmv::RowMajor;
using tmv::ColMajor;
using tmv::Type;

template <class SM1, class SM2> bool CanAdd(const SM1& a, const SM2& b)
{ return a.colsize() == b.colsize() && a.rowsize() == b.rowsize(); }

template <class SM1, class SM2> bool CanAddEq(const SM1& a, const SM2& b)
{ return CanAdd(a,b); }

template <class SM1> bool CanAddX(const SM1& a)
{ return a.IsSquare(); }

template <class SM1> bool CanAddEqX(const SM1& a)
{ return a.IsSquare(); }

template <class SM1> bool CanMultEqX(const SM1& a)
{ return true; }

template <class SM1, class SM2> bool CanMult(const SM1& a, const SM2& b)
{ return a.rowsize() == b.colsize(); }

template <class SM, class T> bool CanMult(const SM& m, const Vector<T>& v)
{ return m.rowsize() == v.size(); }

template <class SM, class T> bool CanMult(const Vector<T>& v, const SM& m)
{ return v.size() == m.colsize(); }

template <class SM1, class SM2> bool CanMultEq(const SM1& a, const SM2& b)
{ return CanMult(a,b) && b.IsSquare(); }

template <class SM1, class SM2> bool CanMultEq2(const SM1& a, const SM2& b)
{ return CanMult(a,b) && a.IsSquare(); }

template <class SM, class T, StorageType S, class V> void DoTestMV1a(
    const SM& a, const Matrix<T,S>& m, const V& v, string label)
{
#ifdef SHOWSTARTDONE
  cout<<"Start MV1a"<<endl;
#endif
  //cerr<<"a = "<<Type(a)<<" = "<<a<<endl;
  //cerr<<"m = "<<Type(m)<<" = "<<m<<endl;
  //cerr<<"v = "<<Type(v)<<" = "<<v<<endl;
  //cerr<<"Norm(a-m) = "<<Norm(a-m)<<endl;
  Assert(Norm(a-m) <= EPS*Norm(m),label+" a != m");
  if (CanMult(m,ColVector(v))) {
    Assert(Norm((a*v)-(m*v)) <= EPS*Norm(m)*Norm(v),label+" a*v");
  }
  if (CanMult(RowVector(v),m)) {
    Assert(Norm((v*a)-(v*m)) <= EPS*Norm(v)*Norm(m),label+" v*a");
  }
#ifdef SHOWSTARTDONE
  cout<<"Done MV1a"<<endl;
#endif
}

template <class SM, class T, class V> void DoTestMV1(
    const SM& a, const Matrix<T,RowMajor>& m, const V& v,
    string label)
{
  Matrix<T,ColMajor> mt = m.Transpose();
  Matrix<T,RowMajor> mc = m.Conjugate();
  Matrix<T,ColMajor> ma = m.Adjoint();

  DoTestMV1a(a.View(),m,v,label);
  DoTestMV1a(Transpose(a),mt,v,label+" Trans");
  DoTestMV1a(Conjugate(a),mc,v,label+" Conj");
  DoTestMV1a(Adjoint(a),ma,v,label+" Adj");
}

template <class SM, class T, StorageType S, class V, class T2> void DoTestMV2a(
    const SM& a, const Matrix<T,S>& m, const V& v1, Vector<T2>& v2,
    string label)
{
#ifdef SHOWSTARTDONE
  cout<<"Start MV2a"<<endl;
#endif
  Vector<T2> v0 = v1;
  if (CanMultEq(v2,m)) {
    v2 = v1;
    v1 *= a;
    v2 = Vector<T2>(v2*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v*=a");
    v1 = v0;
    v2 = v0;
    v1 = v1*a;
    v2 = Vector<T2>(v2*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=v*a");
    v1 = v0;
    v2 = v0;
    v2 = v1;
    v1 = v2*a;
    v2 = Vector<T2>(v2*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=v*a");
    v1 = v0;
    v2 = v0;
  }
  if (CanMultEq2(m,v2)) {
    v2 = v1;
    v1 = a * v1;
    v2 = Vector<T2>(m*v2);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=a*v");
    v1 = v0;
    v2 = v0;
    v2 = v1;
    v1 = a * v2;
    v2 = Vector<T2>(m*v2);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=a*v");
    v1 = v0;
    v2 = v0;
  }
#ifdef SHOWSTARTDONE
  cout<<"Done MV2a"<<endl;
#endif
}

template <class SM, class T, class V, class T2> void DoTestMV2(
    const SM& a, const Matrix<T,RowMajor>& m, const V& v1, const Vector<T2>& v0,
    string label)
{
  Matrix<T,ColMajor> mt = m.Transpose();
  Matrix<T,RowMajor> mc = m.Conjugate();
  Matrix<T,ColMajor> ma = m.Adjoint();

  v1 = v0;
  Vector<T2> v2 = v0;
  DoTestMV2a(a.View(),m,v1,v2,label);
  v1 = v0;
  v2 = v0;
  DoTestMV2a(Transpose(a),mt,v1,v2,label+" Trans");
  v1 = v0;
  v2 = v0;
  DoTestMV2a(Conjugate(a),mc,v1,v2,label+" Conj");
  v1 = v0;
  v2 = v0;
  DoTestMV2a(Adjoint(a),ma,v1,v2,label+" Adj");
}

template <class SM, class V1, class V2, class T> void DoTestMV3a(
    const SM& a, const V1& v1, const V2& v2, Vector<T>& v3,
    string label)
{
#ifdef SHOWSTARTDONE
  cout<<"Start MV3a"<<endl;
#endif
  if (CanMultEq(v3,a)) {
    //cerr<<"a = "<<Type(a)<<"  "<<a<<endl;
    //cerr<<"v1 = "<<v1<<endl;
    v2 = v1*a;
    v3 = v1*a;
    //cerr<<"v2 = "<<v2<<endl;
    //cerr<<"v3 = "<<v3<<endl;
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2=v1*a");
    v2 = v1;
    v3 = v2*a;
    v2 = v2*a;
    //cerr<<"v3=v2*a = "<<v3<<endl;
    //cerr<<"v2=v2*a = "<<v2<<endl;
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2=v2*a");
    v3 = v1;
    v2 = v3*a;
    v3 = Vector<T>(v3*a);
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2=v3*a");
  }
  if (CanMultEq2(a,v3)) {
    v2 = a*v1;
    v3 = a*v1;
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2=a*v1");
    v2 = v1;
    v3 = a*v2;
    v2 = a*v2;
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2=a*v2");
    v3 = v1;
    v2 = a*v3;
    v3 = Vector<T>(a*v3);
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2=a*v3");
  }
#ifdef SHOWSTARTDONE
  cout<<"Done MV3a"<<endl;
#endif
}

template <class SM, class V1, class V2, class T2> void DoTestMV3(
    const SM& a, const V1& v1, const V2& v2, const Vector<T2>& v0,
    string label)
{
#ifdef SHOWSTARTDONE
  cout<<"Start MV3"<<endl;
#endif
  v1 = v0;
  v2 = v0;
  Vector<T2> v3 = v0;

  //cerr<<"v0 = "<<Type(v0)<<v0<<endl;
  //cerr<<"v1 = "<<Type(v1)<<v1<<endl;
  //cerr<<"v2 = "<<Type(v2)<<v2<<endl;
  //cerr<<"v3 = "<<Type(v3)<<v3<<endl;

  DoTestMV3a(a.View(),v1,v2,v3,label);
  v1 = v0;
  v2 = v0;
  v3 = v0;
  DoTestMV3a(Transpose(a),v1,v2,v3,label+" Trans");
  v1 = v0;
  v2 = v0;
  v3 = v0;
  DoTestMV3a(Conjugate(a),v1,v2,v3,label+" Conj");
  v1 = v0;
  v2 = v0;
  v3 = v0;
  DoTestMV3a(Adjoint(a),v1,v2,v3,label+" Adj");
#ifdef SHOWSTARTDONE
  cout<<"Done MV3"<<endl;
#endif
}

template <class SM1, class SM2, class T, StorageType S, class T2, StorageType S2> void DoTestMM1a(
    const SM1& a, const SM2& b, const Matrix<T,S>& m1, const Matrix<T2,S2>& m2,
    string label)
{
#ifdef SHOWSTARTDONE
  cout<<"Start MM1a"<<endl;
#endif
  //cerr<<"a = "<<Type(a)<<" = "<<a<<endl;
  //cerr<<"b = "<<Type(b)<<" = "<<b<<endl;
  //cerr<<"m1 = "<<Type(m1)<<" = "<<m1<<endl;
  //cerr<<"m2 = "<<Type(m2)<<" = "<<m2<<endl;
  //cerr<<"Norm(a-m1) = "<<Norm(a-m1)<<endl;
  //cerr<<"Norm(b-m2) = "<<Norm(b-m2)<<endl;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a != m1");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b != m2");
  if (CanAdd(m1,m2)) {
    //cerr<<"a-m2 = "<<a-m2<<endl;
    //cerr<<"m1-m2 = "<<m1-m2<<endl;
    Assert(Norm((a-m2)-(m1-m2)) <= EPS*(Norm(m1)+Norm(m2)),label+" a-m");
    //cerr<<"m1-b = "<<m1-b<<endl;
    //cerr<<"m1-m2 = "<<m1-m2<<endl;
    Assert(Norm((m1-b)-(m1-m2)) <= EPS*(Norm(m1)+Norm(m2)),label+" m-b");
    //cerr<<"a-b = "<<a-b<<endl;
    //cerr<<"m1-m2 = "<<m1-m2<<endl;
    Assert(Norm((a-b)-(m1-m2)) <= EPS*(Norm(m1)+Norm(m2)),label+" a-b");
    Assert(Norm((a+m2)-(m1+m2)) <= EPS*(Norm(m1)+Norm(m2)),label+" a+m");
    Assert(Norm((m1+b)-(m1+m2)) <= EPS*(Norm(m1)+Norm(m2)),label+" m+b");
    Assert(Norm((a+b)-(m1+m2)) <= EPS*(Norm(m1)+Norm(m2)),label+" a+b");
  }
  if (CanMult(m1,m2)) {
    //cerr<<"m1*b = "<<m1*b<<endl;
    //cerr<<"m1*m2 = "<<m1*m2<<endl;
    Assert(Norm((m1*b)-(m1*m2)) <= EPS*Norm(m1)*Norm(m2),label+" m*b");
    //cerr<<"a*m2 = "<<a*m2<<endl;
    //cerr<<"m1*m2 = "<<m1*m2<<endl;
    Assert(Norm((a*m2)-(m1*m2)) <= EPS*Norm(m1)*Norm(m2),label+" a*m");
    //cerr<<"a*b = "<<a*b<<endl;
    //cerr<<"m1*m2 = "<<m1*m2<<endl;
    //cerr<<"a*b-m1*m2 = "<<a*b-m1*m2<<endl;
    //cerr<<"Norm(a*b-m1*m2) = "<<Norm(a*b-m1*m2)<<endl;
    //cerr<<"EPS*Norm(m1*m2) = "<<EPS*Norm(m1*m2)<<endl;
    Assert(Norm((a*b)-(m1*m2)) <= EPS*Norm(m1)*Norm(m2),label+" a*b");
  }
#ifdef SHOWSTARTDONE
  cout<<"Done MM1a"<<endl;
#endif
}

template <class SM1, class SM2, class T, class T2> void DoTestMM1(
    const SM1& a, const SM2& b, const Matrix<T,RowMajor>& m1,
    const Matrix<T2,RowMajor>& m2, string label)
{
  Matrix<T,ColMajor> m1t = m1.Transpose();
  Matrix<T,RowMajor> m1c = m1.Conjugate();
  Matrix<T,ColMajor> m1a = m1.Adjoint();
  Matrix<T2,ColMajor> m2t = m2.Transpose();
  Matrix<T2,RowMajor> m2c = m2.Conjugate();
  Matrix<T2,ColMajor> m2a = m2.Adjoint();

  DoTestMM1a(a.View(),b.View(),m1,m2,label);
  DoTestMM1a(Transpose(a),b.View(),m1t,m2,label+" TransA");
  DoTestMM1a(Conjugate(a),b.View(),m1c,m2,label+" ConjA");
  DoTestMM1a(Adjoint(a),b.View(),m1a,m2,label+" AdjA");

  DoTestMM1a(a.View(),Transpose(b),m1,m2t,label+" TransB");
  DoTestMM1a(Transpose(a),Transpose(b),m1t,m2t,label+" TransA TransB");
  DoTestMM1a(Conjugate(a),Transpose(b),m1c,m2t,label+" ConjA TransB");
  DoTestMM1a(Adjoint(a),Transpose(b),m1a,m2t,label+" AdjA TransB");

  DoTestMM1a(a.View(),Conjugate(b),m1,m2c,label+" ConjB");
  DoTestMM1a(Transpose(a),Conjugate(b),m1t,m2c,label+" TransA ConjB");
  DoTestMM1a(Conjugate(a),Conjugate(b),m1c,m2c,label+" ConjA ConjB");
  DoTestMM1a(Adjoint(a),Conjugate(b),m1a,m2c,label+" AdjA ConjB");

  DoTestMM1a(a.View(),Adjoint(b),m1,m2a,label+" AdjB");
  DoTestMM1a(Transpose(a),Adjoint(b),m1t,m2a,label+" TransA AdjB");
  DoTestMM1a(Conjugate(a),Adjoint(b),m1c,m2a,label+" ConjA AdjB");
  DoTestMM1a(Adjoint(a),Adjoint(b),m1a,m2a,label+" AdjA AdjB");
}

template <class SM1, class SM2, class T, StorageType S, class T2, StorageType S2> void DoTestMM2a(
    const SM1& a, const SM2& b, Matrix<T,S>& m1, const Matrix<T2,S2>& m2,
    string label)
{
#ifdef SHOWSTARTDONE
  cout<<"Start MM2a"<<endl;
#endif
  //cerr<<"a = "<<Type(a)<<"  "<<a<<endl;
  //cerr<<"b = "<<Type(b)<<"  "<<b<<endl;
  //cerr<<"m1 = "<<Type(m1)<<"  "<<m1<<endl;
  //cerr<<"m2 = "<<Type(m2)<<"  "<<m2<<endl;
  Matrix<T,S> m3 = m1;
  Matrix<T,S> m4 = m1;

  if (CanAddEq(m3,b)) {
    double normm = Norm(m1)+Norm(m2);
    m3 += b;
    m4 = Matrix<T,S>(m4+m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m += b");
    m3 = m4 = m1;
    m3 -= b;
    m4 = Matrix<T,S>(m4-m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m -= b");
    m3 = m4 = m1;
    m3 = m3 + b;
    m4 = Matrix<T,S>(m4+m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m = m+b");
    m3 = m4 = m1;
    m3 = m3 - b;
    m4 = Matrix<T,S>(m4-m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m = m-b");
    m4 = m3;
    m3 = b + m3;
    m4 = Matrix<T,S>(m2+m4);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m = b+m");
    m3 = m4 = m1;
    //cerr<<"m3 = "<<m3<<endl;
    //cerr<<"m4 = "<<m4<<endl;
    m3 = b - m3;
    m4 = Matrix<T,S>(m2-m4);
    //cerr<<"m3 = b-m3 "<<m3<<endl;
    //cerr<<"m4 = m2-m4 "<<m4<<endl;
    Assert(Norm(m3-m4) <= EPS*Norm(m3),label+" m = b-m");
  }
  if (CanAddEq(a,b)) {
    m1 = a;
    double normm = Norm(m1)+Norm(m2);
    //cerr<<"a = "<<Type(a)<<a<<endl;
    //cerr<<"m1 = "<<Type(m1)<<m1<<endl;
    //cerr<<"b = "<<Type(b)<<b<<endl;
    a += b;
    //cerr<<"a+=b = "<<a<<endl;
    //cerr<<"b = "<<b<<endl;
    m1 = Matrix<T,S>(m1+b);
    //cerr<<"m1+=b = "<<m1<<endl;
    Assert(Norm(a-m1) <= EPS*normm,label+" a += b");
    m1 = a;
    normm = Norm(m1)+Norm(m2);
    a -= b;
    m1 = Matrix<T,S>(m1-b);
    Assert(Norm(a-m1) <= EPS*normm,label+" a -= b");
    m1 = a;
    normm = Norm(m1)+Norm(m2);
    a = a+b;
    m1 = Matrix<T,S>(m1+b);
    Assert(Norm(a-m1) <= EPS*normm,label+" a = a+b");
    m1 = a;
    normm = Norm(m1)+Norm(m2);
    a = a-b;
    m1 = Matrix<T,S>(m1-b);
    Assert(Norm(a-m1) <= EPS*normm,label+" a = a-b");
    m1 = a;
    normm = Norm(m1)+Norm(m2);
    a = b+a; 
    m1 = Matrix<T,S>(m2+m1);
    Assert(Norm(a-m1) <= EPS*normm,label+" a = b+a");
    m1 = a;
    normm = Norm(m1)+Norm(m2);
    a = b-a;
    m1 = Matrix<T,S>(m2-m1);
    Assert(Norm(a-m1) <= EPS*normm,label+" a = b-a");
  }
  if (CanMultEq(m4,m2)) {
    double normm = Norm(m1)*Norm(m2);
    m3 = m4 = m1;
    m3 *= b;
    m4 = Matrix<T,S>(m4*m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m *= b");
    m3 = m4 = m1;
    m3 = m3*b;
    m4 = Matrix<T,S>(m4*m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m = m*b");
  }
  if (CanMultEq(a,b)) {
    m1 = a;
    //cerr<<"a = "<<a<<endl;
    //cerr<<"m1 = "<<m1<<endl;
    double normm = Norm(m1)*Norm(m2);
    a *= b;
    m1 = Matrix<T,S>(m1*m2);
    //cerr<<"a*=b = "<<a<<endl;
    //cerr<<"m1*=b = "<<m1<<endl;
    Assert(Norm(a-m1) <= EPS*normm,label+" a *= b");
    a /= Norm1(b);
    m1 = a;
    normm = Norm(m1)*Norm(m2);
    //cerr<<"a = "<<a<<endl;
    //cerr<<"m1 = "<<m1<<endl;
    a = a*b;
    m1 = Matrix<T,S>(m1*m2);
    //cerr<<"a = a*b "<<a<<endl;
    //cerr<<"m1 = m1*m2 "<<m1<<endl;
    //cerr<<"Norm(diff) = "<<Norm(a-m1)<<endl;
    //cerr<<"EPS*Norm(m1) = "<<EPS*Norm(m1)<<endl;
    Assert(Norm(a-m1) <= EPS*normm,label+" a = a*b");
    a /= Norm1(b);
    m1 = a;
  }
  if (CanMultEq2(m2,m4)) {
    m3 = m4 = m1;
    double normm = Norm(m1)*Norm(m2);
    //cerr<<"m3 = "<<m3<<endl;
    //cerr<<"m4 = "<<m4<<endl;
    //cerr<<"b = "<<b<<endl;
    //cerr<<"m2 = "<<m2<<endl;
    m3 = b * m3;
    m4 = Matrix<T,S>(m2*m4);
    //cerr<<"m3 = b*m3 "<<m3<<endl;
    //cerr<<"m4 = m2*m4 "<<m4<<endl;
    //cerr<<"Norm(diff) = "<<Norm(m3-m4)<<"  "<<EPS*Norm(m3)<<endl;
    Assert(Norm(m3-m4) <= EPS*normm,label+" m = b*m");
  }
  if (CanMultEq2(b,a)) {
    m1 = a;
    double normm = Norm(m1)*Norm(m2);
    //cerr<<"a = "<<a<<endl;
    //cerr<<"m1 = "<<m1<<endl;
    //cerr<<"b = "<<b<<endl;
    //cerr<<"m2 = "<<m2<<endl;
    a = b * a;
    m1 = Matrix<T,S>(m2*m1);
    //cerr<<"a = b*a = "<<a<<endl;
    //cerr<<"m1 = m2*m1 = "<<m1<<endl;
    //cerr<<"Norm(diff) = "<<Norm(a-m1)<<"  "<<EPS*Norm(m1)<<endl;
    Assert(Norm(a-m1) <= EPS*normm,label+" a = b*a");
    a /= Norm1(b);
    m1 = a;
  }
#ifdef SHOWSTARTDONE
  cout<<"Done MM2a"<<endl;
#endif
}

template <class SM1, class SM2, class T, class T2, class BaseSM1> void DoTestMM2(
    const SM1& a, const SM2& b, Matrix<T,RowMajor>& m1, 
    Matrix<T2,RowMajor>& m2, const BaseSM1& basea, string label)
{
  m1 = a = basea;
  DoTestMM2a(a.View(),b.View(),m1,m2,label);
  m1 = a = basea;
  Matrix<T,ColMajor> m1t = m1.Transpose();
  DoTestMM2a(Transpose(a),b.View(),m1t,m2,label+" TransA");
  m1 = a = basea;
  Matrix<T,RowMajor> m1c = m1.Conjugate();
  DoTestMM2a(Conjugate(a),b.View(),m1c,m2,label+" ConjA");
  m1 = a = basea;
  Matrix<T,ColMajor> m1a = m1.Adjoint();
  DoTestMM2a(Adjoint(a),b.View(),m1a,m2,label+" AdjA");
  m1 = a = basea;

  Matrix<T2,ColMajor> m2t = m2.Transpose();
  DoTestMM2a(a.View(),Transpose(b),m1,m2t,label+" TransB");
  m1 = a = basea;
  m1t = m1.Transpose();
  DoTestMM2a(Transpose(a),Transpose(b),m1t,m2t,label+" TransA TransB");
  m1 = a = basea;
  m1c = m1.Conjugate();
  DoTestMM2a(Conjugate(a),Transpose(b),m1c,m2t,label+" ConjA TransB");
  m1 = a = basea;
  m1a = m1.Adjoint();
  DoTestMM2a(Adjoint(a),Transpose(b),m1a,m2t,label+" AdjA TransB");
  m1 = a = basea;

  Matrix<T2,RowMajor> m2c = m2.Conjugate();
  DoTestMM2a(a.View(),Conjugate(b),m1,m2c,label+" ConjB");
  m1 = a = basea;
  m1t = m1.Transpose();
  DoTestMM2a(Transpose(a),Conjugate(b),m1t,m2c,label+" TransA ConjB");
  m1 = a = basea;
  m1c = m1.Conjugate();
  DoTestMM2a(Conjugate(a),Conjugate(b),m1c,m2c,label+" ConjA ConjB");
  m1 = a = basea;
  m1a = m1.Adjoint();
  DoTestMM2a(Adjoint(a),Conjugate(b),m1a,m2c,label+" AdjA ConjB");
  m1 = a = basea;

  Matrix<T2,ColMajor> m2a = m2.Adjoint();
  DoTestMM2a(a.View(),Adjoint(b),m1,m2a,label+" AdjB");
  m1 = a = basea;
  m1t = m1.Transpose();
  DoTestMM2a(Transpose(a),Adjoint(b),m1t,m2a,label+" TransA AdjB");
  m1 = a = basea;
  m1c = m1.Conjugate();
  DoTestMM2a(Conjugate(a),Adjoint(b),m1c,m2a,label+" ConjA AdjB");
  m1 = a = basea;
  m1a = m1.Adjoint();
  DoTestMM2a(Adjoint(a),Adjoint(b),m1a,m2a,label+" AdjA AdjB");
}

template <class SM1, class T, StorageType S, class T2> void DoTestMX1a(
    const SM1& a, const Matrix<T,S>& m1, T2 x, string label)
{
#ifdef SHOWSTARTDONE
  cout<<"Start MX1a"<<endl;
#endif
  //cerr<<"a = "<<Type(a)<<" = "<<a<<endl;
  //cerr<<"m1 = "<<Type(m1)<<" = "<<m1<<endl;
  //cerr<<"x = "<<Type(x)<<" = "<<x<<endl;
  //cerr<<"Norm(a-m1) = "<<Norm(a-m1)<<endl;
  double normm = Norm(m1);
  if (CanAddX(m1)) {
    Assert(Norm(a-m1) <= EPS*normm,label+" a != m1");
    Assert(Norm((x-a)-(x-m1)) <= EPS*normm,label+" x-a");
    Assert(Norm((a-x)-(m1-x)) <= EPS*normm,label+" a-x");
    Assert(Norm((x+a)-(x+m1)) <= EPS*normm,label+" x+a");
    Assert(Norm((a+x)-(m1+x)) <= EPS*normm,label+" a+x");
  }
  Assert(Norm((x*a)-(x*m1)) <= EPS*normm*abs(x),label+" x*a");
  Assert(Norm((a*x)-(m1*x)) <= EPS*normm*abs(x),label+" a*x");
  Assert(Norm((a/x)-(m1/x)) <= EPS*normm*abs(x),label+" a/x");
#ifdef SHOWSTARTDONE
  cout<<"Done MX1a"<<endl;
#endif
}

template <class SM1, class T, class T2> void DoTestMX1(
    const SM1& a, Matrix<T,RowMajor>& m1, T2 x, string label)
{
  DoTestMX1a(a.View(),m1,x,label);
  Matrix<T,ColMajor> m1t = m1.Transpose();
  DoTestMX1a(Transpose(a),m1t,x,label+" Trans");
  Matrix<T,RowMajor> m1c = m1.Conjugate();
  DoTestMX1a(Conjugate(a),m1c,x,label+" Conj");
  Matrix<T,ColMajor> m1a = m1.Adjoint();
  DoTestMX1a(Adjoint(a),m1a,x,label+" Adj");
}

template <class SM1, class T, StorageType S, class T2> void DoTestMX2a(
    const SM1& a, Matrix<T,S>& m1, T2 x, string label)
{
#ifdef SHOWSTARTDONE
  cout<<"Start MX2a"<<endl;
#endif
  //cerr<<"a = "<<Type(a)<<"  "<<a<<endl;
  //cerr<<"m1 = "<<Type(m1)<<"  "<<m1<<endl;
  //cerr<<"x = "<<Type(x)<<"  "<<x<<endl;
  if (CanAddEqX(a)) {
    Assert(Norm((a+=x)-(m1+=x)) <= EPS*Norm(m1),label+" a += x");
    Assert(Norm((a-=x)-(m1-=x)) <= EPS*Norm(m1),label+" a -= x");
    a = a+x; 
    m1 = Matrix<T,S>(m1+x);
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = a+x");
    m1 = a;
    a = a-x;
    m1 = Matrix<T,S>(m1-x);
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = a-x");
    m1 = a;
    a = x+a;
    m1 = Matrix<T,S>(x+m1);
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = x+a");
    m1 = a;
    a = x-a; 
    m1 = Matrix<T,S>(x-m1);
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = x-a");
  }
  if (CanMultEqX(a)) {
    m1 = a;
    //cerr<<"a = "<<Type(a)<<"  "<<a<<endl;
    //cerr<<"m1 = "<<Type(m1)<<"  "<<m1<<endl;
    a*=x;
    m1*=x;
    //cerr<<"a*=x = "<<a<<endl;
    //cerr<<"m1*=x = "<<m1<<endl;
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a *= x");
    //Assert(Norm((a*=x)-(m1*=x)) <= EPS*Norm(m1),label+" a *= x");
    m1 = a;
    Assert(Norm((a/=x)-(m1/=x)) <= EPS*Norm(m1),label+" a /= x");
    m1 = a;
    a = a*x;
    m1 = Matrix<T,S>(m1*x);
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = a*x");
    m1 = a;
    a = a/x;
    m1 = Matrix<T,S>(m1/x);
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = a/x");
    m1 = a;
    a = x*a; 
    m1 = Matrix<T,S>(x*m1);
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = x*a");
    m1 = a;
    a /= x;
    m1 = Matrix<T,S>(m1/x);
  }
#ifdef SHOWSTARTDONE
  cout<<"Done MX2a"<<endl;
#endif
}

template <class SM1, class T, class T2, class BaseSM1> void DoTestMX2(
    const SM1& a, Matrix<T,RowMajor>& m1, T2 x,
    const BaseSM1& basea, string label)
{
  m1 = a = basea;
  DoTestMX2a(a.View(),m1,x,label);
  m1 = a = basea;
  Matrix<T,ColMajor> m1t = m1.Transpose();
  DoTestMX2a(Transpose(a),m1t,x,label+" Trans");
  m1 = a = basea;
  Matrix<T,RowMajor> m1c = m1.Conjugate();
  DoTestMX2a(Conjugate(a),m1c,x,label+" Conj");
  m1 = a = basea;
  Matrix<T,ColMajor> m1a = m1.Adjoint();
  DoTestMX2a(Adjoint(a),m1a,x,label+" Adj");
}

template <class SM1, class T, StorageType S> void DoTestMa(
    const SM1& a, const Matrix<T,S>& m, string label)
{
#ifdef SHOWSTARTDONE
  cout<<"Start Ma"<<endl;
#endif
  //cerr<<"a = "<<Type(a)<<" = "<<a<<endl;
  //cerr<<"m = "<<Type(m)<<" = "<<m<<endl;
  //cerr<<"a-m = "<<a-m<<endl;
  //cerr<<"Norm(a-m) = "<<Norm(a-m)<<endl;
  Assert(Norm(a-m) <= EPS*Norm(m),label+" a != m");
  //cerr<<"Trace(a) = "<<Trace(a)<<endl;
  Assert(abs(Trace(a)-Trace(m)) <= EPS*abs(Trace(m)),label+" Trace");
  //cerr<<"NormF(a) = "<<NormF(a)<<endl;
  Assert(abs(NormF(a)-NormF(m)) <= EPS*NormF(m),label+" NormF");
  //cerr<<"Norm(a) = "<<Norm(a)<<endl;
  Assert(abs(Norm(a)-Norm(m)) <= EPS*Norm(m),label+" Norm");
  //cerr<<"Norm1(a) = "<<Norm1(a)<<endl;
  //cerr<<"Norm1(m) = "<<Norm1(m)<<endl;
  Assert(abs(Norm1(a)-Norm1(m)) <= EPS*Norm1(m),label+" Norm1");
#ifdef DONORM2
  a.DivideUsing(tmv::SV);
  m.DivideUsing(tmv::SV);
  //cerr<<"Norm2(a) = "<<Norm2(a)<<endl;
  Assert(abs(Norm2(a)-Norm2(m)) <= EPS*Norm2(m),label+" Norm2");
#endif
  //cerr<<"NormInf(a) = "<<NormInf(a)<<endl;
  //cerr<<"NormInf(m) = "<<NormInf(m)<<endl;
  //cerr<<"abs(diff) = "<<abs(NormInf(a)-NormInf(m))<<endl;
  //cerr<<"eps*norminf = "<<EPS*NormInf(m)<<endl;
  Assert(abs(NormInf(a)-NormInf(m)) <= EPS*NormInf(m),label+" NormInf");
  //cerr<<"Norm(aT-mT) = "<<Norm(Transpose(a)-Transpose(m))<<endl;
  Assert(Norm(Transpose(a)-Transpose(m)) <= EPS*Norm(m),label+" Transpose");
  //cerr<<"Conjugate(a).diag = "<<Conjugate(a).diag()<<endl;
  //cerr<<"Conjugate(m).diag = "<<Conjugate(m).diag()<<endl;
  //cerr<<"Norm(a*-m*) = "<<Norm(Conjugate(a)-Conjugate(m))<<endl;
  Assert(Norm(Conjugate(a)-Conjugate(m)) <= EPS*Norm(m),label+" Conjugate");
  //cerr<<"Norm(at-mt) = "<<Norm(Adjoint(a)-Adjoint(m))<<endl;
  Assert(Norm(Adjoint(a)-Adjoint(m)) <= EPS*Norm(m),label+" Adjoint");
#ifdef SHOWSTARTDONE
  cout<<"Done Ma"<<endl;
#endif
}

template <class SM1, class T> void DoTestM(
    const SM1& a, const Matrix<T,RowMajor>& m, string label)
{
  DoTestMa(a.View(),m,label);
  Matrix<T,ColMajor> mt = m.Transpose();
  DoTestMa(Transpose(a),mt,label+" Trans");
  Matrix<T,RowMajor> mc = m.Conjugate();
  DoTestMa(Conjugate(a),mc,label+" Conj");
  Matrix<T,ColMajor> ma = m.Adjoint();
  DoTestMa(Adjoint(a),ma,label+" Adj");
}

template <class T, class BaseM, class BaseCM, class SM1, class SM2, class CSM1, class CSM2> void TestMatrixArith(
    const SM1& a, const SM2& b, const CSM1& ca, const CSM2& cb, string label)
{
#ifdef SHOWSTARTDONE
  cout<<"Start TestMatrixArith"<<endl;
  cout<<"a = "<<Type(a)<<" "<<a<<endl;
  cout<<"b = "<<Type(b)<<" "<<b<<endl;
  cout<<"ca = "<<Type(ca)<<" "<<ca<<endl;
  cout<<"cb = "<<Type(cb)<<" "<<cb<<endl;
#endif
  T x = 12;
  complex<T> z(9,-2);

  const BaseM a0 = a;
  const BaseCM ca0 = ca;

  Matrix<T,RowMajor> m1 = a;
  Matrix<T,RowMajor> m2 = b;
  Matrix<complex<T>,RowMajor> cm1 = ca;
  Matrix<complex<T>,RowMajor> cm2 = cb;

  Vector<T> v0(a.colsize(),T(4));
  Vector<complex<T> > cv0(a.colsize(),complex<T>(4,-4));

  Vector<T> v = v0;
  Vector<complex<T> > cv = cv0;

  Matrix<T,RowMajor> vx(a.colsize(),5);
  Matrix<complex<T>,RowMajor> cvx(a.colsize(),5);
  vx.col(0) = v;
  cvx.col(0) = cv;
  cvx.col(1) = cv;
  DoTestM(a.View(),m1,label+" R");
  DoTestM(ca.View(),cm1,label+" C");

  DoTestMX1(a.View(),m1,x,label+" R,R");
  DoTestMX2(a.View(),m1,x,a0,label+" R,R");
  m1 = a = a0;
  DoTestMX1(a.View(),m1,z,label+" R,C");
  cm1 = ca = ca0;
  DoTestMX1(ca.View(),cm1,x,label+" C,R");
  DoTestMX2(ca.View(),cm1,x,ca0,label+" C,R");
  cm1 = ca = ca0;
  DoTestMX1(ca.View(),cm1,z,label+" C,C");
  DoTestMX2(ca.View(),cm1,z,ca0,label+" C,C");
  cm1 = ca = ca0;

  DoTestMV1(a.View(),m1,v,label+" R,R");
  DoTestMV2(a.View(),m1,v.View(),v0,label+" R,R");
  DoTestMV1(a.View(),m1,cv,label+" R,C");
  DoTestMV2(a.View(),m1,cv.View(),cv0,label+" R,C");
  DoTestMV1(ca.View(),cm1,v,label+" C,R");
  DoTestMV1(ca.View(),cm1,cv,label+" C,C");
  DoTestMV2(ca.View(),cm1,cv.View(),cv0,label+" C,C");

  DoTestMV1(a.View(),m1,vx.col(0),label+" Step R,R");
  DoTestMV2(a.View(),m1,vx.col(0),v0,label+" Step R,R");
  DoTestMV3(a.View(),vx.col(0),v.View(),v0,label+" Step R,R");
  DoTestMV3(a.View(),vx.col(0),vx.col(1),v0,label+" Step R,R");
  DoTestMV1(a.View(),m1,cvx.col(0),label+" Step R,C");
  DoTestMV2(a.View(),m1,cvx.col(0),cv0,label+" Step R,C");
  DoTestMV3(a.View(),cvx.col(0),cv.View(),cv0,label+" Step R,C");
  DoTestMV3(a.View(),cvx.col(0),cvx.col(1),cv0,label+" Step R,C");
  DoTestMV1(ca.View(),cm1,vx.col(0),label+" Step C,R");
  DoTestMV1(ca.View(),cm1,cvx.col(0),label+" Step C,C");
  DoTestMV2(ca.View(),cm1,cvx.col(0),cv0,label+" Step C,C");
  DoTestMV3(ca.View(),cvx.col(0),cv.View(),cv0,label+" Step C,C");
  DoTestMV3(ca.View(),cvx.col(0),cvx.col(1),cv0,label+" Step C,C");

  DoTestMV1(a.View(),m1,vx.col(0).Reverse(),label+" Rev R,R");
  DoTestMV2(a.View(),m1,vx.col(0).Reverse(),v0,label+" Rev R,R");
  DoTestMV3(a.View(),vx.col(0).Reverse(),v.View(),v0,label+" Rev R,R");
  DoTestMV3(a.View(),vx.col(0).Reverse(),vx.col(1),v0,label+" Rev Step R,R");
  DoTestMV3(a.View(),vx.col(0),vx.col(1).Reverse(),v0,label+" Step Rev R,R");
  DoTestMV3(a.View(),vx.col(0).Reverse(),vx.col(1).Reverse(),v0,label+" Rev Rev R,R");
  DoTestMV1(a.View(),m1,cvx.col(0).Reverse(),label+" Rev R,C");
  DoTestMV2(a.View(),m1,cvx.col(0).Reverse(),cv0,label+" Rev R,C");
  DoTestMV3(a.View(),cvx.col(0).Reverse(),cv.View(),cv0,label+" Rev R,C");
  DoTestMV3(a.View(),cvx.col(0).Reverse(),cvx.col(1),cv0,label+" Rev Step R,C");
  DoTestMV3(a.View(),cvx.col(0),cvx.col(1).Reverse(),cv0,label+" Step Rev R,C");
  DoTestMV3(a.View(),cvx.col(0).Reverse(),cvx.col(1).Reverse(),cv0,label+" Rev Rev R,C");
  DoTestMV1(ca.View(),cm1,vx.col(0).Reverse(),label+" Rev C,R");
  DoTestMV1(ca.View(),cm1,cvx.col(0).Reverse(),label+" Rev C,C");
  DoTestMV2(ca.View(),cm1,cvx.col(0).Reverse(),cv0,label+" Rev C,C");
  DoTestMV3(ca.View(),cvx.col(0).Reverse(),cv.View(),cv0,label+" Rev C,C");
  DoTestMV3(ca.View(),cvx.col(0).Reverse(),cvx.col(1),cv0,label+" Rev Step C,C");
  DoTestMV3(ca.View(),cvx.col(0),cvx.col(1).Reverse(),cv0,label+" Step Rev C,C");
  DoTestMV3(ca.View(),cvx.col(0).Reverse(),cvx.col(1).Reverse(),cv0,label+" Rev Rev C,C");

  if (a.rowsize() != a.colsize()) {
    Vector<T> w0(a.rowsize(),T(4));
    Vector<complex<T> > cw0(a.rowsize(),complex<T>(4,-4));
    Vector<T> w = w0;
    Vector<complex<T> > cw = cw0;
    Matrix<T,RowMajor> wx(a.rowsize(),31);
    Matrix<complex<T>,RowMajor> cwx(a.rowsize(),31);
    wx.col(0) = w;
    cwx.col(0) = cw;
    cwx.col(1) = cw;

    DoTestMV1(a.View(),m1,w,label+" R,R");
    DoTestMV2(a.View(),m1,w.View(),w0,label+" R,R");
    DoTestMV1(a.View(),m1,cw,label+" R,C");
    DoTestMV2(a.View(),m1,cw.View(),cw0,label+" R,C");
    DoTestMV1(ca.View(),cm1,w,label+" C,R");
    DoTestMV1(ca.View(),cm1,cw,label+" C,C");
    DoTestMV2(ca.View(),cm1,cw.View(),cw0,label+" C,C");

    DoTestMV1(a.View(),m1,wx.col(0),label+" Step R,R");
    DoTestMV2(a.View(),m1,wx.col(0),w0,label+" Step R,R");
    DoTestMV3(a.View(),wx.col(0),w.View(),w0,label+" Step R,R");
    DoTestMV3(a.View(),wx.col(0),wx.col(1),w0,label+" Step R,R");
    DoTestMV1(a.View(),m1,cwx.col(0),label+" Step R,C");
    DoTestMV2(a.View(),m1,cwx.col(0),cw0,label+" Step R,C");
    DoTestMV3(a.View(),cwx.col(0),cw.View(),cw0,label+" Step R,C");
    DoTestMV3(a.View(),cwx.col(0),cwx.col(1),cw0,label+" Step R,C");
    DoTestMV1(ca.View(),cm1,wx.col(0),label+" Step C,R");
    DoTestMV1(ca.View(),cm1,cwx.col(0),label+" Step C,C");
    DoTestMV2(ca.View(),cm1,cwx.col(0),cw0,label+" Step C,C");
    DoTestMV3(ca.View(),cwx.col(0),cw.View(),cw0,label+" Step C,C");
    DoTestMV3(ca.View(),cwx.col(0),cwx.col(1),cw0,label+" Step C,C");

    DoTestMV1(a.View(),m1,wx.col(0).Reverse(),label+" Rev R,R");
    DoTestMV2(a.View(),m1,wx.col(0).Reverse(),w0,label+" Rev R,R");
    DoTestMV3(a.View(),wx.col(0).Reverse(),w.View(),w0,label+" Rev R,R");
    DoTestMV3(a.View(),wx.col(0).Reverse(),wx.col(1),w0,label+" Rev Step R,R");
    DoTestMV3(a.View(),wx.col(0),wx.col(1).Reverse(),w0,label+" Step Rev R,R");
    DoTestMV3(a.View(),wx.col(0).Reverse(),wx.col(1).Reverse(),w0,label+" Rev Rev R,R");
    DoTestMV1(a.View(),m1,cwx.col(0).Reverse(),label+" Rev R,C");
    DoTestMV2(a.View(),m1,cwx.col(0).Reverse(),cw0,label+" Rev R,C");
    DoTestMV3(a.View(),cwx.col(0).Reverse(),cw.View(),cw0,label+" Rev R,C");
    DoTestMV3(a.View(),cwx.col(0).Reverse(),cwx.col(1),cw0,label+" Rev Step R,C");
    DoTestMV3(a.View(),cwx.col(0),cwx.col(1).Reverse(),cw0,label+" Step Rev R,C");
    DoTestMV3(a.View(),cwx.col(0).Reverse(),cwx.col(1).Reverse(),cw0,label+" Rev Rev R,C");
    DoTestMV1(ca.View(),cm1,wx.col(0).Reverse(),label+" Rev C,R");
    DoTestMV1(ca.View(),cm1,cwx.col(0).Reverse(),label+" Rev C,C");
    DoTestMV2(ca.View(),cm1,cwx.col(0).Reverse(),cw0,label+" Rev C,C");
    DoTestMV3(ca.View(),cwx.col(0).Reverse(),cw.View(),cw0,label+" Rev C,C");
    DoTestMV3(ca.View(),cwx.col(0).Reverse(),cwx.col(1),cw0,label+" Rev Step C,C");
    DoTestMV3(ca.View(),cwx.col(0),cwx.col(1).Reverse(),cw0,label+" Step Rev C,C");
    DoTestMV3(ca.View(),cwx.col(0).Reverse(),cwx.col(1).Reverse(),cw0,label+" Rev Rev C,C");
  }

  DoTestMM1(a.View(),b,m1,m2,label+" R,R");
  DoTestMM2(a.View(),b,m1,m2,a0,label+" R,R");
  m1 = a = a0;
  DoTestMM1(a.View(),cb,m1,cm2,label+" R,C");
  DoTestMM1(ca.View(),b,cm1,m2,label+" C,R");
  DoTestMM2(ca.View(),b,cm1,m2,ca0,label+" C,R");
  cm1 = ca = ca0;
  DoTestMM1(ca.View(),cb,cm1,cm2,label+" C,C");
  DoTestMM2(ca.View(),cb,cm1,cm2,ca0,label+" C,C");
  cm1 = ca = ca0;
  

#ifdef SHOWSTARTDONE
  cout<<"Done Test"<<endl;
#endif
}

