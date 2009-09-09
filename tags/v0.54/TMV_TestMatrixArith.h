#include "TMV.h"

template <class SM1, class SM2> inline bool CanAdd(const SM1& a, const SM2& b)
{ return a.colsize() == b.colsize() && a.rowsize() == b.rowsize(); }

template <class SM1, class SM2> inline bool CanAddEq(const SM1& a, const SM2& b)
{ return CanAdd(a,b); }

template <class SM1, class T2> inline bool CanAddX(const SM1& a, const T2)
{ return a.IsSquare(); }

template <class SM1, class T2> inline bool CanAddEqX(const SM1& a, const T2 x)
{ return CanAddX(a,x); }

template <class SM1, class T2> inline bool CanMultX(const SM1&, const T2)
{ return true; }

template <class SM1, class T2> inline bool CanMultEqX(const SM1& a, const T2 x)
{ return CanMultX(a,x); }

template <class SM1, class SM2> inline bool CanMult(const SM1& a, const SM2& b)
{ return a.rowsize() == b.colsize(); }

template <class SM, class T> inline bool CanMult(const SM& m, const tmv::Vector<T>& v)
{ return m.rowsize() == v.size(); }

template <class SM, class T> inline bool CanMult(const tmv::Vector<T>& v, const SM& m)
{ return v.size() == m.colsize(); }

template <class SM1, class SM2> inline bool CanMultEq(const SM1& a, const SM2& b)
{ return CanMult(a,b) && b.IsSquare(); }

template <class SM1, class SM2> inline bool CanMultEq2(const SM1& a, const SM2& b)
{ return CanMult(a,b) && a.IsSquare(); }

template <class SM1> inline bool CanDoSV(const SM1&)
{ return true; }

template <class SM, class T, class V> inline void DoTestMV1a(
    const SM& a, const tmv::Matrix<T>& m, const V& v, std::string label)
{
  if (showstartdone)
    std::cout<<"Start MV1a"<<std::endl;
  Assert(Norm(a-m) <= EPS*Norm(m),label+" a != m");
  if (CanMult(m,ColVectorViewOf(v))) {
    Assert(Norm((a*v)-(m*v)) <= EPS*Norm(m)*Norm(v),label+" a*v");
  }
  if (CanMult(RowVectorViewOf(v),m)) {
    Assert(Norm((v*a)-(v*m)) <= EPS*Norm(v)*Norm(m),label+" v*a");
  }
  if (showstartdone)
    std::cout<<"Done MV1a"<<std::endl;
}

template <class SM, class T, class V, class T2> inline void DoTestMV1(
    const SM& a, const tmv::Matrix<T>& m, const V& v, const tmv::Vector<T2>& v0,
    std::string label)
{
  tmv::Matrix<T> mt = m.Transpose();
  tmv::Matrix<T> mc = m.Conjugate();
  tmv::Matrix<T> ma = m.Adjoint();

  v=v0;
  DoTestMV1a(a.View(),m,v,label);
  DoTestMV1a(Transpose(a),mt,v,label+" Trans");
  DoTestMV1a(Conjugate(a),mc,v,label+" Conj");
  DoTestMV1a(Adjoint(a),ma,v,label+" Adj");

  v.Zero();
  DoTestMV1a(a.View(),m,v,label);
  DoTestMV1a(Transpose(a),mt,v,label+" Trans");
  DoTestMV1a(Conjugate(a),mc,v,label+" Conj");
  DoTestMV1a(Adjoint(a),ma,v,label+" Adj");

  v = v0;
  tmv::VectorView<T2>(v).SubVector(0,v0.size()/2).Zero();
  DoTestMV1a(a.View(),m,v,label);
  DoTestMV1a(Transpose(a),mt,v,label+" Trans");
  DoTestMV1a(Conjugate(a),mc,v,label+" Conj");
  DoTestMV1a(Adjoint(a),ma,v,label+" Adj");

  v = v0;
  tmv::VectorView<T2>(v).SubVector(v0.size()/2,v0.size()).Zero();
  DoTestMV1a(a.View(),m,v,label);
  DoTestMV1a(Transpose(a),mt,v,label+" Trans");
  DoTestMV1a(Conjugate(a),mc,v,label+" Conj");
  DoTestMV1a(Adjoint(a),ma,v,label+" Adj");

  v = v0;
  tmv::VectorView<T2>(v).SubVector(0,v0.size()/4).Zero();
  tmv::VectorView<T2>(v).SubVector(3*v0.size()/4,v0.size()).Zero();
  DoTestMV1a(a.View(),m,v,label);
  DoTestMV1a(Transpose(a),mt,v,label+" Trans");
  DoTestMV1a(Conjugate(a),mc,v,label+" Conj");
  DoTestMV1a(Adjoint(a),ma,v,label+" Adj");

  if (v.size() > 1) {
    v = v0;
    tmv::VectorView<T2>(v).SubVector(0,1).Zero();
    DoTestMV1a(a.View(),m,v,label);
    DoTestMV1a(Transpose(a),mt,v,label+" Trans");
    DoTestMV1a(Conjugate(a),mc,v,label+" Conj");
    DoTestMV1a(Adjoint(a),ma,v,label+" Adj");

    v = v0;
    tmv::VectorView<T2>(v).SubVector(v0.size()-1,v0.size()).Zero();
    DoTestMV1a(a.View(),m,v,label);
    DoTestMV1a(Transpose(a),mt,v,label+" Trans");
    DoTestMV1a(Conjugate(a),mc,v,label+" Conj");
    DoTestMV1a(Adjoint(a),ma,v,label+" Adj");

    v = v0;
    tmv::VectorView<T2>(v).SubVector(0,1).Zero();
    tmv::VectorView<T2>(v).SubVector(v0.size()-1,v0.size()).Zero();
    DoTestMV1a(a.View(),m,v,label);
    DoTestMV1a(Transpose(a),mt,v,label+" Trans");
    DoTestMV1a(Conjugate(a),mc,v,label+" Conj");
    DoTestMV1a(Adjoint(a),ma,v,label+" Adj");
  }
}

template <class SM, class T, class V, class T2> inline void DoTestMV2a(
    const SM& a, const tmv::Matrix<T>& m, const V& v1, tmv::Vector<T2>& v2,
    std::string label)
{
  if (showstartdone)
    std::cout<<"Start MV2a"<<std::endl;
  tmv::Vector<T2> v0 = v1;
  if (CanMultEq(v2,m)) {
    v2 = v1 = v0;
    v1 = v0*a;
    v2 = tmv::Vector<T2>(v0*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=v0*a");
    v2 = v1 = v0;
    v1 = T(10)*v0*a;
    v2 = tmv::Vector<T2>(T(10)*v0*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=x*v0*a");
    v2 = v1 = v0;
    v1 = T2(10)*v0*a;
    v2 = tmv::Vector<T2>(T2(10)*v0*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=x2*v0*a");
    v2 = v1 = v0;
    v1 *= a;
    v2 = tmv::Vector<T2>(v2*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v*=a");
    v2 = v1 = v0;
    v1 *= T(10)*a;
    v2 = tmv::Vector<T2>(T(10)*v2*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v*=(x*a)");
    v2 = v1 = v0;
    v1 *= T2(10)*a;
    v2 = tmv::Vector<T2>(T2(10)*v2*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v*=(x2*a)");
    v2 = v1 = v0;
    v1 = v1*a;
    v2 = tmv::Vector<T2>(v2*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=v*a");
    v2 = v1 = v0;
    v1 = T(10)*v1*a;
    v2 = tmv::Vector<T2>(T(10)*v2*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=x*v*a");
    v2 = v1 = v0;
    v1 = T2(10)*v1*a;
    v2 = tmv::Vector<T2>(T2(10)*v2*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=x2*v*a");
    v2 = v1 = v0;
    v1 += v0*a;
    v2 += tmv::Vector<T2>(v0*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v+=v0*a");
    v2 = v1 = v0;
    v1 -= v0*a;
    v2 -= tmv::Vector<T2>(v0*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v-=v0*a");
    v2 = v1 = v0;
    v1 += T(10)*v0*a;
    v2 += tmv::Vector<T2>(T(10)*v0*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v+=x*v0*a");
    v2 = v1 = v0;
    v1 += T2(10)*v0*a;
    v2 += tmv::Vector<T2>(T2(10)*v0*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v+=x2*v0*a");
    v2 = v1 = v0;
    v1 += v1*a;
    v2 += tmv::Vector<T2>(v0*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v+=v*a");
    v2 = v1 = v0;
    v1 -= v1*a;
    v2 -= tmv::Vector<T2>(v0*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v-=v*a");
    v2 = v1 = v0;
    v1 += T(10)*v1*a;
    v2 += tmv::Vector<T2>(T(10)*v0*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v+=x*v*a");
    v2 = v1 = v0;
    v1 += T2(10)*v1*a;
    v2 += tmv::Vector<T2>(T2(10)*v0*m);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v+=x*v*a");
  }
  if (CanMultEq2(m,v2)) {
    v2 = v1 = v0;
    v1 = a*v0;
    v2 = tmv::Vector<T2>(m*v0);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=a*v0");
    v2 = v1 = v0;
    v1 = T(10)*a*v0;
    v2 = tmv::Vector<T2>(T(10)*m*v0);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=x*a*v0");
    v2 = v1 = v0;
    v1 = T2(10)*a*v0;
    v2 = tmv::Vector<T2>(T2(10)*m*v0);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=x2*a*v0");
    v2 = v1 = v0;
    v1 = a*v1;
    v2 = tmv::Vector<T2>(m*v2);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=a*v");
    v2 = v1 = v0;
    v1 = T(10)*a*v1;
    v2 = tmv::Vector<T2>(T(10)*m*v2);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=x*a*v");
    v2 = v1 = v0;
    v1 = T2(10)*a*v1;
    v2 = tmv::Vector<T2>(T2(10)*m*v2);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v=x2*a*v");
    v2 = v1 = v0;
    v1 += a*v0;
    v2 += tmv::Vector<T2>(m*v0);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v+=a*v0");
    v2 = v1 = v0;
    v1 -= a*v0;
    v2 -= tmv::Vector<T2>(m*v0);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v-=a*v0");
    v2 = v1 = v0;
    v1 += T(10)*a*v0;
    v2 += tmv::Vector<T2>(T(10)*m*v0);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v+=x*a*v0");
    v2 = v1 = v0;
    v1 += T2(10)*a*v0;
    v2 += tmv::Vector<T2>(T2(10)*m*v0);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v+=x2*a*v0");
    v2 = v1 = v0;
    v1 += a*v1;
    v2 += tmv::Vector<T2>(m*v0);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v+=a*v");
    v2 = v1 = v0;
    v1 -= a*v1;
    v2 -= tmv::Vector<T2>(m*v0);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v-=a*v");
    v2 = v1 = v0;
    v1 += T(10)*a*v1;
    v2 += tmv::Vector<T2>(T(10)*m*v0);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v+=x*a*v");
    v2 = v1 = v0;
    v1 += T2(10)*a*v1;
    v2 += tmv::Vector<T2>(T2(10)*m*v0);
    Assert(Norm(v1-v2) <= EPS*Norm(v0)*Norm(m),label+" v+=x2*a*v");
  }
  if (showstartdone)
    std::cout<<"Done MV2a"<<std::endl;
}

template <class SM, class T, class V, class T2> inline void DoTestMV2(
    const SM& a, const tmv::Matrix<T>& m, const V& v1, const tmv::Vector<T2>& v0,
    std::string label)
{
  tmv::Matrix<T> mt = m.Transpose();
  tmv::Matrix<T> mc = m.Conjugate();
  tmv::Matrix<T> ma = m.Adjoint();

  tmv::Vector<T2> v00 = v0;

  v1 = v00;
  tmv::Vector<T2> v2 = v00;
  DoTestMV2a(a.View(),m,v1,v2,label);
  v1 = v00;
  v2 = v00;
  DoTestMV2a(Transpose(a),mt,v1,v2,label+" Trans");
  v1 = v00;
  v2 = v00;
  DoTestMV2a(Conjugate(a),mc,v1,v2,label+" Conj");
  v1 = v00;
  v2 = v00;
  DoTestMV2a(Adjoint(a),ma,v1,v2,label+" Adj");

  v00.Zero();
  v1 = v00;
  v2 = v00;
  DoTestMV2a(a.View(),m,v1,v2,label);
  v1 = v00;
  v2 = v00;
  DoTestMV2a(Transpose(a),mt,v1,v2,label+" Trans");
  v1 = v00;
  v2 = v00;
  DoTestMV2a(Conjugate(a),mc,v1,v2,label+" Conj");
  v1 = v00;
  v2 = v00;
  DoTestMV2a(Adjoint(a),ma,v1,v2,label+" Adj");

  v00 = v0;
  v00.SubVector(0,v0.size()/2).Zero();
  v1 = v00;
  v2 = v00;
  DoTestMV2a(a.View(),m,v1,v2,label);
  v1 = v00;
  v2 = v00;
  DoTestMV2a(Transpose(a),mt,v1,v2,label+" Trans");
  v1 = v00;
  v2 = v00;
  DoTestMV2a(Conjugate(a),mc,v1,v2,label+" Conj");
  v1 = v00;
  v2 = v00;
  DoTestMV2a(Adjoint(a),ma,v1,v2,label+" Adj");

  v00 = v0;
  v00.SubVector(v0.size()/2,v0.size()).Zero();
  v1 = v00;
  v2 = v00;
  DoTestMV2a(a.View(),m,v1,v2,label);
  v1 = v00;
  v2 = v00;
  DoTestMV2a(Transpose(a),mt,v1,v2,label+" Trans");
  v1 = v00;
  v2 = v00;
  DoTestMV2a(Conjugate(a),mc,v1,v2,label+" Conj");
  v1 = v00;
  v2 = v00;
  DoTestMV2a(Adjoint(a),ma,v1,v2,label+" Adj");

  v00 = v0;
  v00.SubVector(0,v0.size()/4).Zero();
  v00.SubVector(3*v0.size()/4,v0.size()).Zero();
  v1 = v00;
  v2 = v00;
  DoTestMV2a(a.View(),m,v1,v2,label);
  v1 = v00;
  v2 = v00;
  DoTestMV2a(Transpose(a),mt,v1,v2,label+" Trans");
  v1 = v00;
  v2 = v00;
  DoTestMV2a(Conjugate(a),mc,v1,v2,label+" Conj");
  v1 = v00;
  v2 = v00;
  DoTestMV2a(Adjoint(a),ma,v1,v2,label+" Adj");
}

template <class SM, class V1, class V2, class T> inline void DoTestMV3a(
    const SM& a, const V1& v1, const V2& v2, tmv::Vector<T>& v3,
    std::string label)
{
  if (showstartdone)
    std::cout<<"Start MV3a"<<std::endl;
  if (CanMultEq(v3,a)) {
    v2 = v1*a;
    v3 = tmv::Vector<T>(v1)*tmv::Matrix<T>(a);
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2=v1*a");
    v2 = v3 = v1;
    v2 += v1*a;
    v3 += tmv::Vector<T>(tmv::Vector<T>(v1)*tmv::Matrix<T>(a));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2+=v1*a");
    v2 = v3 = v1;
    v2 -= v1*a;
    v3 -= tmv::Vector<T>(tmv::Vector<T>(v1)*tmv::Matrix<T>(a));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2-=v1*a");
    v2 = v3 = v1;
    v2 += T(10)*v1*a;
    v3 += tmv::Vector<T>(tmv::Vector<T>(T(10)*v1)*tmv::Matrix<T>(a));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2+=x*v1*a");
    v2 = v3 = v1;
    v2 = v2*a;
    v3 = tmv::Vector<T>(v1)*tmv::Matrix<T>(a);
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2=v2*a");
    v2 = v3 = v1;
    v2 += v2*a;
    v3 += tmv::Vector<T>(tmv::Vector<T>(v1)*tmv::Matrix<T>(a));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2+=v2*a");
    v2 = v3 = v1;
    v2 -= v2*a;
    v3 -= tmv::Vector<T>(tmv::Vector<T>(v1)*tmv::Matrix<T>(a));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2-=v2*a");
    v2 = v3 = v1;
    v2 += T(10)*v2*a;
    v3 += tmv::Vector<T>(tmv::Vector<T>(T(10)*v1)*tmv::Matrix<T>(a));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2+=x*v2*a");
    v2 = v3 = v1;
    v2 = v3*a;
    v3 = tmv::Vector<T>(v1)*tmv::Matrix<T>(a);
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2=v3*a");
    v2 = v3 = v1;
    v2 += v3*a;
    v3 += tmv::Vector<T>(tmv::Vector<T>(v1)*tmv::Matrix<T>(a));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2+=v3*a");
    v2 = v3 = v1;
    v2 -= v3*a;
    v3 -= tmv::Vector<T>(tmv::Vector<T>(v1)*tmv::Matrix<T>(a));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2-=v3*a");
    v2 = v3 = v1;
    v2 += T(10)*v3*a;
    v3 += tmv::Vector<T>(tmv::Vector<T>(T(10)*v1)*tmv::Matrix<T>(a));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2+=x*v3*a");
  }
  if (CanMultEq2(a,v3)) {
    v2 = a*v1;
    v3 = tmv::Matrix<T>(a)*tmv::Vector<T>(v1);
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2=a*v1");
    v2 = v3 = v1;
    v2 += a*v1;
    v3 += tmv::Vector<T>(tmv::Matrix<T>(a)*tmv::Vector<T>(v1));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2+=a*v1");
    v2 = v3 = v1;
    v2 -= a*v1;
    v3 -= tmv::Vector<T>(tmv::Matrix<T>(a)*tmv::Vector<T>(v1));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2-=a*v1");
    v2 = v3 = v1;
    v2 += T(10)*a*v1;
    v3 += tmv::Vector<T>(tmv::Matrix<T>(a)*tmv::Vector<T>(T(10)*v1));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2+=x*a*v1");
    v2 = v3 = v1;
    v2 = a*v2;
    v3 = tmv::Matrix<T>(a)*tmv::Vector<T>(v1);
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2=a*v2");
    v2 = v3 = v1;
    v2 += a*v2;
    v3 += tmv::Vector<T>(tmv::Matrix<T>(a)*tmv::Vector<T>(v1));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2+=a*v2");
    v2 = v3 = v1;
    v2 -= a*v2;
    v3 -= tmv::Vector<T>(tmv::Matrix<T>(a)*tmv::Vector<T>(v1));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2-=a*v2");
    v2 = v3 = v1;
    v2 += T(10)*a*v2;
    v3 += tmv::Vector<T>(tmv::Matrix<T>(a)*tmv::Vector<T>(T(10)*v1));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2+=x*a*v2");
    v2 = v3 = v1;
    v2 = a*v3;
    v3 = tmv::Matrix<T>(a)*tmv::Vector<T>(v1);
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2=a*v3");
    v2 = v3 = v1;
    v2 += a*v3;
    v3 += tmv::Vector<T>(tmv::Matrix<T>(a)*tmv::Vector<T>(v1));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2+=a*v3");
    v2 = v3 = v1;
    v2 -= a*v3;
    v3 -= tmv::Vector<T>(tmv::Matrix<T>(a)*tmv::Vector<T>(v1));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2-=a*v3");
    v2 = v3 = v1;
    v2 += T(10)*a*v3;
    v3 += tmv::Vector<T>(tmv::Matrix<T>(a)*tmv::Vector<T>(T(10)*v1));
    Assert(Norm(v2-v3) <= EPS*Norm(v1)*Norm(a),label+" v2+=x*a*v3");
  }
  if (showstartdone)
    std::cout<<"Done MV3a"<<std::endl;
}

template <class SM, class V1, class V2, class T2> inline void DoTestMV3(
    const SM& a, const V1& v1, const V2& v2, const tmv::Vector<T2>& v0,
    std::string label)
{
  if (showstartdone)
    std::cout<<"Start MV3"<<std::endl;
  v1 = v0;
  v2 = v0;
  tmv::Vector<T2> v3 = v0;

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
  if (showstartdone)
    std::cout<<"Done MV3"<<std::endl;
}

template <class SM1, class SM2, class T, class T2> inline void DoTestMM1a(
    const SM1& a, const SM2& b, const tmv::Matrix<T>& m1, const tmv::Matrix<T2>& m2,
    std::string label)
{
  if (showstartdone)
    std::cout<<"Start MM1a"<<std::endl;

//#define XXDEBUG

#ifdef XXDEBUG
  std::cout<<"a = "<<tmv::Type(a)<<" = "<<a<<std::endl;
  std::cout<<"b = "<<tmv::Type(b)<<" = "<<b<<std::endl;
  std::cout<<"m1 = "<<tmv::Type(m1)<<" = "<<m1<<std::endl;
  std::cout<<"m2 = "<<tmv::Type(m2)<<" = "<<m2<<std::endl;
  std::cout<<"Norm(a-m1) = "<<Norm(a-m1)<<std::endl;
  std::cout<<"Norm(b-m2) = "<<Norm(b-m2)<<std::endl;
#endif

  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a != m1 (MM1a)");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b != m2 (MM1a)");

  if (CanAdd(a,b)) {
#ifdef XXDEBUG
    std::cout<<"a-m2 = "<<a-m2<<std::endl;
    std::cout<<"m1-m2 = "<<m1-m2<<std::endl;
    std::cout<<"m1-b = "<<m1-b<<std::endl;
    std::cout<<"m1-m2 = "<<m1-m2<<std::endl;
    std::cout<<"a-b = "<<a-b<<std::endl;
    std::cout<<"m1-m2 = "<<m1-m2<<std::endl;
#endif
    Assert(Norm((a-m2)-(m1-m2)) <= EPS*(Norm(m1)+Norm(m2)),label+" a-m");
    Assert(Norm((m1-b)-(m1-m2)) <= EPS*(Norm(m1)+Norm(m2)),label+" m-b");
    Assert(Norm((a-b)-(m1-m2)) <= EPS*(Norm(m1)+Norm(m2)),label+" a-b");
    Assert(Norm((a+m2)-(m1+m2)) <= EPS*(Norm(m1)+Norm(m2)),label+" a+m");
    Assert(Norm((m1+b)-(m1+m2)) <= EPS*(Norm(m1)+Norm(m2)),label+" m+b");
    Assert(Norm((a+b)-(m1+m2)) <= EPS*(Norm(m1)+Norm(m2)),label+" a+b");
  }
  if (CanMult(a,b)) {
#ifdef XXDEBUG
    std::cout<<"m1*b = "<<m1*b<<std::endl;
    std::cout<<"m1*m2 = "<<m1*m2<<std::endl;
    std::cout<<"a*m2 = "<<a*m2<<std::endl;
    std::cout<<"m1*m2 = "<<m1*m2<<std::endl;
    std::cout<<"a*b = "<<a*b<<std::endl;
    std::cout<<"m1*m2 = "<<m1*m2<<std::endl;
    std::cout<<"a*b-m1*m2 = "<<a*b-m1*m2<<std::endl;
    std::cout<<"Norm(a*b-m1*m2) = "<<Norm(a*b-m1*m2)<<std::endl;
    std::cout<<"EPS*Norm(m1*m2) = "<<EPS*Norm(m1*m2)<<std::endl;
#endif
    Assert(Norm((m1*b)-(m1*m2)) <= EPS*Norm(m1)*Norm(m2),label+" m*b");
    Assert(Norm((a*m2)-(m1*m2)) <= EPS*Norm(m1)*Norm(m2),label+" a*m");
    Assert(Norm((a*b)-(m1*m2)) <= EPS*Norm(m1)*Norm(m2),label+" a*b");
  }

#ifdef XXDEBUG
#undef XXDEBUG
#endif

  if (showstartdone)
    std::cout<<"Done MM1a"<<std::endl;
}

template <class SM1, class SM2, class T, class T2> inline void DoTestMM1(
    const SM1& a, const SM2& b, const tmv::Matrix<T>& m1,
    const tmv::Matrix<T2>& m2, std::string label)
{
  tmv::Matrix<T> m1t = m1.Transpose();
  tmv::Matrix<T> m1c = m1.Conjugate();
  tmv::Matrix<T> m1a = m1.Adjoint();
  tmv::Matrix<T2> m2t = m2.Transpose();
  tmv::Matrix<T2> m2c = m2.Conjugate();
  tmv::Matrix<T2> m2a = m2.Adjoint();

  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a != m1 (MM1)");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b != m2 (MM1)");

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

template <class SM1, class SM2, class T, class T2> inline void DoTestMM2a(
    const SM1& a, const SM2& b, tmv::Matrix<T>& m1, const tmv::Matrix<T2>& m2,
    std::string label)
{
  if (showstartdone)
    std::cout<<"Start MM2a"<<std::endl;

//#define XXDEBUG

#ifdef XXDEBUG
  std::cout<<"a = "<<tmv::Type(a)<<"  "<<a<<std::endl;
  std::cout<<"b = "<<tmv::Type(b)<<"  "<<b<<std::endl;
  std::cout<<"m1 = "<<tmv::Type(m1)<<"  "<<m1<<std::endl;
  std::cout<<"m2 = "<<tmv::Type(m2)<<"  "<<m2<<std::endl;
#endif
  tmv::Matrix<T> m3 = m1;
  tmv::Matrix<T> m4 = m1;

  if (CanAddEq(m3,b)) {
    double normm = Norm(m1)+Norm(m2);
    m3 += b;
    m4 = tmv::Matrix<T>(m4+m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m += b");
    m3 = m4 = m1;
    m3 -= b;
    m4 = tmv::Matrix<T>(m4-m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m -= b");
    m3 = m4 = m1;
    m3 = m3 + b;
    m4 = tmv::Matrix<T>(m4+m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m = m+b");
    m3 = m4 = m1;
    m3 = m3 - b;
    m4 = tmv::Matrix<T>(m4-m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m = m-b");
    m4 = m3;
    m3 = b + m3;
    m4 = tmv::Matrix<T>(m2+m4);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m = b+m");
    m3 = m4 = m1;
    m3 = b - m3;
    m4 = tmv::Matrix<T>(m2-m4);
    Assert(Norm(m3-m4) <= EPS*Norm(m3),label+" m = b-m");
  }
  if (CanAddEq(a,b)) {
    m1 = a;
    double normm = Norm(m1)+Norm(m2);
    a += b;
    m1 = tmv::Matrix<T>(m1+b);
    Assert(Norm(a-m1) <= EPS*normm,label+" a += b");
    m1 = a;
    normm = Norm(m1)+Norm(m2);
    a -= b;
    m1 = tmv::Matrix<T>(m1-b);
    Assert(Norm(a-m1) <= EPS*normm,label+" a -= b");
    m1 = a;
    normm = Norm(m1)+Norm(m2);
    a = a+b;
    m1 = tmv::Matrix<T>(m1+b);
    Assert(Norm(a-m1) <= EPS*normm,label+" a = a+b");
    m1 = a;
    normm = Norm(m1)+Norm(m2);
    a = a-b;
    m1 = tmv::Matrix<T>(m1-b);
    Assert(Norm(a-m1) <= EPS*normm,label+" a = a-b");
    m1 = a;
    normm = Norm(m1)+Norm(m2);
    a = b+a; 
    m1 = tmv::Matrix<T>(m2+m1);
    Assert(Norm(a-m1) <= EPS*normm,label+" a = b+a");
    m1 = a;
    normm = Norm(m1)+Norm(m2);
    a = b-a;
    m1 = tmv::Matrix<T>(m2-m1);
    Assert(Norm(a-m1) <= EPS*normm,label+" a = b-a");
  }
  if (CanMultEq(m1,b)) {
    double normm = Norm(m1)*Norm(m2);
    m3 = m4 = m1;
    m3 *= b;
    m4 = tmv::Matrix<T>(m4*m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m *= b");
    m3 = m4 = m1;
    m3 = m3*b;
    m4 = tmv::Matrix<T>(m4*m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m = m*b");
    m3 = m4 = m1;
    m3 += m1*b;
    m4 += tmv::Matrix<T>(m1*m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m += m0*b");
    m3 = m4 = m1;
    m3 -= m1*b;
    m4 -= tmv::Matrix<T>(m1*m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m -= m0*b");
    m3 = m4 = m1;
    m3 += T(10)*m1*b;
    m4 += tmv::Matrix<T>(T(10)*m1*m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m += x*m0*b");
    m3 = m4 = m1;
    m3 += m3*b;
    m4 += tmv::Matrix<T>(m4*m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m += m*b");
    m3 = m4 = m1;
    m3 -= m3*b;
    m4 -= tmv::Matrix<T>(m4*m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m -= m*b");
    m3 = m4 = m1;
    m3 += T(10)*m3*b;
    m4 += tmv::Matrix<T>(T(10)*m4*m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m += x*m*b");
    m3 = m4 = m1;
    m3 += a*b;
    m4 += tmv::Matrix<T>(m4*m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m += a*b");
    m3 = m4 = m1;
    m3 -= a*b;
    m4 -= tmv::Matrix<T>(m4*m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m -= a*b");
    m3 = m4 = m1;
    m3 += T(10)*a*b;
    m4 += tmv::Matrix<T>(T(10)*m4*m2);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m += x*a*b");
  }
  if (CanMultEq2(b,m1)) {
    m3 = m4 = m1;
    double normm = Norm(m1)*Norm(m2);
    m3 = b * m3;
    m4 = tmv::Matrix<T>(m2*m1);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m = b*m");
    m3 = m4 = m1;
    m3 += b * m3;
    m4 += tmv::Matrix<T>(m2*m1);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m += b*m");
    m3 = m4 = m1;
    m3 -= b * m3;
    m4 -= tmv::Matrix<T>(m2*m1);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m -= b*m");
    m3 = m4 = m1;
    m3 += T(10)*b * m3;
    m4 += tmv::Matrix<T>(T(10)*m2*m1);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m += x*b*m");
    m3 = m4 = m1;
    m3 += b * a;
    m4 += tmv::Matrix<T>(m2*m1);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m += b*a");
    m3 = m4 = m1;
    m3 -= b * a;
    m4 -= tmv::Matrix<T>(m2*m1);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m -= b*a");
    m3 = m4 = m1;
    m3 += T(10)*b*a;
    m4 += tmv::Matrix<T>(T(10)*m2*m1);
    Assert(Norm(m3-m4) <= EPS*normm,label+" m += x*b*a");
  }
  if (CanMultEq(a,b)) {
    m1 = a;
    double normm = Norm(m1)*Norm(m2);
    a *= b;
    m1 = tmv::Matrix<T>(m1*m2);
    Assert(Norm(a-m1) <= EPS*normm,label+" a *= b");
    a /= Norm1(b);
    m1 = a;
    normm = Norm(m1)*Norm(m2);
    a = a*b;
    m1 = tmv::Matrix<T>(m1*m2);
    Assert(Norm(a-m1) <= EPS*normm,label+" a = a*b");
    a /= Norm1(b);
    m1 = a;
  }
  if (CanMultEq2(b,a)) {
    m1 = a;
    double normm = Norm(m1)*Norm(m2);
    a = b * a;
    m1 = tmv::Matrix<T>(m2*m1);
    Assert(Norm(a-m1) <= EPS*normm,label+" a = b*a");
    a /= Norm1(b);
    m1 = a;
  }

#ifdef XXDEBUG
#undef XXDEBUG
#endif

  if (showstartdone)
    std::cout<<"Done MM2a"<<std::endl;
}

template <class SM1, class SM2, class T, class T2, class BaseSM1> inline void DoTestMM2(
    const SM1& a, const SM2& b, tmv::Matrix<T>& m1, 
    tmv::Matrix<T2>& m2, const BaseSM1& basea, std::string label)
{
  m1 = a = basea;
  DoTestMM2a(a.View(),b.View(),m1,m2,label);
  m1 = a = basea;
  tmv::Matrix<T> m1t = m1.Transpose();
  DoTestMM2a(Transpose(a),b.View(),m1t,m2,label+" TransA");
  m1 = a = basea;
  tmv::Matrix<T> m1c = m1.Conjugate();
  DoTestMM2a(Conjugate(a),b.View(),m1c,m2,label+" ConjA");
  m1 = a = basea;
  tmv::Matrix<T> m1a = m1.Adjoint();
  DoTestMM2a(Adjoint(a),b.View(),m1a,m2,label+" AdjA");
  m1 = a = basea;

  tmv::Matrix<T2> m2t = m2.Transpose();
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

  tmv::Matrix<T2> m2c = m2.Conjugate();
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

  tmv::Matrix<T2> m2a = m2.Adjoint();
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

template <class SM1, class T, class T2> inline void DoTestMX1a(
    const SM1& a, const tmv::Matrix<T>& m1, T2 x, std::string label)
{
  if (showstartdone)
    std::cout<<"Start MX1a"<<std::endl;
  double normm = Norm(m1);
  if (CanAddX(a,x)) {
    Assert(Norm(a-m1) <= EPS*normm,label+" a != m1");
    Assert(Norm((x-a)-(x-m1)) <= EPS*normm,label+" x-a");
    Assert(Norm((a-x)-(m1-x)) <= EPS*normm,label+" a-x");
    Assert(Norm((x+a)-(x+m1)) <= EPS*normm,label+" x+a");
    Assert(Norm((a+x)-(m1+x)) <= EPS*normm,label+" a+x");
  }
  if (CanMultX(a,x)) {
    Assert(Norm((x*a)-(x*m1)) <= EPS*normm*std::abs(x),label+" x*a");
    Assert(Norm((a*x)-(m1*x)) <= EPS*normm*std::abs(x),label+" a*x");
    Assert(Norm((a/x)-(m1/x)) <= EPS*normm*std::abs(x),label+" a/x");
  }
  if (showstartdone)
    std::cout<<"Done MX1a"<<std::endl;
}

template <class SM1, class T, class T2> inline void DoTestMX1(
    const SM1& a, tmv::Matrix<T>& m1, T2 x, std::string label)
{
  DoTestMX1a(a.View(),m1,x,label);
  tmv::Matrix<T> m1t = m1.Transpose();
  DoTestMX1a(Transpose(a),m1t,x,label+" Trans");
  tmv::Matrix<T> m1c = m1.Conjugate();
  DoTestMX1a(Conjugate(a),m1c,x,label+" Conj");
  tmv::Matrix<T> m1a = m1.Adjoint();
  DoTestMX1a(Adjoint(a),m1a,x,label+" Adj");
}

template <class SM1, class T, class T2> inline void DoTestMX2a(
    const SM1& a, tmv::Matrix<T>& m1, T2 x, std::string label)
{
  if (showstartdone)
    std::cout<<"Start MX2a"<<std::endl;

//#define XXDEBUG

  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = m1");

  if (CanAddEqX(a,x)) {
    Assert(Norm((a+=x)-(m1+=x)) <= EPS*Norm(m1),label+" a += x");
    Assert(Norm((a-=x)-(m1-=x)) <= EPS*Norm(m1),label+" a -= x");
    a = a+x; 
    m1 = tmv::Matrix<T>(m1+x);
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = a+x");
    m1 = a;
    a = a-x;
    m1 = tmv::Matrix<T>(m1-x);
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = a-x");
    m1 = a;
    a = x+a;
    m1 = tmv::Matrix<T>(x+m1);
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = x+a");
    m1 = a;
    a = x-a; 
    m1 = tmv::Matrix<T>(x-m1);
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = x-a");
  }
  if (CanMultEqX(a,x)) {
    m1 = a;
    a*=x;
    m1*=x;
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a *= x");
    Assert(Norm((a*=x)-(m1*=x)) <= EPS*Norm(m1),label+" a *= x");
    m1 = a;
    Assert(Norm((a/=x)-(m1/=x)) <= EPS*Norm(m1),label+" a /= x");
    m1 = a;
    a = a*x;
    m1 = tmv::Matrix<T>(m1*x);
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = a*x");
    m1 = a;
    a = a/x;
    m1 = tmv::Matrix<T>(m1/x);
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = a/x");
    m1 = a;
    a = x*a; 
    m1 = tmv::Matrix<T>(x*m1);
    Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = x*a");
    m1 = a;
    a /= x;
    m1 = tmv::Matrix<T>(m1/x);
  }
#ifdef XXDEBUG
#undef XXDEBUG
#endif

  if (showstartdone)
    std::cout<<"Done MX2a"<<std::endl;
}

template <class SM1, class T, class T2, class BaseSM1> inline void DoTestMX2(
    const SM1& a, tmv::Matrix<T>& m1, T2 x,
    const BaseSM1& basea, std::string label)
{
  m1 = a = basea;
  DoTestMX2a(a.View(),m1,x,label);
  m1 = a = basea;
  tmv::Matrix<T> m1t = m1.Transpose();
  DoTestMX2a(Transpose(a),m1t,x,label+" Trans");
  m1 = a = basea;
  tmv::Matrix<T> m1c = m1.Conjugate();
  DoTestMX2a(Conjugate(a),m1c,x,label+" Conj");
  m1 = a = basea;
  tmv::Matrix<T> m1a = m1.Adjoint();
  DoTestMX2a(Adjoint(a),m1a,x,label+" Adj");
}

template <class SM1, class T> inline void DoTestMa(
    const SM1& a, const tmv::Matrix<T>& m, std::string label)
{
  if (showstartdone)
    std::cout<<"Start Ma"<<std::endl;

//#define XXDEBUG

#ifdef XXDEBUG
  std::cout<<"a = "<<tmv::Type(a)<<" = "<<a<<std::endl;
  std::cout<<"m = "<<tmv::Type(m)<<" = "<<m<<std::endl;
  std::cout<<"a-m = "<<a-m<<std::endl;
  std::cout<<"Norm(a-m) = "<<Norm(a-m)<<std::endl;
  std::cout<<"Trace(a) = "<<Trace(a)<<std::endl;
  std::cout<<"NormF(a) = "<<NormF(a)<<std::endl;
  std::cout<<"Norm(a) = "<<Norm(a)<<std::endl;
  std::cout<<"Norm1(a) = "<<Norm1(a)<<std::endl;
  std::cout<<"Norm1(m) = "<<Norm1(m)<<std::endl;
#endif
  Assert(Norm(a-m) <= EPS*Norm(m),label+" a != m");
  Assert(std::abs(Trace(a)-Trace(m)) <= EPS*std::abs(Trace(m)),label+" Trace");
  Assert(std::abs(NormF(a)-NormF(m)) <= EPS*NormF(m),label+" NormF");
  Assert(std::abs(Norm(a)-Norm(m)) <= EPS*Norm(m),label+" Norm");
  Assert(std::abs(Norm1(a)-Norm1(m)) <= EPS*Norm1(m),label+" Norm1");
  if (donorm2 && CanDoSV(a)) {
    a.DivideUsing(tmv::SV);
    m.DivideUsing(tmv::SV);
#ifdef XXDEBUG
    std::cout<<"Norm2(a) = "<<Norm2(a)<<std::endl;
    std::cout<<"Norm2(m) = "<<Norm2(m)<<std::endl;
#endif
    Assert(std::abs(Norm2(a)-Norm2(m)) <= EPS*Norm2(m),label+" Norm2");
  }
#ifdef XXDEBUG
  std::cout<<"NormInf(a) = "<<NormInf(a)<<std::endl;
  std::cout<<"NormInf(m) = "<<NormInf(m)<<std::endl;
  std::cout<<"abs(diff) = "<<std::abs(NormInf(a)-NormInf(m))<<std::endl;
  std::cout<<"eps*norminf = "<<EPS*NormInf(m)<<std::endl;
  std::cout<<"Norm(aT-mT) = "<<Norm(Transpose(a)-Transpose(m))<<std::endl;
  std::cout<<"Conjugate(a) = "<<Conjugate(a)<<std::endl;
  std::cout<<"Conjugate(m) = "<<Conjugate(m)<<std::endl;
  std::cout<<"a*-m* = "<<Conjugate(a)-Conjugate(m)<<std::endl;
  std::cout<<"Conjugate(a).diag = "<<Conjugate(a).diag()<<std::endl;
  std::cout<<"Conjugate(m).diag = "<<Conjugate(m).diag()<<std::endl;
  std::cout<<"Norm(a*-m*) = "<<Norm(Conjugate(a)-Conjugate(m))<<std::endl;
  std::cout<<"Norm(at-mt) = "<<Norm(Adjoint(a)-Adjoint(m))<<std::endl;
#endif
  Assert(std::abs(NormInf(a)-NormInf(m)) <= EPS*NormInf(m),label+" NormInf");
  Assert(Norm(Transpose(a)-Transpose(m)) <= EPS*Norm(m),label+" Transpose");
  Assert(Norm(Conjugate(a)-Conjugate(m)) <= EPS*Norm(m),label+" Conjugate");
  Assert(Norm(Adjoint(a)-Adjoint(m)) <= EPS*Norm(m),label+" Adjoint");

#ifdef XXDEBUG
#undef XXDEBUG
#endif

  if (showstartdone)
    std::cout<<"Done Ma"<<std::endl;
}

template <class SM1, class T> inline void DoTestM(
    const SM1& a, const tmv::Matrix<T>& m, std::string label)
{
  DoTestMa(a.View(),m,label);
  tmv::Matrix<T> mt = m.Transpose();
  DoTestMa(Transpose(a),mt,label+" Trans");
  tmv::Matrix<T> mc = m.Conjugate();
  DoTestMa(Conjugate(a),mc,label+" Conj");
  tmv::Matrix<T> ma = m.Adjoint();
  DoTestMa(Adjoint(a),ma,label+" Adj");
}

template <class T, class BaseM, class BaseCM, class SM1, class SM2, class CSM1, class CSM2> inline void TestMatrixArith(
    const SM1& a, const SM2& b, const CSM1& ca, const CSM2& cb, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start TestMatrixArith"<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" "<<b<<std::endl;
    std::cout<<"ca = "<<tmv::Type(ca)<<" "<<ca<<std::endl;
    std::cout<<"cb = "<<tmv::Type(cb)<<" "<<cb<<std::endl;
  }

  const BaseM a0 = a;
  const BaseCM ca0 = ca;

  tmv::Matrix<T> m1 = a;
  tmv::Matrix<T> m2 = b;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a != m1 (basic)");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b != m2 (basic)");
  tmv::Matrix<std::complex<T> > cm1 = ca;
  tmv::Matrix<std::complex<T> > cm2 = cb;

#ifdef XTEST
  T x = 12;
  std::complex<T> z(9,-2);

  tmv::Vector<T> v0(a.colsize(),T(4));
  tmv::Vector<std::complex<T> > cv0(a.colsize(),std::complex<T>(4,-4));

  tmv::Vector<T> v = v0;
  tmv::Vector<std::complex<T> > cv = cv0;

  tmv::Matrix<T> vx(a.colsize(),5);
  tmv::Matrix<std::complex<T> > cvx(a.colsize(),5);
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

  DoTestMV1(a.View(),m1,v.View(),v0,label+" R,R");
  DoTestMV2(a.View(),m1,v.View(),v0,label+" R,R");
  DoTestMV1(a.View(),m1,cv.View(),cv0,label+" R,C");
  DoTestMV2(a.View(),m1,cv.View(),cv0,label+" R,C");
  DoTestMV1(ca.View(),cm1,v.View(),v0,label+" C,R");
  DoTestMV1(ca.View(),cm1,cv.View(),cv0,label+" C,C");
  DoTestMV2(ca.View(),cm1,cv.View(),cv0,label+" C,C");

  DoTestMV1(a.View(),m1,vx.col(0),v0,label+" Step R,R");
  DoTestMV2(a.View(),m1,vx.col(0),v0,label+" Step R,R");
  DoTestMV3(a.View(),vx.col(0),v.View(),v0,label+" Step R,R");
  DoTestMV3(a.View(),vx.col(0),vx.col(1),v0,label+" Step R,R");
  DoTestMV1(a.View(),m1,cvx.col(0),cv0,label+" Step R,C");
  DoTestMV2(a.View(),m1,cvx.col(0),cv0,label+" Step R,C");
  DoTestMV3(a.View(),cvx.col(0),cv.View(),cv0,label+" Step R,C");
  DoTestMV3(a.View(),cvx.col(0),cvx.col(1),cv0,label+" Step R,C");
  DoTestMV1(ca.View(),cm1,vx.col(0),v0,label+" Step C,R");
  DoTestMV1(ca.View(),cm1,cvx.col(0),cv0,label+" Step C,C");
  DoTestMV2(ca.View(),cm1,cvx.col(0),cv0,label+" Step C,C");
  DoTestMV3(ca.View(),cvx.col(0),cv.View(),cv0,label+" Step C,C");
  DoTestMV3(ca.View(),cvx.col(0),cvx.col(1),cv0,label+" Step C,C");

  DoTestMV1(a.View(),m1,vx.col(0).Reverse(),v0,label+" Rev R,R");
  DoTestMV2(a.View(),m1,vx.col(0).Reverse(),v0,label+" Rev R,R");
  DoTestMV3(a.View(),vx.col(0).Reverse(),v.View(),v0,label+" Rev R,R");
  DoTestMV3(a.View(),vx.col(0).Reverse(),vx.col(1),v0,label+" Rev Step R,R");
  DoTestMV3(a.View(),vx.col(0),vx.col(1).Reverse(),v0,label+" Step Rev R,R");
  DoTestMV3(a.View(),vx.col(0).Reverse(),vx.col(1).Reverse(),v0,label+" Rev Rev R,R");
  DoTestMV1(a.View(),m1,cvx.col(0).Reverse(),cv0,label+" Rev R,C");
  DoTestMV2(a.View(),m1,cvx.col(0).Reverse(),cv0,label+" Rev R,C");
  DoTestMV3(a.View(),cvx.col(0).Reverse(),cv.View(),cv0,label+" Rev R,C");
  DoTestMV3(a.View(),cvx.col(0).Reverse(),cvx.col(1),cv0,label+" Rev Step R,C");
  DoTestMV3(a.View(),cvx.col(0),cvx.col(1).Reverse(),cv0,label+" Step Rev R,C");
  DoTestMV3(a.View(),cvx.col(0).Reverse(),cvx.col(1).Reverse(),cv0,label+" Rev Rev R,C");
  DoTestMV1(ca.View(),cm1,vx.col(0).Reverse(),v0,label+" Rev C,R");
  DoTestMV1(ca.View(),cm1,cvx.col(0).Reverse(),cv0,label+" Rev C,C");
  DoTestMV2(ca.View(),cm1,cvx.col(0).Reverse(),cv0,label+" Rev C,C");
  DoTestMV3(ca.View(),cvx.col(0).Reverse(),cv.View(),cv0,label+" Rev C,C");
  DoTestMV3(ca.View(),cvx.col(0).Reverse(),cvx.col(1),cv0,label+" Rev Step C,C");
  DoTestMV3(ca.View(),cvx.col(0),cvx.col(1).Reverse(),cv0,label+" Step Rev C,C");
  DoTestMV3(ca.View(),cvx.col(0).Reverse(),cvx.col(1).Reverse(),cv0,label+" Rev Rev C,C");

  if (a.rowsize() != a.colsize()) {
    tmv::Vector<T> w0(a.rowsize(),T(4));
    tmv::Vector<std::complex<T> > cw0(a.rowsize(),std::complex<T>(4,-4));
    tmv::Vector<T> w = w0;
    tmv::Vector<std::complex<T> > cw = cw0;
    tmv::Matrix<T> wx(a.rowsize(),31);
    tmv::Matrix<std::complex<T> > cwx(a.rowsize(),31);
    wx.col(0) = w;
    cwx.col(0) = cw;
    cwx.col(1) = cw;

    DoTestMV1(a.View(),m1,w.View(),w0,label+" R,R");
    DoTestMV2(a.View(),m1,w.View(),w0,label+" R,R");
    DoTestMV1(a.View(),m1,cw.View(),cw0,label+" R,C");
    DoTestMV2(a.View(),m1,cw.View(),cw0,label+" R,C");
    DoTestMV1(ca.View(),cm1,w.View(),w0,label+" C,R");
    DoTestMV1(ca.View(),cm1,cw.View(),cw0,label+" C,C");
    DoTestMV2(ca.View(),cm1,cw.View(),cw0,label+" C,C");

    DoTestMV1(a.View(),m1,wx.col(0),w0,label+" Step R,R");
    DoTestMV2(a.View(),m1,wx.col(0),w0,label+" Step R,R");
    DoTestMV3(a.View(),wx.col(0),w.View(),w0,label+" Step R,R");
    DoTestMV3(a.View(),wx.col(0),wx.col(1),w0,label+" Step R,R");
    DoTestMV1(a.View(),m1,cwx.col(0),cw0,label+" Step R,C");
    DoTestMV2(a.View(),m1,cwx.col(0),cw0,label+" Step R,C");
    DoTestMV3(a.View(),cwx.col(0),cw.View(),cw0,label+" Step R,C");
    DoTestMV3(a.View(),cwx.col(0),cwx.col(1),cw0,label+" Step R,C");
    DoTestMV1(ca.View(),cm1,wx.col(0),w0,label+" Step C,R");
    DoTestMV1(ca.View(),cm1,cwx.col(0),cw0,label+" Step C,C");
    DoTestMV2(ca.View(),cm1,cwx.col(0),cw0,label+" Step C,C");
    DoTestMV3(ca.View(),cwx.col(0),cw.View(),cw0,label+" Step C,C");
    DoTestMV3(ca.View(),cwx.col(0),cwx.col(1),cw0,label+" Step C,C");

    DoTestMV1(a.View(),m1,wx.col(0).Reverse(),w0,label+" Rev R,R");
    DoTestMV2(a.View(),m1,wx.col(0).Reverse(),w0,label+" Rev R,R");
    DoTestMV3(a.View(),wx.col(0).Reverse(),w.View(),w0,label+" Rev R,R");
    DoTestMV3(a.View(),wx.col(0).Reverse(),wx.col(1),w0,label+" Rev Step R,R");
    DoTestMV3(a.View(),wx.col(0),wx.col(1).Reverse(),w0,label+" Step Rev R,R");
    DoTestMV3(a.View(),wx.col(0).Reverse(),wx.col(1).Reverse(),w0,label+" Rev Rev R,R");
    DoTestMV1(a.View(),m1,cwx.col(0).Reverse(),cw0,label+" Rev R,C");
    DoTestMV2(a.View(),m1,cwx.col(0).Reverse(),cw0,label+" Rev R,C");
    DoTestMV3(a.View(),cwx.col(0).Reverse(),cw.View(),cw0,label+" Rev R,C");
    DoTestMV3(a.View(),cwx.col(0).Reverse(),cwx.col(1),cw0,label+" Rev Step R,C");
    DoTestMV3(a.View(),cwx.col(0),cwx.col(1).Reverse(),cw0,label+" Step Rev R,C");
    DoTestMV3(a.View(),cwx.col(0).Reverse(),cwx.col(1).Reverse(),cw0,label+" Rev Rev R,C");
    DoTestMV1(ca.View(),cm1,wx.col(0).Reverse(),w0,label+" Rev C,R");
    DoTestMV1(ca.View(),cm1,cwx.col(0).Reverse(),cw0,label+" Rev C,C");
    DoTestMV2(ca.View(),cm1,cwx.col(0).Reverse(),cw0,label+" Rev C,C");
    DoTestMV3(ca.View(),cwx.col(0).Reverse(),cw.View(),cw0,label+" Rev C,C");
    DoTestMV3(ca.View(),cwx.col(0).Reverse(),cwx.col(1),cw0,label+" Rev Step C,C");
    DoTestMV3(ca.View(),cwx.col(0),cwx.col(1).Reverse(),cw0,label+" Step Rev C,C");
    DoTestMV3(ca.View(),cwx.col(0).Reverse(),cwx.col(1).Reverse(),cw0,label+" Rev Rev C,C");
  }
#endif

#ifdef XTEST
  DoTestMM1(a.View(),b,m1,m2,label+" R,R");
  DoTestMM2(a.View(),b,m1,m2,a0,label+" R,R");
  m1 = a = a0;
  DoTestMM1(a.View(),cb,m1,cm2,label+" R,C");
  DoTestMM1(ca.View(),b,cm1,m2,label+" C,R");
  DoTestMM2(ca.View(),b,cm1,m2,ca0,label+" C,R");
  cm1 = ca = ca0;
#endif
  DoTestMM1(ca.View(),cb,cm1,cm2,label+" C,C");
  DoTestMM2(ca.View(),cb,cm1,cm2,ca0,label+" C,C");
  cm1 = ca = ca0;
  

  if (showstartdone)
    std::cout<<"Done Test"<<std::endl;
}

