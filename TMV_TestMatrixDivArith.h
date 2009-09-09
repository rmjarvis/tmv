#include "TMV.h"

template <class SM1, class SM2> inline bool CanLDiv(const SM1& a, const SM2& b)
{ return a.colsize() == b.colsize(); }

template <class T, class SM2> inline bool CanLDiv(const tmv::Vector<T>& a, const SM2& b)
{ return a.size() == b.colsize(); }

template <class SM1, class SM2> inline bool CanLDivEq(const SM1& a, const SM2& b)
{ return CanLDiv(a,b) && b.IsSquare(); }

template <class SM1, class SM2> inline bool CanRDiv(const SM1& a, const SM2& b)
{ return a.rowsize() == b.rowsize(); }

template <class T, class SM2> inline bool CanRDiv(const tmv::Vector<T>& a, const SM2& b)
{ return a.size() == b.rowsize(); }

template <class SM1, class SM2> inline bool CanRDivEq(const SM1& a, const SM2& b)
{ return CanRDiv(a,b) && b.IsSquare(); }

template <class SM1, class SM2, class T, tmv::StorageType S, class T2, tmv::StorageType S2> inline void DoTestMatrixDivArithMM1b(
    const SM1& a, const SM2& b, const tmv::Matrix<T,S>& m1,
    const tmv::Matrix<T2,S2>& m2, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM1b: "<<label<<std::endl;
    std::cout<<"b.dt = "<<b.GetDivType()<<std::endl;
    std::cout<<"m2.dt = "<<m2.GetDivType()<<std::endl;
  }
  b.SaveDiv();
  m2.SaveDiv();
  b.SetDiv();
  m2.SetDiv();

//#define XXDEBUG

#ifdef XXDEBUG
  std::cout <<"b.dt = "<<Text(b.GetDivType())<<" = "<<b.GetDiv()->Type()<<std::endl;
  std::cout <<"m.dt = "<<Text(m2.GetDivType())<<" = "<<m2.GetDiv()->Type()<<std::endl;
  std::cout<<"a = "<<tmv::Type(a)<<"  "<<a<<std::endl;
  std::cout<<"b = "<<tmv::Type(b)<<"  "<<b<<std::endl;
  std::cout<<"m1 = "<<tmv::Type(m1)<<"  "<<m1<<std::endl;
  std::cout<<"m2 = "<<tmv::Type(m2)<<"  "<<m2<<std::endl;
#endif
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a!=m1");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b!=m2");

  tmv::Matrix<T2,tmv::ColMajor> m2inv = Inverse(m2);
  double eps = EPS*Norm(m1)*std::max(m2.colsize(),m2.rowsize())*Norm(m2)*Norm(m2inv);

  if (CanLDiv(a,b)) {
#ifdef XXDEBUG
    std::cout<<"LDiv:\n";
    std::cout<<"a = "<<tmv::Type(a)<<a<<std::endl;
    std::cout<<"m1 = "<<tmv::Type(m1)<<m1<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<b<<std::endl;
    std::cout<<"m2 = "<<tmv::Type(m2)<<m2<<std::endl;
    std::cout<<"m1/m2 = "<<m1/m2<<std::endl;
    std::cout<<"a/m2 = "<<a/m2<<std::endl;
    std::cout<<"m1/b = "<<m1/b<<std::endl;
    std::cout<<"a/b = "<<a/b<<std::endl;
    std::cout<<"Norm(diff) = "<<Norm(a/b-m1/m2)<<" <?  "<<eps<<std::endl;
    std::cout<<"Norm(diff) = "<<Norm(a/m2-m1/m2)<<" <?  "<<eps<<std::endl;
    std::cout<<"Norm(diff) = "<<Norm(m1/b-m1/m2)<<" <?  "<<eps<<std::endl;
#endif
    Assert(Norm((a/b)-(m1/m2)) <= eps,label+" a/b");
    Assert(Norm((a/m2)-(m1/m2)) <= eps,label+" a/m");
    Assert(Norm((m1/b)-(m1/m2)) <= eps,label+" m/b");
    Assert(Norm((b.Inverse()*a)-(m1/m2)) <= eps,label+" b^-1*a");
    Assert(Norm((m2.Inverse()*a)-(m1/m2)) <= eps,label+" m^-1*a");
    Assert(Norm((b.Inverse()*m1)-(m1/m2)) <= eps,label+" b^-1*m");
  }
  if (CanRDiv(a,b)) {
#ifdef XXDEBUG
    std::cout<<"RDiv:\n";
    std::cout<<"a = "<<tmv::Type(a)<<a<<std::endl;
    std::cout<<"m1 = "<<tmv::Type(m1)<<m1<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<b<<std::endl;
    std::cout<<"m2 = "<<tmv::Type(m2)<<m2<<std::endl;
    std::cout<<"m1%m2 = "<<m1%m2<<std::endl;
    std::cout<<"a%m2 = "<<a%m2<<std::endl;
    std::cout<<"m1%b = "<<m1%b<<std::endl;
    std::cout<<"a%b = "<<a%b<<std::endl;
    std::cout<<"Norm(diff) = "<<Norm(a%b-m1%m2)<<" <?  "<<eps<<std::endl;
    std::cout<<"Norm(diff) = "<<Norm(a%m2-m1%m2)<<" <?  "<<eps<<std::endl;
    std::cout<<"Norm(diff) = "<<Norm(m1%b-m1%m2)<<" <?  "<<eps<<std::endl;
#endif
    Assert(Norm((a%b)-(m1%m2)) <= eps,label+" a%b");
    Assert(Norm((a%m2)-(m1%m2)) <= eps,label+" a%m");
    Assert(Norm((m1%b)-(m1%m2)) <= eps,label+" m%b");
    Assert(Norm((a*b.Inverse())-(m1%m2)) <= eps,label+" a*b^-1");
    Assert(Norm((a*m2.Inverse())-(m1%m2)) <= eps,label+" a*m^-1");
    Assert(Norm((m1*b.Inverse())-(m1%m2)) <= eps,label+" m*b^-1");
  }

#ifdef XXDEBUG
#undef XXDEBUG
#endif
  if (showstartdone) 
    std::cout<<"Done MM1b"<<std::endl;
}

template <class SM1, class SM2, class T, tmv::StorageType S, class T2, tmv::StorageType S2> inline void DoTestMatrixDivArithMM1a(
    tmv::DivType dt, const SM1& a, const SM2& b, const tmv::Matrix<T,S>& m1,
    const tmv::Matrix<T2,S2>& m2, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM1a: "<<label<<std::endl;
    std::cout<<"b.dt = "<<b.GetDivType()<<std::endl;
    std::cout<<"m2.dt = "<<m2.GetDivType()<<std::endl;
  }
  b.DivideUsing(dt);
  DoTestMatrixDivArithMM1b(a,b,m1,m2,label);
  if (showstartdone) 
    std::cout<<"Done MM1a"<<std::endl;
}

template <class SM1, class SM2, class T, class T2> inline void DoTestMatrixDivArithMM1(
    tmv::DivType dt, const SM1& a, const SM2& b, const tmv::Matrix<T,tmv::RowMajor>& m1,
    const tmv::Matrix<T2,tmv::RowMajor>& m2, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM1: "<<label<<std::endl;
    std::cout<<"b.dt = "<<b.GetDivType()<<std::endl;
    std::cout<<"m2.dt = "<<m2.GetDivType()<<std::endl;
  }
  b.DivideUsing(dt);
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a!=m1");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b!=m2");

  tmv::Matrix<T,tmv::ColMajor> m1t = m1.Transpose();
  tmv::Matrix<T,tmv::RowMajor> m1c = m1.Conjugate();
  tmv::Matrix<T,tmv::ColMajor> m1a = m1.Adjoint();
  tmv::Matrix<T2,tmv::ColMajor> m2t = m2.Transpose();
  tmv::Matrix<T2,tmv::RowMajor> m2c = m2.Conjugate();
  tmv::Matrix<T2,tmv::ColMajor> m2a = m2.Adjoint();

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
    std::cout<<"Done MM1"<<std::endl;
}

template <class SM1, class SM2, class T, tmv::StorageType S, class T2, tmv::StorageType S2> inline void DoTestMatrixDivArithMM2b(
    const SM1& a, const SM2& b, const tmv::Matrix<T,S>& m1,
    const tmv::Matrix<T2,S2>& m2, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM2b: "<<label<<std::endl;
    std::cout<<"b.dt = "<<b.GetDivType()<<std::endl;
    std::cout<<"m2.dt = "<<m2.GetDivType()<<std::endl;
  }

//#define XXDEBUG

#ifdef XXDEBUG
  std::cout<<"a = "<<tmv::Type(a)<<"  "<<a<<std::endl;
  std::cout<<"b = "<<tmv::Type(b)<<"  "<<b<<std::endl;
  std::cout<<"m1 = "<<tmv::Type(m1)<<"  "<<m1<<std::endl;
  std::cout<<"m2 = "<<tmv::Type(m2)<<"  "<<m2<<std::endl;
#endif

  m2.SaveDiv();
  b.SaveDiv();

  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a!=m1");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b!=m2");

  tmv::Matrix<T2,tmv::ColMajor> m2inv = Inverse(m2);
  double eps = EPS*Norm(m1)*std::max(m2.colsize(),m2.rowsize())*Norm(m2)*Norm(m2inv);

  if (CanLDivEq(m1,b)) {
    tmv::Matrix<T,S> m3 = m1;
    tmv::Matrix<T,S> m4 = m1;
    m3 /= b;
    m4 = tmv::Matrix<T,S>(m4/m2);
    Assert(Norm(m3-m4) <= eps,label+" m/=b");
  }
  if (CanRDivEq(m1,b)) {
    tmv::Matrix<T,S> m3 = m1;
    tmv::Matrix<T,S> m4 = m1;
    m3 %= b;
    m4 = tmv::Matrix<T,S>(m4%m2);
#ifdef XXDEBUG
    std::cout<<"m3 = "<<m3<<std::endl;
    std::cout<<"m4 = "<<m4<<std::endl;
    std::cout<<"m3*b = "<<m3*b<<std::endl;
    std::cout<<"m4*m2 = "<<m4*m2<<std::endl;
    std::cout<<"Norm(m3-m4) = "<<Norm(m3-m4)<<std::endl;
    std::cout<<"EPS*Norm(m3)*Norm(m2) = "<<EPS*Norm(m3)*Norm(m2)<<std::endl;
#endif
    Assert(Norm(m3-m4) <= eps,label+" m%=b");
    m3 = m1;
    m4 = m1;
    m3 *= b.Inverse();
    m4 = tmv::Matrix<T,S>(m4%m2);
    Assert(Norm(m3-m4) <= eps,label+" m*=b^-1");
  }

#ifdef XXDEBUG
#undef XXDEBUG
#endif
  if (showstartdone) 
    std::cout<<"Done MM2b"<<std::endl;
}

template <class SM1, class SM2, class T, tmv::StorageType S, class T2, tmv::StorageType S2> inline void DoTestMatrixDivArithMM2a(
    tmv::DivType dt, const SM1& a, const SM2& b, const tmv::Matrix<T,S>& m1,
    const tmv::Matrix<T2,S2>& m2, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM2a: "<<label<<std::endl;
    std::cout<<"b.dt = "<<b.GetDivType()<<std::endl;
    std::cout<<"m2.dt = "<<m2.GetDivType()<<std::endl;
  }
  b.DivideUsing(dt);
  DoTestMatrixDivArithMM2b(a,b,m1,m2,label);
  if (showstartdone) 
    std::cout<<"Done MM2a"<<std::endl;
}

template <class SM1, class SM2, class T, class T2> inline void DoTestMatrixDivArithMM2(
    tmv::DivType dt, const SM1& a, const SM2& b, const tmv::Matrix<T,tmv::RowMajor>& m1,
    const tmv::Matrix<T2,tmv::RowMajor>& m2, std::string label)
{
  b.DivideUsing(dt);
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a!=m1");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b!=m2");

  tmv::Matrix<T,tmv::ColMajor> m1t = m1.Transpose();
  tmv::Matrix<T,tmv::RowMajor> m1c = m1.Conjugate();
  tmv::Matrix<T,tmv::ColMajor> m1a = m1.Adjoint();
  tmv::Matrix<T2,tmv::ColMajor> m2t = m2.Transpose();
  tmv::Matrix<T2,tmv::RowMajor> m2c = m2.Conjugate();
  tmv::Matrix<T2,tmv::ColMajor> m2a = m2.Adjoint();

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

template <class SM, class T, class T2, tmv::StorageType S2> inline void DoTestMatrixDivArithMV1b(
    const tmv::Vector<T>& v, const SM& b, const tmv::Matrix<T2,S2>& m, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV1b: "<<label<<std::endl;
    std::cout<<"b.dt = "<<b.GetDivType()<<std::endl;
    std::cout<<"m.dt = "<<m.GetDivType()<<std::endl;
  }

  m.SaveDiv();
  b.SaveDiv();

  Assert(Norm(b-m) <= EPS*Norm(m),label+" b!=m");

  tmv::Matrix<T2,tmv::ColMajor> minv = Inverse(m);
  double eps = EPS*Norm(v)*std::max(m.colsize(),m.rowsize())*Norm(m)*Norm(minv);

  if (CanLDiv(v,b)) {
//#define XXDEBUG
#ifdef XXDEBUG
    std::cerr<<"v = "<<tmv::Type(v)<<"  step "<<v.step()<<"  "<<v<<std::endl;
    std::cerr<<"b = "<<tmv::Type(b)<<"  "<<b<<std::endl;
    std::cerr<<"m = "<<tmv::Type(m)<<"  "<<m<<std::endl;
    std::cerr<<"v/b = "<<v/b<<std::endl;
    std::cerr<<"v/m = "<<v/m<<std::endl;
    std::cerr<<"Norm(diff) = "<<Norm((v/b)-(v/m))<<std::endl;
    std::cerr<<"eps = "<<EPS<<"*"<<Norm(v)<<"*"<<std::max(m.colsize(),m.rowsize())<<"*"<<Norm(m)<<"*"<<Norm(minv)<<"="<<eps<<std::endl;
#undef XXDEBUG
#endif
    Assert(Norm((v/b)-(v/m)) <= eps,label+" v/m");
    Assert(Norm((b.Inverse()*v)-(v/m)) <= eps,label+" m^-1*v");
  }
  if (CanRDiv(v,b)) {
//#define XXDEBUG
#ifdef XXDEBUG
    std::cerr<<"v = "<<tmv::Type(v)<<"  step "<<v.step()<<"  "<<v<<std::endl;
    std::cerr<<"b = "<<tmv::Type(b)<<"  "<<b<<std::endl;
    std::cerr<<"m = "<<tmv::Type(m)<<"  "<<m<<std::endl;
    std::cerr<<"v%b = "<<v%b<<std::endl;
    std::cerr<<"v%m = "<<v%m<<std::endl;
    std::cerr<<"(v%b)*m = "<<(v%b)*m<<std::endl;
    std::cerr<<"(v%b)*b = "<<(v%b)*b<<std::endl;
    std::cerr<<"(v%m)*m = "<<(v%m)*m<<std::endl;
    std::cerr<<"(v%m)*b = "<<(v%m)*b<<std::endl;
    std::cerr<<"diff = "<<(v%b)-(v%m)<<std::endl;
    std::cerr<<"Norm(diff) = "<<Norm((v%b)-(v%m))<<std::endl;
    std::cerr<<"eps = "<<EPS<<"*"<<Norm(v)<<"*"<<std::max(m.colsize(),m.rowsize())<<"*"<<Norm(m)<<"*"<<Norm(minv)<<"="<<eps<<std::endl;
#undef XXDEBUG
#endif
    Assert(Norm((v%b)-(v%m)) <= eps,label+" v%m");
    Assert(Norm((v*b.Inverse())-(v%m)) <= eps,label+" v*m^-1");
  }

  if (showstartdone) 
    std::cout<<"Done MV1b"<<std::endl;
}

template <class SM, class T, class T2, tmv::StorageType S2> inline void DoTestMatrixDivArithMV1a(
    tmv::DivType dt, const tmv::Vector<T>& v, const SM& b, 
    const tmv::Matrix<T2,S2>& m, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV1a: "<<label<<std::endl;
    std::cout<<"b.dt = "<<b.GetDivType()<<std::endl;
    std::cout<<"m.dt = "<<m.GetDivType()<<std::endl;
  }
  b.DivideUsing(dt);
  DoTestMatrixDivArithMV1b(v,b,m,label);
  if (showstartdone) 
    std::cout<<"Done MV1a"<<std::endl;
}

template <class SM, class T, class T2> inline void DoTestMatrixDivArithMV1(
    tmv::DivType dt, const tmv::Vector<T>& v, const SM& b,
    const tmv::Matrix<T2,tmv::RowMajor>& m, std::string label)
{
  b.DivideUsing(dt);
  Assert(Norm(b-m) <= EPS*Norm(m),label+" b!=m");

  tmv::Matrix<T2,tmv::ColMajor> mt = m.Transpose();
  tmv::Matrix<T2,tmv::RowMajor> mc = m.Conjugate();
  tmv::Matrix<T2,tmv::ColMajor> ma = m.Adjoint();

  DoTestMatrixDivArithMV1a(dt,v,b.View(),m,label);
  DoTestMatrixDivArithMV1a(dt,v,Transpose(b),mt,label+" Trans");
  DoTestMatrixDivArithMV1a(dt,v,Conjugate(b),mc,label+" Conj");
  DoTestMatrixDivArithMV1a(dt,v,Adjoint(b),ma,label+" Adj");
}

template <class SM, class T, class T2, tmv::StorageType S2> inline void DoTestMatrixDivArithMV2b(
    const tmv::Vector<T>& v, const SM& b, const tmv::Matrix<T2,S2>& m, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV2b: "<<label<<std::endl;
    std::cout<<"b.dt = "<<b.GetDivType()<<std::endl;
    std::cout<<"m.dt = "<<m.GetDivType()<<std::endl;
  }
  Assert(Norm(b-m) <= EPS*Norm(m),label+" b!=m");

  m.SaveDiv();
  b.SaveDiv();

  tmv::Matrix<T2,tmv::ColMajor> minv = Inverse(m);
  double eps = EPS*Norm(v)*std::max(m.colsize(),m.rowsize())*Norm(m)*Norm(minv);

  if (CanLDivEq(v,m)) {
    tmv::Vector<T> v1 = v;
    tmv::Vector<T> v2 = v;
    v1 /= b;
    v2 = tmv::Vector<T>(v2/m);
    Assert(Norm(v1-v2) <= eps,label+" v/=m");
  }
  if (CanRDivEq(v,m)) {
    tmv::Vector<T> v1 = v;
    tmv::Vector<T> v2 = v;
    v1 %= b;
    v2 = tmv::Vector<T>(v2%m);
    Assert(Norm(v1-v2) <= eps,label+" v%=m");
    v1 = v;
    v2 = v;
    v1 *= b.Inverse();
    v2 = tmv::Vector<T>(v2%m);
    Assert(Norm(v1-v2) <= eps,label+" v*=m^-1");
  }
  if (showstartdone) 
    std::cout<<"Done MV2b"<<std::endl;
}

template <class SM, class T, class T2, tmv::StorageType S2> inline void DoTestMatrixDivArithMV2a(
    tmv::DivType dt, const tmv::Vector<T>& v, const SM& b,
    const tmv::Matrix<T2,S2>& m, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV2a: "<<label<<std::endl;
    std::cout<<"b.dt = "<<b.GetDivType()<<std::endl;
    std::cout<<"m.dt = "<<m.GetDivType()<<std::endl;
  }
  b.DivideUsing(dt);
  DoTestMatrixDivArithMV2b(v,b,m,label);
  if (showstartdone) 
    std::cout<<"Done MV2a"<<std::endl;
}

template <class SM, class T, class T2> inline void DoTestMatrixDivArithMV2(
    tmv::DivType dt, const tmv::Vector<T>& v, const SM& b,
    const tmv::Matrix<T2,tmv::RowMajor>& m, std::string label)
{
  b.DivideUsing(dt);
  Assert(Norm(b-m) <= EPS*Norm(m),label+" b!=m");

  tmv::Matrix<T2,tmv::ColMajor> mt = m.Transpose();
  tmv::Matrix<T2,tmv::RowMajor> mc = m.Conjugate();
  tmv::Matrix<T2,tmv::ColMajor> ma = m.Adjoint();

  DoTestMatrixDivArithMV2a(dt,v,b.View(),m,label);
  DoTestMatrixDivArithMV2a(dt,v,Transpose(b),mt,label+" TransA");
  DoTestMatrixDivArithMV2a(dt,v,Conjugate(b),mc,label+" ConjA");
  DoTestMatrixDivArithMV2a(dt,v,Adjoint(b),ma,label+" AdjA");
}

template <class SM1, class T, tmv::StorageType S, class T2> inline void DoTestMatrixDivArithMXb(
    const SM1& a, const tmv::Matrix<T,S>& m, T2 x, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MXb: "<<label<<std::endl;
    std::cout<<"a.dt = "<<a.GetDivType()<<std::endl;
    std::cout<<"m.dt = "<<m.GetDivType()<<std::endl;
  }

  m.SaveDiv();
  a.SaveDiv();

  tmv::Matrix<T,tmv::ColMajor> minv = m.Inverse();
  double normm = Norm(m);
  double eps = EPS*std::max(m.colsize(),m.rowsize())*Norm(m)*Norm(minv);

//#define XXDEBUG

#ifdef XXDEBUG
  std::cout<<"eps = "<<eps<<std::endl;
  std::cout<<"x = "<<x<<std::endl;
  std::cout<<"a = "<<tmv::Type(a)<<"  "<<a<<std::endl;
  std::cout<<"m = "<<tmv::Type(m)<<"  "<<m<<std::endl;
  std::cout<<"x/a = "<<x/a<<std::endl;
  std::cout<<"x/m = "<<x/m<<std::endl;
  std::cout<<"a*(x/a) = "<<a*(x/a)<<std::endl;
  std::cout<<"m*(x/m) = "<<m*(x/m)<<std::endl;
  std::cout<<"Norm(diff) = "<<Norm((x/a)-(x/m))<<std::endl;
  std::cout<<"eps*x = "<<eps*std::abs(x)<<std::endl;
#endif
  Assert(Norm((x/a)-(x/m)) <= eps*std::abs(x),label+" x/a");
  Assert(Norm((a.Inverse()*x)-(x/m)) <= eps*std::abs(x),label+" a^-1*x");
#ifdef XXDEBUG
  std::cout<<"x%a = "<<x%a<<std::endl;
  std::cout<<"x%m = "<<x%m<<std::endl;
  std::cout<<"x%a*a = "<<(x%a)*a<<std::endl;
  std::cout<<"x%m*m = "<<(x%m)*m<<std::endl;
  std::cout<<"Norm(diff) = "<<Norm((x%a)-(x%m))<<std::endl;
#endif
  Assert(Norm((x%a)-(x%m)) <= eps*std::abs(x),label+" x%a");
  Assert(Norm((x*a.Inverse())-(x%m)) <= eps*std::abs(x),label+" x*a^-1");
#ifdef XXDEBUG
  std::cout<<"a = "<<tmv::Type(a)<<"  "<<a<<std::endl;
  std::cout<<"1/a = "<<T2(1)/a<<std::endl;
  std::cout<<"a*(1/a) = "<<a*(T2(1)/a)<<std::endl;
  std::cout<<"(a*(1/a))*a = "<<(a*(T2(1)/a))*a<<std::endl;
  std::cout<<"(a*(1/a))*a-a = "<<(a*(T2(1)/a))*a-a<<std::endl;
  std::cout<<Norm(a*((T2(1)/a)*a)-a)<<"  "<<eps*normm<<std::endl;
  std::cout<<Norm((a*(T2(1)/a))*a-a)<<"  "<<eps*normm<<std::endl;
  std::cout<<Norm(a*((T2(1)%a)*a)-a)<<"  "<<eps*normm<<std::endl;
  std::cout<<Norm((a*(T2(1)%a))*a-a)<<"  "<<eps*normm<<std::endl;
  std::cout<<Norm((T2(1)/a)*a-Transpose((T2(1)/a)*a))<<"  "<<eps<<std::endl;
  std::cout<<Norm((T2(1)%a)*a-Transpose((T2(1)%a)*a))<<"  "<<eps<<std::endl;
  std::cout<<Norm(a*(T2(1)/a)-Transpose(a*(T2(1)/a)))<<"  "<<eps<<std::endl;
  std::cout<<Norm(a*(T2(1)%a)-Transpose(a*(T2(1)%a)))<<"  "<<eps<<std::endl;
#endif
  Assert(Norm(a*((T2(1)/a)*a)-a) <= eps*normm,label+" a*(1/a)*a");
  Assert(Norm((a*(T2(1)/a))*a-a) <= eps*normm,label+" a*(1/a)*a");
  Assert(Norm(a*((T2(1)%a)*a)-a) <= eps*normm,label+" a*(1%a)*a");
  Assert(Norm((a*(T2(1)%a))*a-a) <= eps*normm,label+" a*(1%a)*a");
  Assert(Norm((T2(1)/a)*a-Transpose((T2(1)/a)*a)) <= 2*eps,label+" (1/a)*a-((1/a)*a)T");
  Assert(Norm((T2(1)%a)*a-Transpose((T2(1)%a)*a)) <= 2*eps,label+" (1/a)*a-((1/a)*a)T");
  Assert(Norm(a*(T2(1)/a)-Transpose(a*(T2(1)/a))) <= 2*eps,label+" a*(1/a)-(a*(1/a))T");
  Assert(Norm(a*(T2(1)%a)-Transpose(a*(T2(1)%a))) <= 2*eps,label+" a*(1/a)-(a*(1/a))T");

#ifdef XXDEBUG
#undef XXDEBUG
#endif
  if (showstartdone) 
    std::cout<<"Done MXb"<<std::endl;
}

template <class SM1, class T, tmv::StorageType S, class T2> inline void DoTestMatrixDivArithMXa(
    tmv::DivType dt, const SM1& a, const tmv::Matrix<T,S>& m, T2 x, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MXa: "<<label<<std::endl;
    std::cout<<"a.dt = "<<a.GetDivType()<<std::endl;
    std::cout<<"m.dt = "<<m.GetDivType()<<std::endl;
  }
  a.DivideUsing(dt);
  DoTestMatrixDivArithMXb(a,m,x,label);
  if (showstartdone) 
    std::cout<<"Done MXa"<<std::endl;
}

template <class SM1, class T, class T2> inline void DoTestMatrixDivArithMX(
    tmv::DivType dt, const SM1& a, const tmv::Matrix<T,tmv::RowMajor>& m, T2 x, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MX: "<<label<<std::endl;
    std::cout<<"a.dt = "<<a.GetDivType()<<std::endl;
    std::cout<<"m.dt = "<<m.GetDivType()<<std::endl;
  }
  tmv::Matrix<T,tmv::ColMajor> mt = m.Transpose();
  tmv::Matrix<T,tmv::RowMajor> mc = m.Conjugate();
  tmv::Matrix<T,tmv::ColMajor> ma = m.Adjoint();
  DoTestMatrixDivArithMXa(dt,a.View(),m,x,label);
  DoTestMatrixDivArithMXa(dt,Transpose(a),mt,x,label+" Trans");
  DoTestMatrixDivArithMXa(dt,Conjugate(a),mc,x,label+" Conj");
  DoTestMatrixDivArithMXa(dt,Adjoint(a),ma,x,label+" Adj");
  if (showstartdone) 
    std::cout<<"Done MX"<<std::endl;
}

template <class T, class SM1, class SM2, class CSM1, class CSM2> inline void TestMatrixDivArith(
    tmv::DivType dt, const SM1& a, const SM2& b, const CSM1& ca, const CSM2& cb, 
    std::string label)
{
  if (showstartdone) {
    std::cout<<"Start Test Div: "<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<"  "<<b<<std::endl;
    std::cout<<"ca = "<<tmv::Type(ca)<<"  "<<ca<<std::endl;
    std::cout<<"cb = "<<tmv::Type(cb)<<"  "<<cb<<std::endl;
    std::cout<<"a.dt = "<<a.GetDivType()<<std::endl;
    std::cout<<"b.dt = "<<b.GetDivType()<<std::endl;
    std::cout<<"ca.dt = "<<ca.GetDivType()<<std::endl;
    std::cout<<"cb.dt = "<<cb.GetDivType()<<std::endl;
  }
  
  a.SaveDiv();
  b.SaveDiv();
  ca.SaveDiv();
  cb.SaveDiv();

  tmv::Matrix<T,tmv::RowMajor> m1(a);
  tmv::Matrix<T,tmv::RowMajor> m2(b);
  tmv::Matrix<std::complex<T>,tmv::RowMajor> cm1(ca);
  tmv::Matrix<std::complex<T>,tmv::RowMajor> cm2(cb);

  m1.SaveDiv();
  m2.SaveDiv();
  cm1.SaveDiv();
  cm2.SaveDiv();

  std::complex<T> z1(9,-2);
#ifdef XTEST
  T x = 12;

  DoTestMatrixDivArithMX(dt,a.View(),m1,x,label+" R,R");
  DoTestMatrixDivArithMX(dt,a.View(),m1,z1,label+" R,C");
  DoTestMatrixDivArithMX(dt,ca.View(),cm1,x,label+" C,R");
#endif
  DoTestMatrixDivArithMX(dt,ca.View(),cm1,z1,label+" C,C");

#ifdef XTEST
  DoTestMatrixDivArithMM1(dt,b.View(),a.View(),m2,m1,label+" R,R");
  DoTestMatrixDivArithMM2(dt,b.View(),a.View(),m2,m1,label+" R,R");
  DoTestMatrixDivArithMM1(dt,b.View(),ca.View(),m2,cm1,label+" R,C");
  DoTestMatrixDivArithMM1(dt,cb.View(),a.View(),cm2,m1,label+" C,R");
  DoTestMatrixDivArithMM2(dt,cb.View(),a.View(),cm2,m1,label+" C,R");
#endif
  DoTestMatrixDivArithMM1(dt,cb.View(),ca.View(),cm2,cm1,label+" C,C");
  DoTestMatrixDivArithMM2(dt,cb.View(),ca.View(),cm2,cm1,label+" C,C");

#ifdef XTEST
  tmv::Vector<T> v(b.colsize());
  tmv::Vector<std::complex<T> > cv(b.colsize());
  if (b.rowsize() > 0) {
    v = m2.col(0);
    cv = cm2.col(0);
  }
  else {
    v.SetAllTo(1);
    cv.SetAllTo(std::complex<T>(1,2));
  }

  tmv::Vector<T> w(b.rowsize());
  tmv::Vector<std::complex<T> > cw(b.rowsize());
  if (b.colsize() > 0) {
    w = m2.row(0);
    cw = cm2.row(0);
  }
  else {
    w.SetAllTo(1);
    cw.SetAllTo(std::complex<T>(1,2));
  }

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

