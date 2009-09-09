#define CT std::complex<T>

template <class M1, class M2> inline bool CanAdd(const M1& a, const M2& b)
{ return a.colsize() == b.colsize() && a.rowsize() == b.rowsize(); }

template <class M1, class M2> inline bool CanAddEq(const M1& a, const M2& b)
{ return CanAdd(a,b); }

template <class M, class T2> inline bool CanAddX(const M& a, const T2)
{ return a.IsSquare(); }

template <class M, class T2> inline bool CanAddEqX(const M& a, const T2 x)
{ return CanAddX(a,x); }

template <class M, class T2> inline bool CanMultX(const M&, const T2)
{ return true; }

template <class M, class T2> inline bool CanMultEqX(const M& a, const T2 x)
{ return CanMultX(a,x); }

template <class M1, class M2> inline bool CanMult(const M1& a, const M2& b)
{ return a.rowsize() == b.colsize(); }

template <class M, class T> inline bool CanMult(const M& m, const tmv::Vector<T>& v)
{ return m.rowsize() == v.size(); }

template <class M, class T> inline bool CanMult(const tmv::Vector<T>& v, const M& m)
{ return v.size() == m.colsize(); }

template <class M1, class M2, class M3> inline bool CanMult(
    const M1& a, const M2& b, const M3& c)
{
  return CanMult(a,b) && c.colsize() == a.colsize() && 
    c.rowsize() == b.rowsize(); 
}

template <class M1, class M2, class M3> inline bool CanMultXM(const M1& a, const M2& b, const M3& c)
{ return CanMult(a,b,c); }

template <class T, class M> inline void DoTestMa(const M& a, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start Ma "<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<std::endl;
  }

  tmv::Matrix<T> m = a;
  double eps = EPS * Norm(m);

  if (XXDEBUG1) {
    std::cout<<"a = "<<tmv::Type(a)<<" = "<<a<<std::endl;
    std::cout<<"m = "<<tmv::Type(m)<<" = "<<m<<std::endl;
    std::cout<<"a-m = "<<a-m<<std::endl;
    std::cout<<"Norm(a-m) = "<<Norm(a-m)<<std::endl;
    std::cout<<"Trace(a) = "<<Trace(a)<<"  "<<Trace(m)<<std::endl;
    std::cout<<"NormF(a) = "<<NormF(a)<<"  "<<NormF(m)<<std::endl;
    std::cout<<"Norm(a) = "<<Norm(a)<<"  "<<Norm(m)<<std::endl;
    std::cout<<"Norm1(a) = "<<Norm1(a)<<"  "<<Norm1(m)<<std::endl;
    std::cout<<"Norm1(m) = "<<Norm1(m)<<std::endl;
    std::cout<<"NormInf(a) = "<<NormInf(a)<<"  "<<NormInf(m)<<std::endl;
    std::cout<<"abs(diff) = "<<std::abs(NormInf(a)-NormInf(m))<<std::endl;
    std::cout<<"eps*norminf = "<<EPS*NormInf(m)<<std::endl;
  }

  Assert(Norm(a-m) <= eps,label+" a != m");
  Assert(std::abs(Trace(a)-Trace(m)) <= eps,label+" Trace");
  Assert(std::abs(NormF(a)-NormF(m)) <= eps,label+" NormF");
  Assert(std::abs(Norm(a)-Norm(m)) <= eps,label+" Norm");
  Assert(std::abs(Norm1(a)-Norm1(m)) <= eps,label+" Norm1");
  Assert(std::abs(NormInf(a)-NormInf(m)) <= eps,label+" NormInf");
#ifndef NOSV
  if (donorm2) {
    if (XXDEBUG1) {
      std::cout<<"Norm2(a) = "<<a.DoNorm2()<<"  "<<m.DoNorm2()<<std::endl;
      std::cout<<"abs(diff) = "<<std::abs(a.DoNorm2()-m.DoNorm2())<<std::endl;
      std::cout<<"eps*kappa = "<<eps*m.DoCondition()<<std::endl;
    }
    Assert(std::abs(a.DoNorm2()-m.DoNorm2()) <= eps*m.DoCondition(),label+" DoNorm2");
    a.DivideUsing(tmv::SV);
    m.DivideUsing(tmv::SV);
    a.SetDiv();
    m.SetDiv();
    Assert(std::abs(Norm2(a)-m.DoNorm2()) <= eps*m.DoCondition(),label+" Norm2");
  }
#endif
  if (XXDEBUG1) {
    std::cout<<"Norm(aT-mT) = "<<Norm(Transpose(a)-Transpose(m))<<std::endl;
    std::cout<<"Conjugate(a) = "<<Conjugate(a)<<std::endl;
    std::cout<<"Conjugate(m) = "<<Conjugate(m)<<std::endl;
    std::cout<<"a*-m* = "<<Conjugate(a)-Conjugate(m)<<std::endl;
    std::cout<<"Conjugate(a).diag = "<<Conjugate(a).diag()<<std::endl;
    std::cout<<"Conjugate(m).diag = "<<Conjugate(m).diag()<<std::endl;
    std::cout<<"Norm(a*-m*) = "<<Norm(Conjugate(a)-Conjugate(m))<<std::endl;
    std::cout<<"Norm(at-mt) = "<<Norm(Adjoint(a)-Adjoint(m))<<std::endl;
  }
  Assert(Norm(Transpose(a)-Transpose(m)) <= eps,label+" Transpose");
  Assert(Norm(Conjugate(a)-Conjugate(m)) <= eps,label+" Conjugate");
  Assert(Norm(Adjoint(a)-Adjoint(m)) <= eps,label+" Adjoint");

  if (showstartdone)
    std::cout<<"Done Ma"<<std::endl;
}

template <class T, class M> inline void DoTestMR(
    const M& a, std::string label)
{
  DoTestMa<T>(a.View(),label);
  DoTestMa<T>(Transpose(a),label+" Trans");
}

template <class T, class M> inline void DoTestMC(
    const M& a, std::string label)
{
  DoTestMa<T>(a.View(),label);
  DoTestMa<T>(Transpose(a),label+" Trans");
  DoTestMa<T>(Conjugate(a),label+" Conj");
  DoTestMa<T>(Adjoint(a),label+" Adj");
}

template <class T, class M, class T2> inline void DoTestMX1a(
    const M& a, T2 x, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MX1a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"x = "<<tmv::Type(x)<<" "<<x<<std::endl;
  }

  tmv::Matrix<T> m = a;

  double eps = EPS*Norm(m);
  if (CanAddX(a,x)) {
    if (XXDEBUG2) {
      std::cout<<"CanAddX("<<tmv::Type(a)<<","<<tmv::Type(x)<<")\n";
      std::cout<<"x = "<<tmv::Type(x)<<"  "<<x<<std::endl;
      std::cout<<"a = "<<tmv::Type(a)<<"  "<<a<<std::endl;
      std::cout<<"x-a = "<<(x-a)<<std::endl;
      std::cout<<"x-m = "<<(x-m)<<std::endl;
      std::cout<<"(x-a)-(x-m) = "<<(x-a)-(x-m)<<std::endl;
    }
    Assert(Norm(a-m) <= eps,label+" a != m1");
    Assert(Norm((x-a)-(x-m)) <= eps,label+" x-a");
#ifdef XTEST
    Assert(Norm((a-x)-(m-x)) <= eps,label+" a-x");
    Assert(Norm((x+a)-(x+m)) <= eps,label+" x+a");
    Assert(Norm((a+x)-(m+x)) <= eps,label+" a+x");
#endif
  }
  if (CanMultX(a,x)) {
    if (XXDEBUG2) {
      std::cout<<"CanMultX("<<tmv::Type(a)<<","<<tmv::Type(x)<<")\n";
    }
    Assert(Norm((x*a)-(x*m)) <= eps*std::abs(x),label+" x*a");
#ifdef XTEST
    Assert(Norm((a*x)-(x*m)) <= eps*std::abs(x),label+" a*x");
    if (tmv::Epsilon<T>()  != T(0)) {
      Assert(Norm((a/x)-(m/x)) <= eps/std::abs(x),label+" a/x");
    }
#endif
  }
  if (showstartdone)
    std::cout<<"Done MX1a"<<std::endl;
}

template <class T, class M, class T2> inline void DoTestMX1R(
    const M& a, T2 x, std::string label)
{
  DoTestMX1a<T>(a.View(),x,label);
}

template <class T, class M, class T2> inline void DoTestMX1C(
    const M& a, T2 x, std::string label)
{
  DoTestMX1a<T>(a.View(),x,label);
  DoTestMX1a<T>(Conjugate(a),x,label+" Conj");
}

template <class T, class BaseM, class M, class T2> inline void DoTestMX2a(
    BaseM& a0, const M& a, T2 x, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MX2a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"a0 = "<<tmv::Type(a0)<<" "<<a0<<std::endl;
    std::cout<<"x = "<<tmv::Type(x)<<" "<<x<<std::endl;
  }

  a0 = a;
  tmv::Matrix<T> m1 = a;
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a = m1");

  double eps = EPS * Norm(m1);
#ifndef NONSQUARE
  if (CanAddEqX(a,x)) {
    if (XXDEBUG3) {
      std::cout<<"CanAddEqX("<<tmv::Type(a)<<","<<tmv::Type(x)<<")\n";
    }
    tmv::Matrix<T> m2 = m1;
    a += x;
    m2 = m1+x;
    Assert(Norm(a-m2) <= eps,label+" a += x");
    Assert(Norm((a+=x)-(m2+=x)) <= eps,label+" a += x (2)");
    a = a0;
    a = a+x; 
    m2 = m1+x;
    Assert(Norm(a-m2) <= eps,label+" a = a+x");
    a = a0;
#ifdef XTEST
    a += -x;
    m2 = m1-x;
    Assert(Norm(a-m2) <= eps,label+" a += x");
    Assert(Norm((a+=-x)-(m2+=-x)) <= eps,label+" a += -x (2)");
    a = a0;
    a -= x;
    m2 = m1-x;
    Assert(Norm(a-m2) <= eps,label+" a -= x");
    Assert(Norm((a-=x)-(m2-=x)) <= eps,label+" a -= x (2)");
    a = a0;
    a -= -x;
    m2 = m1+x;
    Assert(Norm(a-m2) <= eps,label+" a -= x");
    Assert(Norm((a-=-x)-(m2-=-x)) <= eps,label+" a -= -x (2)");
    a = a0;
    a = a-x;
    m2 = m1-x;
    Assert(Norm(a-m2) <= eps,label+" a = a-x");
    a = a0;
    a = x+a;
    m2 = x+m1;
    Assert(Norm(a-m2) <= eps,label+" a = x+a");
    a = a0;
    a = x-a; 
    m2 = x-m1;
    Assert(Norm(a-m2) <= eps,label+" a = x-a");
    a = a0;
#endif
  }
#endif
  if (CanMultEqX(a,x)) {
    if (XXDEBUG3) {
      std::cout<<"CanMultEqX("<<tmv::Type(a)<<","<<tmv::Type(x)<<")\n";
    }
    tmv::Matrix<T> m2 = m1;
    a *= x;
    m2 = m1*x;
    Assert(Norm(a-m2) <= eps*std::abs(x),label+" a *= x");
    Assert(Norm((a*=x)-(m2*=x)) <= eps*std::abs(x*x),label+" a *= x");
    a = a0;
    a = a*x;
    m2 = m1*x;
    Assert(Norm(a-m2) <= eps*std::abs(x),label+" a = a*x");
    a = a0;
#ifdef XTEST
    a *= -x;
    m2 = -m1*x;
    Assert(Norm(a-m2) <= eps*std::abs(x),label+" a *= -x");
    Assert(Norm((a*=-x)-(m2*=-x)) <= eps*std::abs(x*x),label+" a *= -x");
    a = a0;
    if (tmv::Epsilon<T>() != T(0)) {
      a /= x;
      m2 = m1/x;
      Assert(Norm(a-m2) <= eps*std::abs(x),label+" a /= x");
      Assert(Norm((a/=x)-(m2/=x)) <= eps,label+" a /= x");
      a = a0;
      a /= -x;
      m2 = -m1/x;
      Assert(Norm(a-m2) <= eps*std::abs(x),label+" a /= -x");
      Assert(Norm((a/=-x)-(m2/=-x)) <= eps,label+" a /= -x");
      a = a0;
      a = a/x;
      m2 = m1/x;
      Assert(Norm(a-m2) <= eps,label+" a = a/x");
      a = a0;
    }
    a = x*a; 
    m2 = x*m1;
    Assert(Norm(a-m2) <= eps*std::abs(x),label+" a = x*a");
    a = a0;
#endif
  }

  if (showstartdone)
    std::cout<<"Done MX2a"<<std::endl;
}

template <class T, class BaseM, class M, class T2> inline void DoTestMX2R(
    BaseM& a0, const M& a, T2 x, std::string label)
{
  DoTestMX2a<T>(a0,a,x,label);
}

template <class T, class BaseM, class M, class T2> inline void DoTestMX2C(
    BaseM& a0, const M& a, T2 x, std::string label)
{
  DoTestMX2a<T>(a0,a,x,label);
  DoTestMX2a<T>(a0,Conjugate(a),x,label+" Conj");
}

// m*v, v*m
template <class Ta, class T, class M, class V> inline void DoTestMV1a(
    const M& a, const V& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV1a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" "<<b<<std::endl;
  }

  tmv::Matrix<Ta> m = a;
  tmv::Vector<T> v = b;
  double eps = EPS * Norm(m) * Norm(v);
  if (CanMult(m,ColVectorViewOf(b))) {
    if (XXDEBUG4) {
      std::cout<<"CanMult("<<tmv::Type(m)<<","<<tmv::Type(b)<<")\n";
    }
    Assert(Norm((a*b)-(m*v)) <= eps,label+" a*v");
  }
  if (showstartdone)
    std::cout<<"Done MV1a"<<std::endl;
}

template <class Ta, class T, class M, class V> inline void DoTestVM1a(
    const M& a, const V& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV1a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" "<<b<<std::endl;
  }

  tmv::Matrix<Ta> m = a;
  tmv::Vector<T> v = b;
  double eps = EPS * Norm(m) * Norm(v);
  if (CanMult(RowVectorViewOf(v),m)) {
    if (XXDEBUG4) {
      std::cout<<"CanMult("<<tmv::Type(v)<<","<<tmv::Type(m)<<")\n";
    }
    Assert(Norm((b*a)-(v*m)) <= eps,label+" v*a");
  }
  if (showstartdone)
    std::cout<<"Done VM1a"<<std::endl;
}

template <class Ta, class T, class M, class V> inline void DoTestMV1R(
    const M& a, const V& b, std::string label)
{
  DoTestMV1a<Ta,T>(a,b,label);
  DoTestMV1a<Ta,T>(a,b.Reverse(),label);
  DoTestVM1a<Ta,T>(a,b,label);
  DoTestVM1a<Ta,T>(a,b.Reverse(),label);

#ifdef XTEST
  tmv::Vector<T> b0 = b;

  b.Zero();
  DoTestMV1a<Ta,T>(a,b,label);
  DoTestVM1a<Ta,T>(a,b,label);

  b = b0;
  b.SubVector(0,b.size()/2).Zero();
  DoTestMV1a<Ta,T>(a,b,label);
  DoTestVM1a<Ta,T>(a,b,label);

  b = b0;
  b.SubVector(b.size()/2,b.size()).Zero();
  DoTestMV1a<Ta,T>(a,b,label);
  DoTestVM1a<Ta,T>(a,b,label);

  b = b0;
  b.SubVector(0,b.size()/4).Zero();
  b.SubVector(3*b.size()/4,b.size()).Zero();
  DoTestMV1a<Ta,T>(a,b,label);
  DoTestVM1a<Ta,T>(a,b,label);

  if (b.size() > 1) {
    b = b0;
    b.SubVector(0,1).Zero();
    DoTestMV1a<Ta,T>(a.View(),b,label);
    DoTestVM1a<Ta,T>(a.View(),b,label);

    b = b0;
    b.SubVector(b.size()-1,b.size()).Zero();
    DoTestMV1a<Ta,T>(a.View(),b,label);
    DoTestVM1a<Ta,T>(a.View(),b,label);

    b = b0;
    b.SubVector(0,1).Zero();
    b.SubVector(b.size()-1,b.size()).Zero();
    DoTestMV1a<Ta,T>(a.View(),b,label);
    DoTestVM1a<Ta,T>(a.View(),b,label);
  }
  b=b0;
#endif
}
template <class Ta, class T, class M, class V> inline void DoTestMV1C(
    const M& a, const V& b, std::string label)
{
  DoTestMV1a<Ta,T>(a,b,label);
  DoTestMV1a<Ta,T>(a,b.Reverse(),label);
  DoTestVM1a<Ta,T>(a,b,label);
  DoTestVM1a<Ta,T>(a,b.Reverse(),label);

  DoTestMV1a<Ta,T>(Conjugate(a),b,label+" Conj");
  DoTestVM1a<Ta,T>(Conjugate(a),b,label+" Conj");
#ifdef XTEST
  DoTestMV1a<Ta,T>(Conjugate(a),b.Reverse(),label+" Conj");
  DoTestVM1a<Ta,T>(Conjugate(a),b.Reverse(),label+" Conj");

  tmv::Vector<T> b0 = b;

  b.Zero();
  DoTestMV1a<Ta,T>(a,b,label);
  DoTestVM1a<Ta,T>(a,b,label);
  DoTestMV1a<Ta,T>(Conjugate(a),b,label+" Conj");
  DoTestVM1a<Ta,T>(Conjugate(a),b,label+" Conj");

  b = b0;
  b.SubVector(0,b.size()/2).Zero();
  DoTestMV1a<Ta,T>(a,b,label);
  DoTestVM1a<Ta,T>(a,b,label);
  DoTestMV1a<Ta,T>(Conjugate(a),b,label+" Conj");
  DoTestVM1a<Ta,T>(Conjugate(a),b,label+" Conj");

  b = b0;
  b.SubVector(b.size()/2,b.size()).Zero();
  DoTestMV1a<Ta,T>(a,b,label);
  DoTestVM1a<Ta,T>(a,b,label);
  DoTestMV1a<Ta,T>(Conjugate(a),b,label+" Conj");
  DoTestVM1a<Ta,T>(Conjugate(a),b,label+" Conj");

  b = b0;
  b.SubVector(0,b.size()/4).Zero();
  b.SubVector(3*b.size()/4,b.size()).Zero();
  DoTestMV1a<Ta,T>(a,b,label);
  DoTestVM1a<Ta,T>(a,b,label);
  DoTestMV1a<Ta,T>(Conjugate(a),b,label+" Conj");
  DoTestVM1a<Ta,T>(Conjugate(a),b,label+" Conj");

  if (b.size() > 1) {
    b = b0;
    b.SubVector(0,1).Zero();
    DoTestMV1a<Ta,T>(a.View(),b,label);
    DoTestVM1a<Ta,T>(a.View(),b,label);
    DoTestMV1a<Ta,T>(Conjugate(a),b,label+" Conj");
    DoTestVM1a<Ta,T>(Conjugate(a),b,label+" Conj");

    b = b0;
    b.SubVector(b.size()-1,b.size()).Zero();
    DoTestMV1a<Ta,T>(a.View(),b,label);
    DoTestVM1a<Ta,T>(a.View(),b,label);
    DoTestMV1a<Ta,T>(Conjugate(a),b,label+" Conj");
    DoTestVM1a<Ta,T>(Conjugate(a),b,label+" Conj");

    b = b0;
    b.SubVector(0,1).Zero();
    b.SubVector(b.size()-1,b.size()).Zero();
    DoTestMV1a<Ta,T>(a.View(),b,label);
    DoTestVM1a<Ta,T>(a.View(),b,label);
    DoTestMV1a<Ta,T>(Conjugate(a),b,label+" Conj");
    DoTestVM1a<Ta,T>(Conjugate(a),b,label+" Conj");
  }
  b=b0;
#endif
}

template <class Ta, class T, class M, class V> inline void DoTestMV2a(
    const M& a, const V& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV2a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" "<<b<<std::endl;
  }

  tmv::Matrix<Ta> m = a;
  tmv::Vector<T> v = b;

  double eps = EPS * Norm(v) * Norm(m);
  double eps2 = EPS * Norm(v) * (1.+Norm(m));

  if (CanMult(m,ColVectorViewOf(v),ColVectorViewOf(v))) {
    if (XXDEBUG5) {
      std::cout<<"CanMult("<<tmv::Type(v)<<","<<tmv::Type(m)<<","<<tmv::Type(v)<<")\n";
    }
    tmv::Vector<T> v0 = b;
    b = a*v0;
    v = m*v0;
    Assert(Norm(b-v) <= eps,label+" v=a*v0");
    v = b = v0;
    b = T(10)*a*v0;
    v = T(10)*m*v0;
    Assert(Norm(b-v) <= eps*10.,label+" v=x*a*v0");
    v = b = v0;
    b += a*v0;
    v = v0 + m*v0;
    Assert(Norm(b-v) <= eps2,label+" v+=a*v0");
    v = b = v0;
    b += T(10)*a*v0;
    v = v0 + T(10)*m*v0;
    Assert(Norm(b-v) <= eps2*10.,label+" v+=x*a*v0");
    v = b = v0;
#ifdef XTEST
    b = -a*v0;
    v = -m*v0;
    Assert(Norm(b-v) <= eps,label+" v=-a*v0");
    v = b = v0;
    b = -T(10)*a*v0;
    v = -T(10)*m*v0;
    Assert(Norm(b-v) <= eps*10.,label+" v=-x*a*v0");
    v = b = v0;
    b += -a*v0;
    v = v0 - m*v0;
    Assert(Norm(b-v) <= eps2,label+" v+=-a*v0");
    v = b = v0;
    b -= a*v0;
    v = v0 - m*v0;
    Assert(Norm(b-v) <= eps2,label+" v-=a*v0");
    v = b = v0;
    b -= -a*v0;
    v = v0 + m*v0;
    Assert(Norm(b-v) <= eps2,label+" v-=-a*v0");
    v = b = v0;
    b += -T(10)*a*v0;
    v = v0 - T(10)*m*v0;
    Assert(Norm(b-v) <= eps2*10.,label+" v+=-x*a*v0");
    v = b = v0;
    b -= T(10)*a*v0;
    v = v0 - T(10)*m*v0;
    Assert(Norm(b-v) <= eps2*10.,label+" v-=x*a*v0");
    v = b = v0;
    b -= -T(10)*a*v0;
    v = v0 + T(10)*m*v0;
    Assert(Norm(b-v) <= eps2*10.,label+" v-=-x*a*v0");
    v = b = v0;
#endif

#ifdef ALIASOK
    b = a*b;
    v = m*v0;
    Assert(Norm(b-v) <= eps,label+" v=a*v");
    v = b = v0;
    b = T(10)*a*b;
    v = T(10)*m*v0;
    Assert(Norm(b-v) <= eps*10.,label+" v=x*a*v");
    v = b = v0;
    b += a*b;
    v = v0 + m*v0;
    Assert(Norm(b-v) <= eps2,label+" v+=a*v");
    v = b = v0;
    b += T(10)*a*b;
    v = v0 + T(10)*m*v0;
    Assert(Norm(b-v) <= eps2*10.,label+" v+=x*a*v");
    v = b = v0;
#ifdef XTEST
    b = -a*b;
    v = -m*v0;
    Assert(Norm(b-v) <= eps,label+" v=-a*v");
    v = b = v0;
    b = -T(10)*a*b;
    v = -T(10)*m*v0;
    Assert(Norm(b-v) <= eps*10.,label+" v=-x*a*v");
    v = b = v0;
    b += -a*b;
    v = v0 - m*v0;
    Assert(Norm(b-v) <= eps2,label+" v+=-a*v");
    v = b = v0;
    b -= a*b;
    v = v0 - m*v0;
    Assert(Norm(b-v) <= eps2,label+" v-=a*v");
    v = b = v0;
    b -= -a*b;
    v = v0 + m*v0;
    Assert(Norm(b-v) <= eps2,label+" v-=-a*v");
    v = b = v0;
    b += -T(10)*a*b;
    v = v0 - T(10)*m*v0;
    Assert(Norm(b-v) <= eps2*10.,label+" v+=-x*a*v");
    v = b = v0;
    b -= T(10)*a*b;
    v = v0 - T(10)*m*v0;
    Assert(Norm(b-v) <= eps2*10.,label+" v-=x*a*v");
    v = b = v0;
    b -= -T(10)*a*b;
    v = v0 + T(10)*m*v0;
    Assert(Norm(b-v) <= eps2*10.,label+" v-=-x*a*v");
    v = b = v0;
#endif
#endif // ALIAS
  }
  if (CanMult(RowVectorViewOf(v),m,RowVectorViewOf(v))) {
    if (XXDEBUG5) {
      std::cout<<"CanMult("<<tmv::Type(v)<<","<<tmv::Type(m)<<","<<tmv::Type(v)<<")\n";
    }
    tmv::Vector<T> v0 = b;
    b = v0*a;
    v = v0*m;
    Assert(Norm(b-v) <= eps,label+" v=v0*a");
    v = b = v0;
    b = T(10)*v0*a;
    v = T(10)*v0*m;
    Assert(Norm(b-v) <= eps*10.,label+" v=x*v0*a");
    v = b = v0;
    b *= a;
    v = v0*m;
    Assert(Norm(b-v) <= eps,label+" v*=a");
    v = b = v0;
    b *= T(10)*a;
    v = T(10)*v0*m;
    Assert(Norm(b-v) <= eps*10.,label+" v*=(x*a)");
    v = b = v0;
    b += v0*a;
    v = v0 + v0*m;
    Assert(Norm(b-v) <= eps2,label+" v+=v0*a");
    v = b = v0;
    b += T(10)*v0*a;
    v = v0 + T(10)*v0*m;
    Assert(Norm(b-v) <= eps2*10.,label+" v+=x*v0*a");
    v = b = v0;
#ifdef XTEST
    b = -v0*a;
    v = -v0*m;
    Assert(Norm(b-v) <= eps,label+" v=-v0*a");
    v = b = v0;
    b = -T(10)*v0*a;
    v = -T(10)*v0*m;
    Assert(Norm(b-v) <= eps*10.,label+" v=-x*v0*a");
    v = b = v0;
    b *= -a;
    v = -v0*m;
    Assert(Norm(b-v) <= eps,label+" v*=-a");
    v = b = v0;
    b *= -T(10)*a;
    v = -T(10)*v0*m;
    Assert(Norm(b-v) <= eps*10.,label+" v*=-(x*a)");
    v = b = v0;
    b += -v0*a;
    v = v0 - v0*m;
    Assert(Norm(b-v) <= eps2,label+" v+=-v0*a");
    v = b = v0;
    b -= v0*a;
    v = v0 - v0*m;
    Assert(Norm(b-v) <= eps2,label+" v-=v0*a");
    v = b = v0;
    b -= -v0*a;
    v = v0 + v0*m;
    Assert(Norm(b-v) <= eps2,label+" v-=-v0*a");
    v = b = v0;
    b += -T(10)*v0*a;
    v = v0 - T(10)*v0*m;
    Assert(Norm(b-v) <= eps2*10.,label+" v+=-x*v0*a");
    v = b = v0;
    b -= T(10)*v0*a;
    v = v0 - T(10)*v0*m;
    Assert(Norm(b-v) <= eps2*10.,label+" v-=x*v0*a");
    v = b = v0;
    b -= -T(10)*v0*a;
    v = v0 + T(10)*v0*m;
    Assert(Norm(b-v) <= eps2*10.,label+" v-=-x*v0*a");
    v = b = v0;
#endif

#ifdef ALIASOK
    b = b*a;
    v = v0*m;
    Assert(Norm(b-v) <= eps,label+" v=v*a");
    v = b = v0;
    b = T(10)*b*a;
    v = T(10)*v0*m;
    Assert(Norm(b-v) <= eps*10.,label+" v=x*v*a");
    v = b = v0;
    b += b*a;
    v = v0 + v0*m;
    Assert(Norm(b-v) <= eps2,label+" v+=v*a");
    v = b = v0;
    b += T(10)*b*a;
    v = v0 + T(10)*v0*m;
    Assert(Norm(b-v) <= eps2*10.,label+" v+=x*v*a");
    v = b = v0;
#ifdef XTEST
    b = -b*a;
    v = -v0*m;
    Assert(Norm(b-v) <= eps,label+" v=-v*a");
    v = b = v0;
    b = -T(10)*b*a;
    v = -T(10)*v0*m;
    Assert(Norm(b-v) <= eps*10.,label+" v=-x*v*a");
    v = b = v0;
    b += -b*a;
    v = v0 - v0*m;
    Assert(Norm(b-v) <= eps2,label+" v+=-v*a");
    v = b = v0;
    b -= b*a;
    v = v0 - v0*m;
    Assert(Norm(b-v) <= eps2,label+" v-=v*a");
    v = b = v0;
    b -= -b*a;
    v = v0 + v0*m;
    Assert(Norm(b-v) <= eps2,label+" v-=-v*a");
    v = b = v0;
    b += -T(10)*b*a;
    v = v0 - T(10)*v0*m;
    Assert(Norm(b-v) <= eps2*10.,label+" v+=-x*v*a");
    v = b = v0;
    b -= T(10)*b*a;
    v = v0 - T(10)*v0*m;
    Assert(Norm(b-v) <= eps2*10.,label+" v-=x*v*a");
    v = b = v0;
    b -= -T(10)*b*a;
    v = v0 + T(10)*v0*m;
    Assert(Norm(b-v) <= eps2*10.,label+" v-=-x*v*a");
    v = b = v0;
#endif
#endif
  }
  if (showstartdone)
    std::cout<<"Done MV2a"<<std::endl;
}

template <class Ta, class T, class M, class V> inline void DoTestMV2R(
    const M& a, const V& b, std::string label)
{
  DoTestMV2a<Ta,T>(a,b,label);
  DoTestMV2a<Ta,T>(a,b.Reverse(),label+" Rev");

#ifdef XTEST
  tmv::Vector<T> b0 = b;
  b.Zero();
  DoTestMV2a<Ta,T>(a,b,label+" 1");
  b = b0;

  b.SubVector(0,b.size()/2).Zero();
  DoTestMV2a<Ta,T>(a,b,label+" 2");
  b = b0;

  b.SubVector(b.size()/2,b.size()).Zero();
  DoTestMV2a<Ta,T>(a,b,label+" 3");
  b = b0;

  b.SubVector(0,b.size()/4).Zero();
  b.SubVector(3*b.size()/4,b.size()).Zero();
  DoTestMV2a<Ta,T>(a,b,label+" 4");
  b = b0;
#endif
}

template <class Ta, class T, class M, class V> inline void DoTestMV2C(
    const M& a, const V& b, std::string label)
{
  DoTestMV2a<Ta,T>(a,b,label);
  DoTestMV2a<Ta,T>(a,b.Reverse(),label+" Rev");

  DoTestMV2a<Ta,T>(Conjugate(a),b,label+" Conj");
#ifdef XTEST
  DoTestMV2a<Ta,T>(Conjugate(a),b.Reverse(),label+" Rev,Conj");

  tmv::Vector<T> b0 = b;
  b.Zero();
  DoTestMV2a<Ta,T>(a,b,label+" 1");
  DoTestMV2a<Ta,T>(Conjugate(a),b,label+" Conj1");
  b = b0;

  b.SubVector(0,b.size()/2).Zero();
  DoTestMV2a<Ta,T>(a,b,label+" 2");
  DoTestMV2a<Ta,T>(Conjugate(a),b,label+" Conj2");
  b = b0;

  b.SubVector(b.size()/2,b.size()).Zero();
  DoTestMV2a<Ta,T>(a,b,label+" 3");
  DoTestMV2a<Ta,T>(Conjugate(a),b,label+" Conj3");
  b = b0;

  b.SubVector(0,b.size()/4).Zero();
  b.SubVector(3*b.size()/4,b.size()).Zero();
  DoTestMV2a<Ta,T>(a,b,label+" 4");
  DoTestMV2a<Ta,T>(Conjugate(a),b,label+" Conj4");
  b = b0;
#endif
}

template <class Ta, class Tb, class T, class M, class V1, class V2> 
inline void DoTestMV3a(
    const M& a, const V1& b, const V2& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV3a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::Type(c)<<" "<<c<<std::endl;
  }

  tmv::Matrix<Ta> m = a;
  tmv::Matrix<Ta> mt = Transpose(a);
  tmv::Vector<T> c0 = c;
  tmv::Vector<Tb> v1 = b;
  tmv::Vector<T> v2 = c;

  double eps = EPS * Norm(b) * Norm(a);
  double eps2 = EPS * (Norm(c0) + Norm(b) * Norm(a));

  if (XXDEBUG6) {
    std::cout<<"a = "<<Type(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<Type(b)<<"  "<<b.step()<<"  "<<b<<std::endl;
    std::cout<<"c = "<<Type(c)<<"  "<<c.step()<<"  "<<c<<std::endl;
    std::cout<<"a*b = "<<a*b<<std::endl;
  }

  c = a*b;
  v2 = m*v1;
  Assert(Norm(c-v2) <= eps,label+" c=a*b");
  c = v2 = c0;
  c += a*b;
  v2 = c0 + m*v1;
  Assert(Norm(c-v2) <= eps2,label+" c+=a*b");
  c = v2 = c0;
  c += T(10)*a*b;
  v2 = c0 + T(10)*m*v1;
  Assert(Norm(c-v2) <= eps2*10.,label+" c+=x*a*b");
  c = v2 = c0;
#ifdef XTEST
  c = -a*b;
  v2 = -m*v1;
  Assert(Norm(c-v2) <= eps,label+" c=-a*b");
  c = v2 = c0;
  c += -a*b;
  v2 = c0 - m*v1;
  Assert(Norm(c-v2) <= eps2,label+" c+=-a*b");
  c = v2 = c0;
  c -= a*b;
  v2 = c0 - m*v1;
  Assert(Norm(c-v2) <= eps2,label+" c-=a*b");
  c = v2 = c0;
  c -= -a*b;
  v2 = c0 + m*v1;
  Assert(Norm(c-v2) <= eps2,label+" c-=-a*b");
  c = v2 = c0;
  c += -T(10)*a*b;
  v2 = c0 - T(10)*m*v1;
  Assert(Norm(c-v2) <= eps2*10.,label+" c+=-x*a*b");
  c = v2 = c0;
  c -= T(10)*a*b;
  v2 = c0 - T(10)*m*v1;
  Assert(Norm(c-v2) <= eps2*10.,label+" c-=x*a*b");
  c = v2 = c0;
  c -= -T(10)*a*b;
  v2 = c0 + T(10)*m*v1;
  Assert(Norm(c-v2) <= eps2*10.,label+" c-=-x*a*b");
#endif

  c = b*Transpose(a);
  v2 = v1*mt;
  Assert(Norm(c-v2) <= eps,label+" c=b*a");
  c = v2 = c0;
  c += b*Transpose(a);
  v2 = c0 + v1*mt;
  Assert(Norm(c-v2) <= eps2,label+" c+=b*a");
  c = v2 = c0;
  c += T(10)*b*Transpose(a);
  v2 = c0 + T(10)*v1*mt;
  Assert(Norm(c-v2) <= eps2*10.,label+" c+=x*b*a");
  c = v2 = c0;
#ifdef XTEST
  c = -b*Transpose(a);
  v2 = -v1*mt;
  Assert(Norm(c-v2) <= eps,label+" c=-b*a");
  c = v2 = c0;
  c += -b*Transpose(a);
  v2 = c0 - v1*mt;
  Assert(Norm(c-v2) <= eps2,label+" c+=-b*a");
  c = v2 = c0;
  c -= b*Transpose(a);
  v2 = c0 - v1*mt;
  Assert(Norm(c-v2) <= eps2,label+" c-=b*a");
  c = v2 = c0;
  c -= -b*Transpose(a);
  v2 = c0 + v1*mt;
  Assert(Norm(c-v2) <= eps2,label+" c-=-b*a");
  c = v2 = c0;
  c += -T(10)*b*Transpose(a);
  v2 = c0 - T(10)*v1*mt;
  Assert(Norm(c-v2) <= eps2*10.,label+" c+=-x*b*a");
  c = v2 = c0;
  c -= T(10)*b*Transpose(a);
  v2 = c0 - T(10)*v1*mt;
  Assert(Norm(c-v2) <= eps2*10.,label+" c-=x*b*a");
  c = v2 = c0;
  c -= -T(10)*b*Transpose(a);
  v2 = c0 + T(10)*v1*mt;
  Assert(Norm(c-v2) <= eps2*10.,label+" c-=-x*b*a");
#endif

  if (showstartdone)
    std::cout<<"Done MV3a"<<std::endl;
}

template <class Ta, class Tb, class T, class M, class V1, class V2> 
inline void DoTestMV3R(
    const M& a, const V1& b, const V2& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV3"<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::Type(c)<<" "<<c<<std::endl;
  }

  DoTestMV3a<Ta,Tb,T>(a,b,c,label);
  DoTestMV3a<Ta,Tb,T>(a,b.Reverse(),c,label);
  DoTestMV3a<Ta,Tb,T>(a,b,c.Reverse(),label);
  DoTestMV3a<Ta,Tb,T>(a,b.Reverse(),c.Reverse(),label);

  if (showstartdone)
    std::cout<<"Done MV3"<<std::endl;
}

template <class Ta, class Tb, class T, class M, class V1, class V2> 
inline void DoTestMV3C(
    const M& a, const V1& b, const V2& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV3"<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::Type(c)<<" "<<c<<std::endl;
  }

  DoTestMV3a<Ta,Tb,T>(a,b,c,label);
  DoTestMV3a<Ta,Tb,T>(a,b.Reverse(),c,label);
  DoTestMV3a<Ta,Tb,T>(a,b,c.Reverse(),label);
  DoTestMV3a<Ta,Tb,T>(a,b.Reverse(),c.Reverse(),label);

  DoTestMV3a<Ta,Tb,T>(Conjugate(a),b,c,label+" Conj");
#ifdef XTEST
  DoTestMV3a<Ta,Tb,T>(Conjugate(a),b.Reverse(),c,label+" Conj");
  DoTestMV3a<Ta,Tb,T>(Conjugate(a),b,c.Reverse(),label+" Conj");
  DoTestMV3a<Ta,Tb,T>(Conjugate(a),b.Reverse(),c.Reverse(),label+" Conj");
#endif

  if (showstartdone)
    std::cout<<"Done MV3"<<std::endl;
}

template <class T, class Tb, class Tsum, class M1, class M2> inline void DoTestMM1a(
    const M1& a, const M2& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM1a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" "<<b<<std::endl;
  }

  if (XXDEBUG7) {
    std::cout<<"a = "<<tmv::Type(a)<<" = "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" = "<<b<<std::endl;
  }

  if (CanAdd(a,b)) {
    tmv::Matrix<T> m1 = a;
    tmv::Matrix<Tb> m2 = b;
    tmv::Matrix<Tsum> sum = m1+m2;
    tmv::Matrix<Tsum> diff = m1-m2;

    double eps = EPS*(Norm(m1)+Norm(m2));
    if (XXDEBUG7) {
      std::cout<<"CanAdd("<<tmv::Type(a)<<","<<tmv::Type(b)<<")\n";
      std::cout<<"m1-m2 = "<<m1-m2<<std::endl;
      std::cout<<"a-b = "<<a-b<<std::endl;
      std::cout<<"m1+m2 = "<<m1+m2<<std::endl;
      std::cout<<"a+b = "<<a+b<<std::endl;
    }
    Assert(Norm((a-m2)-diff) <= eps,label+" a-m");
    Assert(Norm((m1-b)-diff) <= eps,label+" m-b");
    Assert(Norm((a-b)-diff) <= eps,label+" a-b");
    Assert(Norm((a+m2)-sum) <= eps,label+" a+m");
    Assert(Norm((m1+b)-sum) <= eps,label+" m+b");
    Assert(Norm((a+b)-sum) <= eps,label+" a+b");
  }
  if (showstartdone)
    std::cout<<"Done MM1a"<<std::endl;
}

template <class T, class M1, class M2> inline void DoTestMM1RR(
    const M1& a, const M2& b, std::string label)
{
  DoTestMM1a<T,T,T>(a,b,label);

#ifdef XTEST
  DoTestMM1a<T,T,T>(Transpose(b),Transpose(a),label+" TransB TransA");
#endif
}

template <class T, class M1, class M2> inline void DoTestMM1RC(
    const M1& a, const M2& b, std::string label)
{
  DoTestMM1a<T,CT,CT>(a,b,label);

#ifdef XTEST
  DoTestMM1a<CT,T,CT>(Transpose(b),Transpose(a),
      label+" TransB TransA");

  DoTestMM1a<T,CT,CT>(a,Conjugate(b),label+" ConjB");
  DoTestMM1a<CT,T,CT>(Adjoint(b),Transpose(a),label+" AdjB TransA");
#endif
}

template <class T, class M1, class M2> inline void DoTestMM1CR(
    const M1& a, const M2& b, std::string label)
{
  DoTestMM1a<CT,T,CT>(a,b,label);

#ifdef XTEST
  DoTestMM1a<T,CT,CT>(Transpose(b),Transpose(a),label+" TransB TransA");

  DoTestMM1a<CT,T,CT>(Conjugate(a),b,label+" ConjA");
  DoTestMM1a<T,CT,CT>(Transpose(b),Adjoint(a),label+" TransB AdjA");
#endif
}

template <class T, class M1, class M2> inline void DoTestMM1CC(
    const M1& a, const M2& b, std::string label)
{
  DoTestMM1a<CT,CT,CT>(a,b,label);

#ifdef XTEST
  DoTestMM1a<CT,CT,CT>(Transpose(b),Transpose(a),label+" TransB TransA");

  DoTestMM1a<CT,CT,CT>(Conjugate(a),b,label+" ConjA");
  DoTestMM1a<CT,CT,CT>(Transpose(b),Adjoint(a),label+" TransB AdjA");

  DoTestMM1a<CT,CT,CT>(a,Conjugate(b),label+" ConjB");
  DoTestMM1a<CT,CT,CT>(Adjoint(b),Transpose(a),label+" AdjB TransA");

  DoTestMM1a<CT,CT,CT>(Conjugate(a),Conjugate(b),label+" ConjA ConjB");
  DoTestMM1a<CT,CT,CT>(Adjoint(b),Adjoint(a),label+" AdjA ConjB");
#endif
}

template <class T, class Tb, class BaseM, class M1, class M2> inline void DoTestMM2a(
    BaseM& a0, const M1& a, const M2& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM2a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"a0 = "<<tmv::Type(a0)<<" "<<a0<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" "<<b<<std::endl;
  }

  if (XXDEBUG8) {
    std::cout<<"a = "<<tmv::Type(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<"  "<<b<<std::endl;
  }
  a0 = a;
  const tmv::Matrix<T> m1 = a;
  const tmv::Matrix<Tb> m2 = b;
  if (XXDEBUG8) {
    std::cout<<"m1 = "<<tmv::Type(m1)<<"  "<<m1<<std::endl;
    std::cout<<"m2 = "<<tmv::Type(m2)<<"  "<<m2<<std::endl;
  }
  double eps = EPS*(Norm(m1)+Norm(m2));

  if (CanAddEq(m1,b)) {
    if (XXDEBUG8) {
      std::cout<<"CanAddEq("<<tmv::Type(m1)<<","<<tmv::Type(b)<<")\n";
    }
    tmv::Matrix<T> m3 = m1;
    tmv::Matrix<T> m4 = m1;
    m3 += b;
    m4 = m1+m2;
    Assert(Norm(m3-m4) <= eps,label+" m += b");
    m3 = m1;
    m3 = m3+b;
    m4 = m1+m2;
    Assert(Norm(m3-m4) <= eps,label+" m = m+b");
    m3 = m1;
#ifdef XTEST
    m3 += -b;
    m4 = m1-m2;
    Assert(Norm(m3-m4) <= eps,label+" m += -b");
    m3 = m1;
    m3 -= b;
    m4 = m1-m2;
    Assert(Norm(m3-m4) <= eps,label+" m -= b");
    m3 = m1;
    m3 -= -b;
    m4 = m1+m2;
    Assert(Norm(m3-m4) <= eps,label+" m -= -b");
    m3 = m1;
    m3 = m3-b;
    m4 = m1-m2;
    Assert(Norm(m3-m4) <= eps,label+" m = m-b");
    m3 = m1;
    m3 = b+m3;
    m4 = m2+m1;
    Assert(Norm(m3-m4) <= eps,label+" m = b+m");
    m3 = m1;
    m3 = b-m3;
    m4 = m2-m1;
    Assert(Norm(m3-m4) <= eps,label+" m = b-m");
#endif
  }
#ifndef NOADDEQ
  if (CanAddEq(a,b)) {
    if (XXDEBUG8) {
      std::cout<<"CanAddEq("<<tmv::Type(a)<<","<<tmv::Type(b)<<")\n";
    }
    tmv::Matrix<T> m4 = a = a0;
    a += b;
    m4 = m1+m2;
    if (XXDEBUG8) {
      std::cout<<"a += b = "<<a<<std::endl;
      std::cout<<"m4 = "<<m4<<std::endl;
    }
    Assert(Norm(a-m4) <= eps,label+" a += b");
    a = a0;
#ifdef XTEST
    a += -b;
    m4 = m1-m2;
    Assert(Norm(a-m4) <= eps,label+" a += -b");
    a = a0;
    a -= b;
    m4 = m1-m2;
    Assert(Norm(a-m4) <= eps,label+" a -= b");
    a = a0;
    a -= -b;
    m4 = m1+m2;
    Assert(Norm(a-m4) <= eps,label+" a -= -b");
    a = a0;
#endif
#ifdef ALIASOK
    a = a+b;
    m4 = m1+m2;
    Assert(Norm(a-m4) <= eps,label+" a = a+b");
    a = a0;
#ifdef XTEST
    a = a-b;
    m4 = m1-m2;
    Assert(Norm(a-m4) <= eps,label+" a = a-b");
    a = a0;
    a = b+a; 
    m4 = m2+m1;
    Assert(Norm(a-m4) <= eps,label+" a = b+a");
    a = a0;
    a = b-a;
    m4 = m2-m1;
    Assert(Norm(a-m4) <= eps,label+" a = b-a");
    a = a0;
#endif
#endif
  }
#endif // NOADDEQ

  if (showstartdone)
    std::cout<<"Done MM2a"<<std::endl;
}

template <class T, class BaseM, class M1, class M2> inline void DoTestMM2RR(
    BaseM& a0, const M1& a, const M2& b, std::string label)
{
  DoTestMM2a<T,T>(a0,a,b,label);
}

template <class T, class BaseM, class M1, class M2> inline void DoTestMM2RC(
    BaseM& a0, const M1& a, const M2& b, std::string label)
{
  DoTestMM2a<T,CT>(a0,a,b,label);
#ifdef XTEST
  DoTestMM2a<T,CT>(a0,a,Conjugate(b),label+" ConjB");
#endif
}

template <class T, class BaseM, class M1, class M2> inline void DoTestMM2CR(
    BaseM& a0, const M1& a, const M2& b, std::string label)
{
  DoTestMM2a<CT,T>(a0,a,b,label);
  DoTestMM2a<CT,T>(a0,Conjugate(a),b,label+" ConjA");
}

template <class T, class BaseM, class M1, class M2> inline void DoTestMM2CC(
    BaseM& a0, const M1& a, const M2& b, std::string label)
{
  DoTestMM2a<CT,CT>(a0,a,b,label);
  DoTestMM2a<CT,CT>(a0,Conjugate(a),b,label+" ConjA");

#ifdef XTEST
  DoTestMM2a<CT,CT>(a0,a,Conjugate(b),label+" ConjB");
  DoTestMM2a<CT,CT>(a0,Conjugate(a),Conjugate(b),label+" ConjA ConjB");
#endif
}

template <class T, class Tb, class Tprod, class M1, class M2> inline void DoTestMM3a(
    const M1& a, const M2& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM3a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" "<<b<<std::endl;
  }

  if (XXDEBUG7) {
    std::cout<<"a = "<<tmv::Type(a)<<" = "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" = "<<b<<std::endl;
  }

  if (CanMult(a,b)) {
    tmv::Matrix<T> m1 = a;
    tmv::Matrix<Tb> m2 = b;
    tmv::Matrix<Tprod> mm = m1*m2;
    double eps = EPS*Norm(m1)*Norm(m2);
    if (XXDEBUG7) {
      std::cout<<"CanMult("<<tmv::Type(a)<<","<<tmv::Type(b)<<")\n";
      std::cout<<"m1*m2 = "<<mm<<std::endl;
      std::cout<<"m1*b = "<<m1*b<<std::endl;
      std::cout<<"a*m2 = "<<a*m2<<std::endl;
      std::cout<<"a*b = "<<a*b<<std::endl;
    }
    Assert(Norm((m1*b)-mm) <= eps,label+" m*b");
    Assert(Norm((a*m2)-mm) <= eps,label+" a*m");
    Assert(Norm((a*b)-mm) <= eps,label+" a*b");
  }

  if (showstartdone)
    std::cout<<"Done MM3a"<<std::endl;
}

template <class T, class M1, class M2> inline void DoTestMM3RR(
    const M1& a, const M2& b, std::string label)
{
  DoTestMM3a<T,T,T>(a,b,label);

#ifdef XTEST
  DoTestMM3a<T,T,T>(Transpose(b),Transpose(a),label+" TransB TransA");
#endif
}

template <class T, class M1, class M2> inline void DoTestMM3RC(
    const M1& a, const M2& b, std::string label)
{
  DoTestMM3a<T,CT,CT>(a,b,label);

#ifdef XTEST
  DoTestMM3a<CT,T,CT>(Transpose(b),Transpose(a),label+" TransB TransA");

  DoTestMM3a<T,CT,CT>(a,Conjugate(b),label+" ConjB");
  DoTestMM3a<CT,T,CT>(Adjoint(b),Transpose(a),label+" AdjB TransA");
#endif
}

template <class T, class M1, class M2> inline void DoTestMM3CR(
    const M1& a, const M2& b, std::string label)
{
  DoTestMM3a<CT,T,CT>(a,b,label);

#ifdef XTEST
  DoTestMM3a<T,CT,CT>(Transpose(b),Transpose(a),label+" TransB TransA");

  DoTestMM3a<CT,T,CT>(Conjugate(a),b,label+" ConjA");
  DoTestMM3a<T,CT,CT>(Transpose(b),Adjoint(a),label+" TransB AdjA");
#endif
}

template <class T, class M1, class M2> inline void DoTestMM3CC(
    const M1& a, const M2& b, std::string label)
{
  DoTestMM3a<CT,CT,CT>(a,b,label);

#ifdef XTEST
  DoTestMM3a<CT,CT,CT>(Transpose(b),Transpose(a),label+" TransB TransA");

  DoTestMM3a<CT,CT,CT>(Conjugate(a),b,label+" ConjA");
  DoTestMM3a<CT,CT,CT>(Transpose(b),Adjoint(a),label+" TransB AdjA");

  DoTestMM3a<CT,CT,CT>(a,Conjugate(b),label+" ConjB");
  DoTestMM3a<CT,CT,CT>(Adjoint(b),Transpose(a),label+" AdjB TransA");

  DoTestMM3a<CT,CT,CT>(Conjugate(a),Conjugate(b),label+" ConjA ConjB");
  DoTestMM3a<CT,CT,CT>(Adjoint(b),Adjoint(a),label+" AdjA ConjB");
#endif
}

template <class T, class Tb, class Tprod, class BaseM, class M1, class M2> inline void DoTestMM4a(
    BaseM& a0, const M1& a, const M2& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM4a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"a0 = "<<tmv::Type(a0)<<" "<<a0<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" "<<b<<std::endl;
  }

  if (XXDEBUG8) {
    std::cout<<"a = "<<tmv::Type(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<"  "<<b<<std::endl;
  }
  a0 = a;
  const tmv::Matrix<T> m1 = a;
  const tmv::Matrix<Tb> m2 = b;
  if (XXDEBUG8) {
    std::cout<<"m1 = "<<tmv::Type(m1)<<"  "<<m1<<std::endl;
    std::cout<<"m2 = "<<tmv::Type(m2)<<"  "<<m2<<std::endl;
  }

  double eps = EPS*Norm(m1)*Norm(m2);
  double eps2 = EPS*Norm(m1)*(1.+Norm(m2));

  if (m1.rowsize() == m2.colsize()) {
    tmv::Matrix<Tprod> mm = m1*m2;
    if (CanMult(m1,b,m1)) {
      if (XXDEBUG8) {
	std::cout<<"CanMult("<<tmv::Type(m1)<<","<<tmv::Type(b)<<","<<tmv::Type(m1)<<")\n";
      }
      tmv::Matrix<T> m3 = m1;
      tmv::Matrix<T> m4 = m1;
      m3 *= b;
      m4 = mm;
      Assert(Norm(m3-m4) <= eps,label+" m *= b");
      m3 = m1;
      m3 *= T(10)*b;
      m4 = T(10)*mm;
      Assert(Norm(m3-m4) <= eps*10.,label+" m *= x*b");
      m3 = m1;
      m3 += m1*b;
      m4 = m1 + mm;
      Assert(Norm(m3-m4) <= eps2,label+" m += m1*b");
      m3 = m1;
      m3 += T(10)*m1*b;
      m4 = m1 + T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m += x*m1*b");
      m3 = m1;
      m3 += m3*b;
      m4 = m1 + mm;
      Assert(Norm(m3-m4) <= eps2,label+" m += m*b");
      m3 = m1;
      m3 = a*b;
      m4 = mm;
      Assert(Norm(m3-m4) <= eps,label+" m = a*b");
      m3 = m1;
      m3 += T(10)*m3*b;
      m4 = m1 + T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m += x*m*b");
      m3 = m1;
      m3 += a*b;
      m4 = m1 + mm;
      Assert(Norm(m3-m4) <= eps2,label+" m += a*b");
      m3 = m1;
      m3 = T(10)*a*b;
      m4 = T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m = x*a*b");
      m3 = m1;
      m3 += T(10)*a*b;
      m4 = m1 + T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m += x*a*b");
      m3 = m1;
#ifdef XTEST
      m3 *= -b;
      m4 = -mm;
      Assert(Norm(m3-m4) <= eps,label+" m *= -b");
      m3 = m1;
      m3 *= -T(10)*b;
      m4 = -T(10)*mm;
      Assert(Norm(m3-m4) <= eps*10.,label+" m *= -x*b");
      m3 = m1;
      m3 *= Tb(10)*b;
      m4 = Tb(10)*mm;
      Assert(Norm(m3-m4) <= eps*10.,label+" m *= x*b");
      m3 = m1;
      m3 *= -Tb(10)*b;
      m4 = -Tb(10)*mm;
      Assert(Norm(m3-m4) <= eps*10.,label+" m *= -x*b");
      m3 = m1;
      m3 = m3*b;
      m4 = mm;
      Assert(Norm(m3-m4) <= eps,label+" m = m*b");
      m3 = m1;
      m3 = -m3*b;
      m4 = -mm;
      Assert(Norm(m3-m4) <= eps,label+" m = -m*b");
      m3 = m1;
      m3 = m3*(T(10)*b);
      m4 = T(10)*mm;
      Assert(Norm(m3-m4) <= eps*10.,label+" m = m*(x*b)");
      m3 = m1;
      m3 = (T(10)*m3)*b;
      m4 = T(10)*mm;
      Assert(Norm(m3-m4) <= eps*10.,label+" m = (x*m)*b");
      m3 = m1;
      m3 = (T(10)*m3)*(T(10)*b);
      m4 = T(100)*mm;
      Assert(Norm(m3-m4) <= eps*100.,label+" m = (x*m)*(x*b)");
      m3 = m1;
      m3 += -m1*b;
      m4 = m1 - mm;
      Assert(Norm(m3-m4) <= eps2,label+" m += -m1*b");
      m3 = m1;
      m3 -= m1*b;
      m4 = m1 - mm;
      Assert(Norm(m3-m4) <= eps2,label+" m -= m1*b");
      m3 = m1;
      m3 -= -m1*b;
      m4 = m1 + mm;
      Assert(Norm(m3-m4) <= eps2,label+" m -= -m1*b");
      m3 = m1;
      m3 -= T(10)*m1*b;
      m4 = m1 - T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m -= x*m1*b");
      m3 = m1;
      m3 += -m3*b;
      m4 = m1 - mm;
      Assert(Norm(m3-m4) <= eps2,label+" m += -m*b");
      m3 = m1;
      m3 -= m3*b;
      m4 = m1 - mm;
      Assert(Norm(m3-m4) <= eps2,label+" m -= m*b");
      m3 = m1;
      m3 -= -m3*b;
      m4 = m1 + mm;
      Assert(Norm(m3-m4) <= eps2,label+" m -= -m*b");
      m3 = m1;
      m3 -= T(10)*m3*b;
      m4 = m1 - T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m -= x*m*b");
      m3 = m1;
      m3 = -a*b;
      m4 = -mm;
      Assert(Norm(m3-m4) <= eps,label+" m = -a*b");
      m3 = m1;
      m3 += -a*b;
      m4 = m1 - mm;
      Assert(Norm(m3-m4) <= eps2,label+" m += -a*b");
      m3 = m1;
      m3 -= a*b;
      m4 = m1 - mm;
      Assert(Norm(m3-m4) <= eps2,label+" m -= a*b");
      m3 = m1;
      m3 -= -a*b;
      m4 = m1 + mm;
      Assert(Norm(m3-m4) <= eps2,label+" m -= -a*b");
      m3 = m1;
      m3 = -T(10)*a*b;
      m4 = -T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m = -x*a*b");
      m3 = m1;
      m3 += -T(10)*a*b;
      m4 = m1 - T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m += -x*a*b");
      m3 = m1;
      m3 -= T(10)*a*b;
      m4 = m1 - T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m -= x*a*b");
      m3 = m1;
      m3 -= -T(10)*a*b;
      m4 = m1 + T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m -= -x*a*b");
#endif
    }
#ifndef NOMULTEQ
    if (CanMult(a,b,a)) {
      if (XXDEBUG8) {
	std::cout<<"CanMult("<<tmv::Type(a)<<","<<tmv::Type(b)<<","<<tmv::Type(a)<<")\n";
      }
      tmv::Matrix<T> m4 = m1;
      a *= b;
      m4 = mm;
      if (XXDEBUG8) {
	std::cout<<"a *= b = "<<a<<std::endl;
	std::cout<<"m4 = "<<m4<<std::endl;
      }
      Assert(Norm(a-m4) <= eps,label+" a *= b");
      a = a0;
#ifdef ALIASOK
      a = a*b;
      m4 = mm;
      Assert(Norm(a-m4) <= eps,label+" a = a*b");
      a = a0;
#endif
    }
    if (CanMultXM(a,b,a)) {
      if (XXDEBUG8) {
	std::cout<<"CanMultXM("<<tmv::Type(a)<<","<<tmv::Type(b)<<","<<tmv::Type(a)<<")\n";
      }
      tmv::Matrix<T> m4 = m1;
      a *= -b;
      m4 = -mm;
      Assert(Norm(a-m4) <= eps,label+" a *= -b");
      a = a0;
      a *= T(10)*b;
      m4 = T(10)*mm;
      Assert(Norm(a-m4) <= eps*10.,label+" a *= x*b");
      a = a0;
#ifdef XTEST
      a *= -T(10)*b;
      m4 = -T(10)*mm;
      Assert(Norm(a-m4) <= eps*10.,label+" a *= -x*b");
      a = a0;
      a *= Tb(10)*b;
      m4 = Tb(10)*mm;
      Assert(Norm(a-m4) <= eps*10.,label+" a *= x*b");
      a = a0;
      a *= -Tb(10)*b;
      m4 = -Tb(10)*mm;
      Assert(Norm(a-m4) <= eps*10.,label+" a *= -x*b");
      a = a0;
#endif
#ifdef ALIASOK
      a = -a*b;
      m4 = -mm;
      Assert(Norm(a-m4) <= eps,label+" a = -a*b");
      a = a0;
      a += a*b;
      m4 = m1 + mm;
      Assert(Norm(a-m4) <= eps2,label+" a += a*b");
      a = a0;
      a = (T(10)*a)*b;
      m4 = T(10)*mm;
      Assert(Norm(a-m4) <= eps*10.,label+" a = (x*a)*b");
      a = a0;
      a += T(10)*a*b;
      m4 = m1 + T(10)*mm;
      Assert(Norm(a-m4) <= eps2*10.,label+" a += x*a*b");
      a = a0;
#ifdef XTEST
      a += -a*b;
      m4 = m1 - mm;
      Assert(Norm(a-m4) <= eps2,label+" a += -a*b");
      a = a0;
      a -= a*b;
      m4 = m1 - mm;
      Assert(Norm(a-m4) <= eps2,label+" a -= a*b");
      a = a0;
      a -= -a*b;
      m4 = m1 + mm;
      Assert(Norm(a-m4) <= eps2,label+" a -= -a*b");
      a = a0;
      a = a*(T(10)*b);
      m4 = T(10)*mm;
      Assert(Norm(a-m4) <= eps*10.,label+" a = a*(x*b)");
      a = a0;
      a = (T(10)*a)*(T(10)*b);
      m4 = T(100)*mm;
      Assert(Norm(a-m4) <= eps*100.,label+" a = (x*a)*(x*b)");
      a = a0;
      a -= T(10)*a*b;
      m4 = m1 - T(10)*mm;
      Assert(Norm(a-m4) <= eps2*10.,label+" a -= x*a*b");
      a = a0;
#endif
#endif
    }
#endif // NOMULTEQ
  }
  if (m2.rowsize() == m1.colsize()) {
    tmv::Matrix<T> mm = m2*m1;
    if (CanMult(b,m1,m1)) {
      if (XXDEBUG8) {
	std::cout<<"CanMult("<<tmv::Type(b)<<","<<tmv::Type(m1)<<","<<tmv::Type(m1)<<")\n";
      }
      tmv::Matrix<T> m3 = m1;
      tmv::Matrix<T> m4 = m1;
      m3 = b * m3;
      m4 = mm;
      Assert(Norm(m3-m4) <= eps,label+" m = b*m");
      m3 = m1;
      m3 = (T(10)*b) * m3;
      m4 = T(10)*mm;
      Assert(Norm(m3-m4) <= eps*10.,label+" m = (x*b)*m");
      m3 = m1;
      m3 += b * m3;
      m4 = m1 + mm;
      Assert(Norm(m3-m4) <= eps2,label+" m += b*m");
      m3 = m1;
      m3 += T(10)*b * m3;
      m4 = m1 + T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m += x*b*m");
      m3 = m1;
      m3 = b * a;
      m4 = mm;
      Assert(Norm(m3-m4) <= eps,label+" m = b*a");
      m3 = m1;
      m3 += b * a;
      m4 = m1 + mm;
      Assert(Norm(m3-m4) <= eps2,label+" m += b*a");
      m3 = m1;
      m3 = T(10)*b*a;
      m4 = T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m = x*b*a");
      m3 = m1;
      m3 += T(10)*b*a;
      m4 = m1 + T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m += x*b*a");
      m3 = m1;
#ifdef XTEST
      m3 = -b * m3;
      m4 = -mm;
      Assert(Norm(m3-m4) <= eps,label+" m = -b*m");
      m3 = m1;
      m3 = b * (T(10)*m3);
      m4 = T(10)*mm;
      Assert(Norm(m3-m4) <= eps*10.,label+" m = b*(x*m)");
      m3 = m1;
      m3 = (T(10)*b) * (T(10)*m3);
      m4 = T(100)*mm;
      Assert(Norm(m3-m4) <= eps*100.,label+" m = (x*b)*(x*m)");
      m3 = m1;
      m3 += -b * m3;
      m4 = m1 - mm;
      Assert(Norm(m3-m4) <= eps2,label+" m += -b*m");
      m3 = m1;
      m3 -= b * m3;
      m4 = m1 - mm;
      Assert(Norm(m3-m4) <= eps2,label+" m -= b*m");
      m3 = m1;
      m3 -= -b * m3;
      m4 = m1 + mm;
      Assert(Norm(m3-m4) <= eps2,label+" m -= -b*m");
      m3 = m1;
      m3 -= T(10)*b * m3;
      m4 = m1 - T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m -= x*b*m");
      m3 = m1;
      m3 = -b * a;
      m4 = -mm;
      Assert(Norm(m3-m4) <= eps,label+" m = -b*a");
      m3 = m1;
      m3 += -b * a;
      m4 = m1 - mm;
      Assert(Norm(m3-m4) <= eps2,label+" m += -b*a");
      m3 = m1;
      m3 -= b * a;
      m4 = m1 - mm;
      Assert(Norm(m3-m4) <= eps2,label+" m -= b*a");
      m3 = m1;
      m3 -= -b * a;
      m4 = m1 + mm;
      Assert(Norm(m3-m4) <= eps2,label+" m -= -b*a");
      m3 = m1;
      m3 = -T(10)*b*a;
      m4 = -T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m = -x*b*a");
      m3 = m1;
      m3 += -T(10)*b*a;
      m4 = m1 - T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m += -x*b*a");
      m3 = m1;
      m3 -= T(10)*b*a;
      m4 = m1 - T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m -= x*b*a");
      m3 = m1;
      m3 -= -T(10)*b*a;
      m4 = m1 + T(10)*mm;
      Assert(Norm(m3-m4) <= eps2*10.,label+" m -= -x*b*a");
#endif
    }
#ifndef NOMULTEQ
#ifdef ALIASOK
    if (CanMult(b,a,a)) {
      if (XXDEBUG8) {
	std::cout<<"CanMult("<<tmv::Type(b)<<","<<tmv::Type(a)<<","<<tmv::Type(a)<<")\n";
      }
      a = b * a;
      Assert(Norm(a-mm) <= eps,label+" a = b*a");
      a = a0;
    }
    if (CanMultXM(b,a,a)) {
      if (XXDEBUG8) {
	std::cout<<"CanMultXM("<<tmv::Type(b)<<","<<tmv::Type(a)<<","<<tmv::Type(a)<<")\n";
      }
      tmv::Matrix<T> m4 = m1;
      a = -b * a;
      m4 = -mm;
      Assert(Norm(a-m4) <= eps,label+" a = -b*a");
      a = a0;
      a += b*a;
      m4 = m1 + mm;
      Assert(Norm(a-m4) <= eps2,label+" a += b*a");
      a = a0;
      a = (T(10)*b) * a;
      m4 = T(10)*mm;
      Assert(Norm(a-m4) <= eps*10.,label+" a = (x*b)*a");
      a = a0;
      a += T(10)*b*a;
      m4 = m1 + T(10)*mm;
      Assert(Norm(a-m4) <= eps2*10.,label+" a += x*b*a");
      a = a0;
#ifdef XTEST
      a -= b*a;
      m4 = m1 - mm;
      Assert(Norm(a-m4) <= eps2,label+" a -= b*a");
      a = a0;
      a = b*(T(10)*a);
      m4 = T(10)*mm;
      Assert(Norm(a-m4) <= eps*10.,label+" a = b*(x*a)");
      a = a0;
      a = (T(10)*b)*(T(10)*a);
      m4 = T(100)*mm;
      Assert(Norm(a-m4) <= eps*100.,label+" a = (x*b)*(x*a)");
      a = a0;
      a -= T(10)*b*a;
      m4 = m1 - T(10)*mm;
      Assert(Norm(a-m4) <= eps2*10.,label+" a -= x*b*a");
      a = a0;
#endif
    }
#endif // ALIASOK
#endif // NOMULTEQ
  }
#ifndef NOMULTEQ
  if (m2.rowsize() == m2.colsize()) {
    eps = EPS*Norm(m2)*Norm(m2);
    eps2 = EPS*Norm(m2)*(1.+Norm(m2));
    tmv::Matrix<Tb> mm = m2*m2;
    if (CanMult(b,b,a)) {
      if (XXDEBUG8) {
	std::cout<<"CanMult("<<tmv::Type(b)<<","<<tmv::Type(b)<<","<<tmv::Type(a)<<")\n";
      }
      tmv::Matrix<T> m4 = m1;
      a = b*b;
      m4 = mm;
      Assert(Norm(a-m4) <= eps,label+" a = b*b");
      a = a0;
    }
    if (CanMultXM(b,b,a)) {
      if (XXDEBUG8) {
	std::cout<<"CanMultXM("<<tmv::Type(b)<<","<<tmv::Type(b)<<","<<tmv::Type(a)<<")\n";
      }
      tmv::Matrix<T> m4 = m1;
      a = -b*b;
      m4 = -mm;
      Assert(Norm(a-m4) <= eps,label+" a = -b*b");
      a = a0;
      a += b*b;
      m4 = m1 + mm;
      Assert(Norm(a-m4) <= eps2,label+" a += b*b");
      a = a0;
      a = T(10)*b*b;
      m4 = T(10)*mm;
      Assert(Norm(a-m4) <= eps*10.,label+" a = x*b*b");
      a = a0;
      a += T(10)*b*b;
      m4 = m1 + T(10)*mm;
      Assert(Norm(a-m4) <= eps2*10.,label+" a += x*b*b");
      a = a0;
#ifdef XTEST
      a += -b*b;
      m4 = m1 - mm;
      Assert(Norm(a-m4) <= eps2,label+" a += -b*b");
      a = a0;
      a -= b*b;
      m4 = m1 - mm;
      Assert(Norm(a-m4) <= eps2,label+" a -= b*b");
      a = a0;
      a -= -b*b;
      m4 = m1 + mm;
      Assert(Norm(a-m4) <= eps2,label+" a -= -b*b");
      a = a0;
      a = -T(10)*b*b;
      m4 = -T(10)*mm;
      Assert(Norm(a-m4) <= eps*10.,label+" a = -x*b*b");
      a = a0;
      a += -T(10)*b*b;
      m4 = m1 - T(10)*mm;
      Assert(Norm(a-m4) <= eps2*10.,label+" a += -x*b*b");
      a = a0;
      a += Tb(10)*b*b;
      m4 = m1 + Tb(10)*mm;
      Assert(Norm(a-m4) <= eps2*10.,label+" a += x*b*b");
      a = a0;
      a += -Tb(10)*b*b;
      m4 = m1 - Tb(10)*mm;
      Assert(Norm(a-m4) <= eps2*10.,label+" a += -x*b*b");
      a = a0;
      a -= T(10)*b*b;
      m4 = m1 - T(10)*mm;
      Assert(Norm(a-m4) <= eps2*10.,label+" a -= x*b*b");
      a = a0;
      a -= -T(10)*b*b;
      m4 = m1 + T(10)*mm;
      Assert(Norm(a-m4) <= eps2*10.,label+" a -= -x*b*b");
      a = a0;
      a -= Tb(10)*b*b;
      m4 = m1 - Tb(10)*mm;
      Assert(Norm(a-m4) <= eps2*10.,label+" a -= x*b*b");
      a = a0;
      a -= -Tb(10)*b*b;
      m4 = m1 + Tb(10)*mm;
      Assert(Norm(a-m4) <= eps2*10.,label+" a -= -x*b*b");
      a = a0;
#endif
    }
  }
#endif // NOMULTEQ

  if (showstartdone)
    std::cout<<"Done MM4a"<<std::endl;
}

template <class T, class BaseM, class M1, class M2> inline void DoTestMM4RR(
    BaseM& a0, const M1& a, const M2& b, std::string label)
{
  DoTestMM4a<T,T,T>(a0,a,b,label);
}

template <class T, class BaseM, class M1, class M2> inline void DoTestMM4RC(
    BaseM& a0, const M1& a, const M2& b, std::string label)
{
  DoTestMM4a<T,CT,CT>(a0,a,b,label);
#ifdef XTEST
  DoTestMM4a<T,CT,CT>(a0,a,Conjugate(b),label+" ConjB");
#endif
}

template <class T, class BaseM, class M1, class M2> inline void DoTestMM4CR(
    BaseM& a0, const M1& a, const M2& b, std::string label)
{
  DoTestMM4a<CT,T,CT>(a0,a,b,label);
  DoTestMM4a<CT,T,CT>(a0,Conjugate(a),b,label+" ConjA");
}

template <class T, class BaseM, class M1, class M2> inline void DoTestMM4CC(
    BaseM& a0, const M1& a, const M2& b, std::string label)
{
  DoTestMM4a<CT,CT,CT>(a0,a,b,label);
  DoTestMM4a<CT,CT,CT>(a0,Conjugate(a),b,label+" ConjA");

#ifdef XTEST
  DoTestMM4a<CT,CT,CT>(a0,a,Conjugate(b),label+" ConjB");
  DoTestMM4a<CT,CT,CT>(a0,Conjugate(a),Conjugate(b),label+" ConjA ConjB");
#endif
}

template <class T, class BaseM, class M, class V1, class V2> inline void DoTestOProda(
    BaseM& a0, const M& a, const V1& v1, const V2& v2, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start OProd"<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"v1 = "<<tmv::Type(v1)<<" "<<v1<<std::endl;
    std::cout<<"v2 = "<<tmv::Type(v2)<<" "<<v2<<std::endl;
  }

  a0 = a;
  tmv::Matrix<T> m1 = a;
  tmv::Matrix<T> m2 = v1^v2;

  double eps1 = EPS*Norm(m1);
  double eps2 = EPS*Norm(v1)*Norm(v2);

  a += v1^v2;
  m1 += m2;
  Assert(Norm(a-m1) <= eps1+eps2,label+" a += v1^v2");
  m1 = a = a0;
  a += T(7) * (v1^v2);
  m1 += T(7) * m2;
  Assert(Norm(a-m1) <= eps1+7.*eps2,label+" a += 7 * (v1^v2)");
  m1 = a = a0;
#ifdef XTEST
  a += (T(7) * v1)^v2;
  m1 += T(7) * m2;
  Assert(Norm(a-m1) <= eps1+7.*eps2,label+" a += (7*v1) ^ v2)");
  m1 = a = a0;
  a += v1 ^ (T(7) * v2);
  m1 += T(7) * m2;
  Assert(Norm(a-m1) <= eps1+7.*eps2,label+" a += v1 ^ (7*v2)");
  m1 = a = a0;
  a += -v1^v2;
  m1 -= m2;
  Assert(Norm(a-m1) <= eps1+eps2,label+" a += -v1^v2");
  m1 = a = a0;
  a -= v1^v2;
  m1 -= m2;
  Assert(Norm(a-m1) <= eps1+eps2,label+" a -= v1^v2");
  m1 = a = a0;
  a -= -v1^v2;
  m1 -= -m2;
  Assert(Norm(a-m1) <= eps1+eps2,label+" a -= -v1^v2");
  m1 = a = a0;
  a -= T(7) * (v1^v2);
  m1 -= T(7) * m2;
  Assert(Norm(a-m1) <= eps1+7.*eps2,label+" a -= 7 * (v1^v2)");
  m1 = a = a0;
  a -= (T(7) * v1)^v2;
  m1 -= T(7) * m2;
  Assert(Norm(a-m1) <= eps1+7.*eps2,label+" a -= (7*v1) ^ v2");
  m1 = a = a0;
  a -= v1 ^ (T(7) * v2);
  m1 -= T(7) * m2;
  Assert(Norm(a-m1) <= eps1+7.*eps2,label+" a -= v1 ^ (7*v2)");
  m1 = a = a0;
#endif

  if (showstartdone)
    std::cout<<"Done OProd"<<std::endl;
}

template <class T, class BaseM, class M, class V1, class V2> inline void DoTestOProdR(
    BaseM& a0, const M& a, const V1& v1, const V2& v2, std::string label)
{
  DoTestOProda<T>(a0,a,v1,v2,label);
  DoTestOProda<T>(a0,a,v1.Reverse(),v2.Reverse(),label+" RevBC");
#ifndef SYMOPROD
  DoTestOProda<T>(a0,a,v1.Reverse(),v2,label+" RevB");
  DoTestOProda<T>(a0,a,v1,v2.Reverse(),label+" RevC");
#endif
}

template <class T, class BaseM, class M, class V1, class V2> inline void DoTestOProdC(
    BaseM& a0, const M& a, const V1& v1, const V2& v2, std::string label)
{
  DoTestOProda<T>(a0,a,v1,v2,label);
  DoTestOProda<T>(a0,a,v1.Reverse(),v2.Reverse(),label+" RevBC");
  DoTestOProda<T>(a0,Conjugate(a),v1,v2,label+" ConjA");
  DoTestOProda<T>(a0,Conjugate(a),v1.Reverse(),v2.Reverse(),
      label+" ConjA RevBC");
#ifndef SYMOPROD
  DoTestOProda<T>(a0,a,v1.Reverse(),v2,label+" RevB");
  DoTestOProda<T>(a0,a,v1,v2.Reverse(),label+" RevC");
  DoTestOProda<T>(a0,Conjugate(a),v1.Reverse(),v2,label+" ConjA RevB");
  DoTestOProda<T>(a0,Conjugate(a),v1,v2.Reverse(),label+" ConjA RevC");
#endif
}

template <class T, class BaseM, class BaseCM, class M, class CM> 
inline void TestMatrixArith1(
    BaseM& a0, BaseCM& ca0, const M& a, const CM& ca,
    std::string label)
{
  if (showstartdone) {
    std::cout<<"Start TestMatrixArith1 "<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"a0 = "<<tmv::Type(a0)<<" "<<a0<<std::endl;
    std::cout<<"ca = "<<tmv::Type(ca)<<" "<<ca<<std::endl;
    std::cout<<"ca0 = "<<tmv::Type(ca0)<<" "<<ca0<<std::endl;
  }

  CT z(9,-2);
#ifdef XTEST
  T x = 12;
  DoTestMR<T>(a,label+" R");
#endif
  DoTestMC<CT>(ca,label+" C");

#ifdef XTEST
  DoTestMX1R<T>(a,x,label+" R,R");
  DoTestMX2R<T>(a0,a,x,label+" R,R");
  DoTestMX1R<T>(a,z,label+" R,C");
  DoTestMX1C<CT>(ca,x,label+" C,R");
  DoTestMX2C<CT>(ca0,ca,x,label+" C,R");
#endif
  DoTestMX1C<CT>(ca,z,label+" C,C");
  DoTestMX2C<CT>(ca0,ca,z,label+" C,C");

#ifdef ALIASOK
#ifdef XTEST
  DoTestMM2RR<T>(a0,a,a,label+" self_arith");
  DoTestMM4RR<T>(a0,a,a,label+" self_arith");
#endif
  DoTestMM2CC<T>(ca0,ca,ca,label+" self_arith");
  DoTestMM4CC<T>(ca0,ca,ca,label+" self_arith");
#endif

  if (showstartdone)
    std::cout<<"Done Test1"<<std::endl;
}

template <class T, class M, class CM, class V, class CV> 
inline void TestMatrixArith2(
    const M& a, const CM& ca, const V& b, const CV& cb, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start TestMatrixArith2 "<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"ca = "<<tmv::Type(ca)<<" "<<ca<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" "<<b<<std::endl;
    std::cout<<"cb = "<<tmv::Type(cb)<<" "<<cb<<std::endl;
  }

#ifdef XTEST
  DoTestMV1R<T,T>(a,b,label+" R,R");
  DoTestMV2R<T,T>(a,b,label+" R,R");
  DoTestMV1R<T,CT>(a,cb,label+" R,C");
  DoTestMV2R<T,CT>(a,cb,label+" R,C");
  DoTestMV1C<CT,T>(ca,b,label+" C,R");
#endif
  DoTestMV1C<CT,CT>(ca,cb,label+" C,C");
  DoTestMV2C<CT,CT>(ca,cb,label+" C,C");

  if (showstartdone)
    std::cout<<"Done Test2"<<std::endl;
}

template <class T, class M, class CM, class V1, class CV1, class V2, class CV2> 
inline void TestMatrixArith3(
    const M& a, const CM& ca, const V1& b, const CV1& cb,
    const V2& c, const CV2& cc, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start TestMatrixArith3 "<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"ca = "<<tmv::Type(ca)<<" "<<ca<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" "<<b<<std::endl;
    std::cout<<"cb = "<<tmv::Type(cb)<<" "<<cb<<std::endl;
    std::cout<<"c = "<<tmv::Type(c)<<" "<<c<<std::endl;
    std::cout<<"cc = "<<tmv::Type(cc)<<" "<<cc<<std::endl;
  }

#ifdef XTEST
  DoTestMV3R<T,T,T>(a,b,c,label+" R,R,R");
  DoTestMV3R<T,T,CT>(a,b,cc,label+" C,R,R");
  DoTestMV3R<T,CT,CT>(a,cb,cc,label+" C,R,C");
  DoTestMV3C<CT,T,CT>(ca,b,cc,label+" C,C,R");
#endif
  DoTestMV3C<CT,CT,CT>(
      ca,cb,cc,label+" C,C,C");

  if (showstartdone)
    std::cout<<"Done Test2"<<std::endl;
}

template <class T, class BaseM, class BaseCM, class M1, class CM1> 
inline void TestMatrixArith123(
    BaseM& a0, BaseCM& ca0, const M1& a, const CM1& ca, 
    std::string label)
{
  TestMatrixArith1<T>(a0,ca0,a,ca,label);

  tmv::Vector<T> v(a.rowsize());
  for(size_t i=0;i<a.rowsize();i++) v(i) = T(i+3);
  tmv::Vector<CT> cv = CT(4,5) * v;
  tmv::VectorView<T> vv = v.View();
  tmv::VectorView<CT> cvv = cv.View();

  tmv::Vector<T> v5(5*a.rowsize());
  tmv::Vector<CT> cv5(5*a.rowsize());
  tmv::VectorView<T> vs = v5.SubVector(0,5*a.rowsize(),5);
  tmv::VectorView<CT> cvs = cv5.SubVector(0,5*a.rowsize(),5);
  vs = vv;
  cvs = cvv;

  tmv::Vector<T> w(a.colsize());
  for(size_t i=0;i<a.colsize();i++) w(i) = T(2*i-6);
  tmv::Vector<CT> cw = CT(-1,2) * w;
  tmv::VectorView<T> wv = w.View();
  tmv::VectorView<CT> cwv = cw.View();

  tmv::Vector<T> w5(5*a.colsize());
  tmv::Vector<CT> cw5(5*a.colsize());
  tmv::VectorView<T> ws = w5.SubVector(0,5*a.colsize(),5);
  tmv::VectorView<CT> cws = cw5.SubVector(0,5*a.colsize(),5);
  ws = wv;
  cws = cwv;

  TestMatrixArith2<T>(a,ca,vv,cvv,label);
  TestMatrixArith2<T>(a,ca,vs,cvs,label);

  TestMatrixArith3<T>(a,ca,vv,cvv,wv,cwv,label);
  TestMatrixArith3<T>(a,ca,vv,cvv,ws,cws,label);
  TestMatrixArith3<T>(a,ca,vs,cvs,wv,cwv,label);
  TestMatrixArith3<T>(a,ca,vs,cvs,ws,cws,label); 
}

template <class T, class BaseM, class BaseCM, class M1, class CM1, class M2, class CM2> 
inline void TestMatrixArith4(
    BaseM& a0, BaseCM& ca0,
    const M1& a, const CM1& ca, const M2& b, const CM2& cb, 
    std::string label)
{
  if (showstartdone) {
    std::cout<<"Start TestMatrixArith4 "<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" "<<b<<std::endl;
    std::cout<<"ca = "<<tmv::Type(ca)<<" "<<ca<<std::endl;
    std::cout<<"cb = "<<tmv::Type(cb)<<" "<<cb<<std::endl;
    std::cout<<"a0 = "<<tmv::Type(a0)<<std::endl;
    std::cout<<"ca0 = "<<tmv::Type(ca0)<<std::endl;
  }

#ifdef XTEST
  DoTestMM1RR<T>(a,b,label+" R,R");
  DoTestMM2RR<T>(a0,a,b,label+" R,R");
  DoTestMM1RC<T>(a,cb,label+" R,C");
  DoTestMM1CR<T>(ca,b,label+" C,R");
  DoTestMM2CR<T>(ca0,ca,b,label+" C,R");
#endif
  DoTestMM1CC<T>(ca,cb,label+" C,C");
  DoTestMM2CC<T>(ca0,ca,cb,label+" C,C");

  if (showstartdone)
    std::cout<<"Done Test4"<<std::endl;
}

template <class T, class BaseM, class BaseCM, class M1, class CM1, class M2, class CM2> 
inline void TestMatrixArith5(
    BaseM& a0, BaseCM& ca0,
    const M1& a, const CM1& ca, const M2& b, const CM2& cb, 
    std::string label)
{
  if (showstartdone) {
    std::cout<<"Start TestMatrixArith5 "<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<" "<<b<<std::endl;
    std::cout<<"ca = "<<tmv::Type(ca)<<" "<<ca<<std::endl;
    std::cout<<"cb = "<<tmv::Type(cb)<<" "<<cb<<std::endl;
    std::cout<<"a0 = "<<tmv::Type(a0)<<std::endl;
    std::cout<<"ca0 = "<<tmv::Type(ca0)<<std::endl;
  }

#ifdef XTEST
  DoTestMM3RR<T>(a,b,label+" R,R");
  DoTestMM4RR<T>(a0,a,b,label+" R,R");
  DoTestMM3RC<T>(a,cb,label+" R,C");
  DoTestMM3CR<T>(ca,b,label+" C,R");
  DoTestMM4CR<T>(ca0,ca,b,label+" C,R");
#endif
  DoTestMM3CC<T>(ca,cb,label+" C,C");
  DoTestMM4CC<T>(ca0,ca,cb,label+" C,C");

  if (showstartdone)
    std::cout<<"Done Test5"<<std::endl;
}

template <class T, class BaseM, class BaseCM, class M1, class CM1, class M2, class CM2> 
inline void TestMatrixArith45(
    BaseM& a0, BaseCM& ca0,
    const M1& a, const CM1& ca, const M2& b, const CM2& cb, 
    std::string label)
{
  TestMatrixArith4<T>(a0,ca0,a,ca,b,cb,label);
  TestMatrixArith5<T>(a0,ca0,a,ca,b,cb,label);
}

template <class T, class BaseM, class BaseCM, class M, class CM, class V1, class CV1, class V2, class CV2> 
inline void TestMatrixArith6(
    BaseM& a0, BaseCM& ca0, const M& a, const CM& ca,
    const V1& v1, const CV1& cv1, const V2& v2, const CV2& cv2, 
    std::string label)
{
  if (showstartdone) {
    std::cout<<"Start TestMatrixArith6 "<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<" "<<a<<std::endl;
    std::cout<<"ca = "<<tmv::Type(ca)<<" "<<ca<<std::endl;
    std::cout<<"v1 = "<<tmv::Type(v1)<<" "<<v1<<std::endl;
    std::cout<<"cv1 = "<<tmv::Type(cv1)<<" "<<cv1<<std::endl;
    std::cout<<"v2 = "<<tmv::Type(v2)<<" "<<v2<<std::endl;
    std::cout<<"cv2 = "<<tmv::Type(cv2)<<" "<<cv2<<std::endl;
    std::cout<<"a0 = "<<tmv::Type(a0)<<std::endl;
    std::cout<<"ca0 = "<<tmv::Type(ca0)<<std::endl;
  }

#ifdef XTEST
  DoTestOProdR<T>(a0,a,v1,v2,label+" R,R,R");
  DoTestOProdC<CT>(ca0,ca,v1,v2,label+" C,R,R");
#ifndef SYMOPROD
  DoTestOProdC<CT>(ca0,ca,cv1,v2,label+" C,C,R");
  DoTestOProdC<CT>(ca0,ca,v1,cv2,label+" C,R,C");
#endif
#endif
  DoTestOProdC<CT>(ca0,ca,cv1,cv2,label+" C,C,C");

  if (showstartdone)
    std::cout<<"Done Test6"<<std::endl;
}

#undef CT
