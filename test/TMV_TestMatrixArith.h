// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
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

template <class M1, class M2> inline bool CanMultMM(const M1& a, const M2& b)
{ return a.rowsize() == b.colsize(); }

template <class M, class V> inline bool CanMultMV(const M& m, const V& v)
{ return m.rowsize() == v.size(); }

template <class M, class V> inline bool CanMultVM(const V& v, const M& m)
{ return m.colsize() == v.size(); }

template <class M1, class M2, class M3> inline bool CanMultMM(
    const M1& a, const M2& b, const M3& c)
{
  return CanMultMM(a,b) && c.colsize() == a.colsize() && 
  c.rowsize() == b.rowsize(); 
}

template <class M1, class V2, class V3> inline bool CanMultMV(
    const M1& a, const V2& b, const V3& c)
{ return CanMultMV(a,b) && c.size() == a.colsize(); }

template <class V1, class M2, class V3> inline bool CanMultVM(
    const V1& a, const M2& b, const V3& c)
{ return CanMultVM(a,b) && c.size() == b.rowsize(); }

template <class M1, class M2, class M3> inline bool CanMultXMM(const M1& a, const M2& b, const M3& c)
{ return CanMultMM(a,b,c); }

#ifndef TMV_TESTVECTOR_H
template <class T1, class T2> struct ProdType
{ typedef T1 Tprod; };

template <class T> struct ProdType<T,std::complex<T> >
{ typedef std::complex<T> Tprod; };

#define ProductType(T1,T2) typename ProdType<T1,T2>::Tprod
#endif

#ifdef NOVIEWS
#define CONST
#endif

#ifndef CONST
#define CONST const
#endif

#ifdef USETEMP
#define IFTEMP(x) x
#define IFTEMP1(x) x,
#else
#define IFTEMP(x) 
#define IFTEMP1(x) 
#endif

#ifdef NOMIX
#define VEC(T,v) (tmv::Vector<T>(v))
#define VEC2(T1,T2,v) (tmv::Vector<ProductType(T1,T2)>(v))
#define MAT(T,m) (tmv::Matrix<T>(m))
#define MAT2(T1,T2,m) (tmv::Matrix<ProductType(T1,T2)>(m))
#else
#define VEC(T,v) (v)
#define VEC2(T1,T2,v) (v)
#define MAT(T,m) (m)
#define MAT2(T1,T2,m) (m)
#endif

template <class T, class MM> static void DoTestMa(const MM& a, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start Ma "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<std::endl;
  }

  tmv::Matrix<T> m = a;
  double eps = EPS * Norm(m);

  if (XXDEBUG1) {
    std::cout<<"a = "<<tmv::TypeText(a)<<" = "<<a<<std::endl;
    std::cout<<"m = "<<tmv::TypeText(m)<<" = "<<m<<std::endl;
    std::cout<<"a-m = "<<a-m<<std::endl;
    std::cout<<"Norm(a-m) = "<<Norm(MAT(T,a)-m)<<std::endl;
#ifndef NONSQUARE
    if (m.IsSquare()) {
      std::cout<<"Trace(a) = "<<Trace(a)<<"  "<<Trace(m)<<std::endl;
      std::cout<<"Det(a) = "<<Det(a)<<"  "<<Det(m)<<std::endl;
      std::cout<<"diff = "<<std::abs(Det(a)-Det(m))<<"  "<<eps*std::abs(Det(m)+Norm(m.Inverse()))<<std::endl;
      std::cout<<"LogDet(a) = "<<LogDet(a)<<"  "<<LogDet(m)<<std::endl;
      std::cout<<"diff = "<<std::abs(LogDet(a)-LogDet(m))<<"  "<<m.colsize()*eps*Norm(m.Inverse())<<std::endl;
    }
#endif
    std::cout<<"NormF(a) = "<<NormF(a)<<"  "<<NormF(m)<<std::endl;
    std::cout<<"Norm(a) = "<<Norm(a)<<"  "<<Norm(m)<<std::endl;
    std::cout<<"Norm1(a) = "<<Norm1(a)<<"  "<<Norm1(m)<<std::endl;
    std::cout<<"Norm1(m) = "<<Norm1(m)<<std::endl;
    std::cout<<"NormInf(a) = "<<NormInf(a)<<"  "<<NormInf(m)<<std::endl;
    std::cout<<"abs(diff) = "<<std::abs(NormInf(a)-NormInf(m))<<std::endl;
    std::cout<<"eps*norminf = "<<EPS*NormInf(m)<<std::endl;
  }

  Assert(Norm(MAT(T,a)-m) <= eps,label+" a != m");
#ifndef NONSQUARE
  if (m.IsSquare()) {
    Assert(std::abs(Trace(a)-Trace(m)) <= eps,label+" Trace");
    T d = Det(m);
    if (std::abs(d) > 0.5) {
      double eps1 = eps * Norm(m.Inverse());
      Assert(std::abs(Det(a)-d) <= eps1*std::abs(d),label+" Det");
      Assert(std::abs(LogDet(a)-LogDet(m)) <= m.colsize()*eps1,label+" LogDet");
    } else if (std::abs(d) != 0.0) {
      double eps1 = eps * Norm(m.Inverse());
      Assert(std::abs(Det(a)-d) <= eps1*(1.+std::abs(d)),label+" Det");
    } else {
      Assert(std::abs(Det(a)) <= eps,label+" Det");
    }
  }
#endif
  Assert(std::abs(NormF(a)-NormF(m)) <= eps,label+" NormF");
  Assert(std::abs(Norm(a)-Norm(m)) <= eps,label+" Norm");
  Assert(std::abs(Norm1(a)-Norm1(m)) <= eps,label+" Norm1");
  Assert(std::abs(NormInf(a)-NormInf(m)) <= eps,label+" NormInf");
#ifdef XTEST
  Assert(std::abs(a.NormF()-NormF(m)) <= eps,label+" NormF");
  Assert(std::abs(a.Norm()-Norm(m)) <= eps,label+" Norm");
  Assert(std::abs(a.Norm1()-Norm1(m)) <= eps,label+" Norm1");
  Assert(std::abs(a.NormInf()-NormInf(m)) <= eps,label+" NormInf");
#endif
  if (donorm2) {
    if (XXDEBUG1) {
      std::cout<<"Norm2(a) = "<<a.DoNorm2()<<"  "<<m.DoNorm2()<<std::endl;
      std::cout<<"abs(diff) = "<<std::abs(a.DoNorm2()-m.DoNorm2())<<std::endl;
      std::cout<<"eps*kappa = "<<eps*m.DoCondition()<<std::endl;
    }
    Assert(std::abs(a.DoNorm2()-m.DoNorm2()) <= eps*a.DoCondition(),
        label+" DoNorm2");
#ifndef NOSV
    a.DivideUsing(tmv::SV);
    a.SetDiv();
    m.DivideUsing(tmv::SV);
    m.SetDiv();
    Assert(std::abs(Norm2(a)-Norm2(m)) <= eps*a.Condition(),label+" Norm2");
#ifdef XTEST
    Assert(std::abs(a.Norm2()-m.DoNorm2()) <= eps*a.Condition(),
        label+" Norm2");
#endif
#endif
  }
  if (XXDEBUG1) {
    std::cout<<"Norm(aT-mT) = "<<Norm(Transpose(a)-Transpose(m))<<std::endl;
    std::cout<<"Conjugate(a) = "<<Conjugate(a)<<std::endl;
    std::cout<<"Conjugate(m) = "<<Conjugate(m)<<std::endl;
    std::cout<<"a*-m* = "<<Conjugate(a)-Conjugate(m)<<std::endl;
    std::cout<<"Conjugate(a).diag = "<<Conjugate(a).diag()<<std::endl;
    std::cout<<"Conjugate(m).diag = "<<Conjugate(m).diag()<<std::endl;
    std::cout<<"Norm(a*-m*) = "<<Norm(MAT(T,Conjugate(a))-Conjugate(m))<<std::endl;
    std::cout<<"Norm(at-mt) = "<<Norm(MAT(T,Adjoint(a))-Adjoint(m))<<std::endl;
  }
  Assert(Norm(MAT(T,Transpose(a))-Transpose(m)) <= eps,label+" Transpose");
  Assert(Norm(MAT(T,Conjugate(a))-Conjugate(m)) <= eps,label+" Conjugate");
  Assert(Norm(MAT(T,Adjoint(a))-Adjoint(m)) <= eps,label+" Adjoint");

  if (showstartdone)
    std::cout<<"Done Ma"<<std::endl;
}

template <class T, class MM> static void DoTestMR(
    const MM& a, std::string label)
{
  DoTestMa<T>(a,label);
#ifndef NOVIEWS
  DoTestMa<T>(Transpose(a),label+" Trans");
#endif
}

template <class T, class MM> static void DoTestMC(
    const MM& a, std::string label)
{
  DoTestMa<T>(a,label);
#ifndef NOVIEWS
  DoTestMa<T>(Transpose(a),label+" Trans");
  DoTestMa<T>(Conjugate(a),label+" Conj");
  DoTestMa<T>(Adjoint(a),label+" Adj");
#endif
}

template <class T, IFTEMP1(class M0) class MM, class T2> static void DoTestMX1a(
    IFTEMP1(M0& temp) const MM& a, T2 x, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MX1a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"x = "<<tmv::TypeText(x)<<" "<<x<<std::endl;
  }

  tmv::Matrix<T> m = a;

  double eps = EPS*Norm(m);
#ifndef NONSQUARE
  if (CanAddX(a,x)) {
    if (XXDEBUG2) {
      std::cout<<"CanAddX("<<tmv::TypeText(a)<<","<<tmv::TypeText(x)<<")\n";
      std::cout<<"x = "<<tmv::TypeText(x)<<"  "<<x<<std::endl;
      std::cout<<"a = "<<tmv::TypeText(a)<<"  "<<a<<std::endl;
      std::cout<<"x-a = "<<(IFTEMP(temp=)x-a)<<std::endl;
      std::cout<<"x-m = "<<x-m<<std::endl;
      std::cout<<"(x-a)-(x-m) = "<<MAT2(T,T2,IFTEMP(temp=)x-a)-(x-m)<<std::endl;
    }
    Assert(Norm(MAT(T,a)-m) <= eps,label+" a != m1");
    Assert(Norm(MAT2(T,T2,IFTEMP(temp=)x-a)-(x-m)) <= eps,label+" x-a");
#ifdef XTEST
    Assert(Norm(MAT2(T,T2,IFTEMP(temp=)a-x)-(m-x)) <= eps,label+" a-x");
    Assert(Norm(MAT2(T,T2,IFTEMP(temp=)x+a)-(x+m)) <= eps,label+" x+a");
    Assert(Norm(MAT2(T,T2,IFTEMP(temp=)a+x)-(m+x)) <= eps,label+" a+x");
#endif
  }
#endif
  if (CanMultX(a,x)) {
    if (XXDEBUG2) {
      std::cout<<"CanMultX("<<tmv::TypeText(a)<<","<<tmv::TypeText(x)<<")\n";
    }
    Assert(Norm(MAT2(T,T2,IFTEMP(temp=)x*a)-(x*m)) <= eps*std::abs(x),label+" x*a");
#ifdef XTEST
    Assert(Norm(MAT2(T,T2,IFTEMP(temp=)a*x)-(x*m)) <= eps*std::abs(x),label+" a*x");
    if (tmv::Epsilon<T>()  != T(0)) {
      Assert(Norm(MAT2(T,T2,IFTEMP(temp=)a/x)-(m/x)) <= eps/std::abs(x),label+" a/x");
    }
#endif
  }
  if (showstartdone)
    std::cout<<"Done MX1a"<<std::endl;
}

template <class T, IFTEMP1(class M0) class MM, class T2> static void DoTestMX1R(
    IFTEMP1(M0& temp) const MM& a, T2 x, std::string label)
{
  DoTestMX1a<T>(IFTEMP1(temp) a,x,label);
}

template <class T, IFTEMP1(class M0) class MM, class T2> static void DoTestMX1C(
    IFTEMP1(M0& temp) const MM& a, T2 x, std::string label)
{
  DoTestMX1a<T>(IFTEMP1(temp) a,x,label);
#ifndef NOVIEWS
  DoTestMX1a<T>(IFTEMP1(temp) Conjugate(a),x,label+" Conj");
#endif
}

template <class T, class BaseM, class MM, class T2> static void DoTestMX2a(
    BaseM& a0, CONST MM& a, T2 x, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MX2a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"a0 = "<<tmv::TypeText(a0)<<" "<<a0<<std::endl;
    std::cout<<"x = "<<tmv::TypeText(x)<<" "<<x<<std::endl;
  }

  a0 = a;
  tmv::Matrix<T> m1 = a;
  tmv::Matrix<T> temp = a;
  Assert(Norm(MAT(T,a)-m1) <= EPS*Norm(m1),label+" a = m1");

  double eps = EPS * Norm(m1);
#ifndef NONSQUARE
  if (CanAddEqX(a,x)) {
    if (XXDEBUG3) {
      std::cout<<"CanAddEqX("<<tmv::TypeText(a)<<","<<tmv::TypeText(x)<<")\n";
    }
    tmv::Matrix<T> m2 = m1;
    a += x;
    m2 = m1+x;
    Assert(Norm(MAT(T,a)-m2) <= eps,label+" a += x");
    Assert(Norm(MAT(T,a+=x)-(m2+=x)) <= eps,label+" a += x (2)");
    a = a0;
    a = a+x; 
    m2 = m1+x;
    Assert(Norm(MAT(T,a)-m2) <= eps,label+" a = a+x");
    a = a0;
#ifdef XTEST
    a += -x;
    m2 = m1-x;
    Assert(Norm(MAT(T,a)-m2) <= eps,label+" a += x");
    Assert(Norm(MAT(T,a+=-x)-(m2+=-x)) <= eps,label+" a += -x (2)");
    a = a0;
    a -= x;
    m2 = m1-x;
    Assert(Norm(MAT(T,a)-m2) <= eps,label+" a -= x");
    Assert(Norm(MAT(T,a-=x)-(m2-=x)) <= eps,label+" a -= x (2)");
    a = a0;
    a -= -x;
    m2 = m1+x;
    Assert(Norm(MAT(T,a)-m2) <= eps,label+" a -= x");
    Assert(Norm(MAT(T,a-=-x)-(m2-=-x)) <= eps,label+" a -= -x (2)");
    a = a0;
    a = a-x;
    m2 = m1-x;
    Assert(Norm(MAT(T,a)-m2) <= eps,label+" a = a-x");
    a = a0;
    a = x+a;
    m2 = x+m1;
    Assert(Norm(MAT(T,a)-m2) <= eps,label+" a = x+a");
    a = a0;
    a = x-a; 
    m2 = x-m1;
    Assert(Norm(MAT(T,a)-m2) <= eps,label+" a = x-a");
    a = a0;
#endif
  }
#endif
  if (CanMultEqX(a,x)) {
    if (XXDEBUG3) {
      std::cout<<"CanMultEqX("<<tmv::TypeText(a)<<","<<tmv::TypeText(x)<<")\n";
    }
    tmv::Matrix<T> m2 = m1;
    a *= x;
    m2 = m1*x;
    Assert(Norm(MAT(T,a)-m2) <= eps*std::abs(x),label+" a *= x");
    Assert(Norm(MAT(T,a*=x)-(m2*=x)) <= eps*std::abs(x*x),label+" a *= x");
    a = a0;
    a = a*x;
    m2 = m1*x;
    Assert(Norm(MAT(T,a)-m2) <= eps*std::abs(x),label+" a = a*x");
    a = a0;
#ifdef XTEST
    a *= -x;
    m2 = -m1*x;
    Assert(Norm(MAT(T,a)-m2) <= eps*std::abs(x),label+" a *= -x");
    Assert(Norm(MAT(T,a*=-x)-(m2*=-x)) <= eps*std::abs(x*x),label+" a *= -x");
    a = a0;
    if (tmv::Epsilon<T>() != T(0)) {
      a /= x;
      m2 = m1/x;
      Assert(Norm(MAT(T,a)-m2) <= eps*std::abs(x),label+" a /= x");
      Assert(Norm(MAT(T,a/=x)-(m2/=x)) <= eps,label+" a /= x");
      a = a0;
      a /= -x;
      m2 = -m1/x;
      Assert(Norm(MAT(T,a)-m2) <= eps*std::abs(x),label+" a /= -x");
      Assert(Norm(MAT(T,a/=-x)-(m2/=-x)) <= eps,label+" a /= -x");
      a = a0;
      a = a/x;
      m2 = m1/x;
      Assert(Norm(MAT(T,a)-m2) <= eps,label+" a = a/x");
      a = a0;
    }
    a = x*a; 
    m2 = x*m1;
    Assert(Norm(MAT(T,a)-m2) <= eps*std::abs(x),label+" a = x*a");
    a = a0;
#endif
  }

  if (showstartdone)
    std::cout<<"Done MX2a"<<std::endl;
}

template <class T, class BaseM, class MM, class T2> static void DoTestMX2R(
    BaseM& a0, CONST MM& a, T2 x, std::string label)
{
  DoTestMX2a<T>(a0,a,x,label);
}

template <class T, class BaseM, class MM, class T2> static void DoTestMX2C(
    BaseM& a0, CONST MM& a, T2 x, std::string label)
{
  DoTestMX2a<T>(a0,a,x,label);
#ifndef NOVIEWS
  DoTestMX2a<T>(a0,Conjugate(a),x,label+" Conj");
#endif
}

// m*v, v*m
template <class Ta, class T, IFTEMP1(class V0) IFTEMP1(class CV0) class MM, class V> static void DoTestMV1a(
    IFTEMP1(V0& temp) IFTEMP1(CV0& ctemp) const MM& a, const V& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV1a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
  }

  tmv::Matrix<Ta> m = a;
  tmv::Vector<T> v = b;
  double eps = EPS * Norm(m) * Norm(v);
  if (CanMultMV(a,b)) {
    if (XXDEBUG4) {
      std::cout<<"CanMult("<<tmv::TypeText(m)<<","<<tmv::TypeText(b)<<")\n";
      std::cout<<"a*b = "<<(IFTEMP(temp=)a*b)<<std::endl;
      std::cout<<"m*v = "<<m*v<<std::endl;
      std::cout<<"a*b-m*v = "<<(VEC2(T,Ta,IFTEMP(temp=)a*b)-(m*v))<<std::endl;
      std::cout<<"Norm(a*b-m*v) = "<<Norm(VEC2(T,Ta,IFTEMP(temp=)a*b)-(m*v))<<std::endl;
      std::cout<<"eps = "<<eps<<std::endl;
    }
    Assert(Norm(VEC2(T,Ta,IFTEMP(temp=)a*b)-(m*v)) <= eps,label+" a*b");
    RealType(T) x(5);
    ComplexType(T) z(3,4);
    Assert(Norm(VEC2(T,Ta,IFTEMP(temp=)x*a*b)-(x*m*v)) <= x*eps,label+" x*a*b");
    Assert(Norm(VEC(ComplexType(T),IFTEMP(ctemp=)z*a*b)-(z*m*v)) <= x*eps,label+" z*a*b");
#ifdef XTEST
    Assert(Norm(VEC2(T,Ta,IFTEMP(temp=)(x*a)*b)-(x*m*v)) <= x*eps,label+" (x*a)*b");
    Assert(Norm(VEC2(T,Ta,IFTEMP(temp=)x*(a*b))-(x*m*v)) <= x*eps,label+" x*(a*b)");
    Assert(Norm(VEC2(T,Ta,IFTEMP(temp=)a*(x*b))-(x*m*v)) <= x*eps,label+" a*(x*b)");
    Assert(Norm(VEC(ComplexType(T),IFTEMP(ctemp=)(z*a)*b)-(z*m*v)) <= x*eps,label+" (z*a)*b");
    Assert(Norm(VEC(ComplexType(T),IFTEMP(ctemp=)z*(a*b))-(z*m*v)) <= x*eps,label+" z*(a*b)");
    Assert(Norm(VEC(ComplexType(T),IFTEMP(ctemp=)a*(z*b))-(z*m*v)) <= x*eps,label+" a*(z*b)");
#endif
  }
  if (showstartdone)
    std::cout<<"Done MV1a"<<std::endl;
}

template <class Ta, class T, IFTEMP1(class V0) IFTEMP1(class CV0) class MM, class V> static void DoTestVM1a(
    IFTEMP1(V0& temp) IFTEMP1(CV0& ctemp) const MM& a, const V& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start VM1a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
  }

  tmv::Matrix<Ta> m = a;
  tmv::Vector<T> v = b;
  double eps = EPS * Norm(m) * Norm(v);
  if (CanMultVM(v,m)) {
    RealType(T) x(5);
    ComplexType(T) z(3,4);
    if (XXDEBUG4) {
      std::cout<<"CanMult("<<tmv::TypeText(v)<<","<<tmv::TypeText(m)<<")\n";
      std::cout<<"v*m = "<<v*m<<std::endl;
      std::cout<<"b*a = "<<(IFTEMP(temp=)b*a)<<std::endl;
      std::cout<<"x*b*a = "<<(IFTEMP(temp=)x*b*a)<<std::endl;
      std::cout<<"z*b*a = "<<(IFTEMP(ctemp=)z*b*a)<<std::endl;
      std::cout<<"b*a Norm(diff) = "<<Norm(VEC2(T,Ta,IFTEMP(temp=)b*a)-(v*m))<<std::endl;
      std::cout<<"x*b*a Norm(diff) = "<<Norm(VEC2(T,Ta,IFTEMP(temp=)x*b*a)-(x*v*m))<<std::endl;
      std::cout<<"z*b*a Norm(diff) = "<<Norm(VEC(ComplexType(T),IFTEMP(ctemp=)z*b*a)-(z*v*m))<<std::endl;
      std::cout<<"eps = "<<eps<<std::endl;
      std::cout<<"x*eps = "<<x*eps<<std::endl;
    }
    Assert(Norm(VEC2(T,Ta,IFTEMP(temp=)b*a)-(v*m)) <= eps,label+" b*a");
    Assert(Norm(VEC2(T,Ta,IFTEMP(temp=)x*b*a)-(x*v*m)) <= x*eps,label+" x*b*a");
    Assert(Norm(VEC(ComplexType(T),IFTEMP(ctemp=)z*b*a)-(z*v*m)) <= x*eps,label+" z*b*a");
#ifdef XTEST
    Assert(Norm(VEC2(T,Ta,IFTEMP(temp=)(x*b)*a)-(x*v*m)) <= x*eps,label+" (x*b)*a");
    Assert(Norm(VEC2(T,Ta,IFTEMP(temp=)x*(b*a))-(x*v*m)) <= x*eps,label+" x*(b*a)");
    Assert(Norm(VEC2(T,Ta,IFTEMP(temp=)b*(x*a))-(x*v*m)) <= x*eps,label+" b*(x*a)");
    Assert(Norm(VEC(ComplexType(T),IFTEMP(ctemp=)(z*b)*a)-(z*v*m)) <= x*eps,label+" (z*b)*a");
    Assert(Norm(VEC(ComplexType(T),IFTEMP(ctemp=)z*(b*a))-(z*v*m)) <= x*eps,label+" z*(b*a)");
    Assert(Norm(VEC(ComplexType(T),IFTEMP(ctemp=)b*(z*a))-(z*v*m)) <= x*eps,label+" b*(z*a)");
#endif
  }
  if (showstartdone)
    std::cout<<"Done VM1a"<<std::endl;
}

template <class Ta, class T, IFTEMP1(class V0) IFTEMP1(class CV0) class MM, class V> static void DoTestMV1R(
    IFTEMP1(V0& temp) IFTEMP1(CV0& ctemp) const MM& a, CONST V& b, std::string label)
{
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
#ifndef NOVIEWS
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b.Reverse(),label);
#endif

#ifdef XTEST
  tmv::Vector<T> b0 = b;

  b.Zero();
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);

  b = b0;
#ifndef NOVIEWS
  b.SubVector(0,b.size()/2).Zero();
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);

  b = b0;
  b.SubVector(b.size()/2,b.size()).Zero();
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);

  b = b0;
  b.SubVector(0,b.size()/4).Zero();
  b.SubVector(3*b.size()/4,b.size()).Zero();
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);

  if (b.size() > 1) {
    b = b0;
    b.SubVector(0,1).Zero();
    DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);

    b = b0;
    b.SubVector(b.size()-1,b.size()).Zero();
    DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);

    b = b0;
    b.SubVector(0,1).Zero();
    b.SubVector(b.size()-1,b.size()).Zero();
    DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
  }
  b = b0; 
#endif
#endif
}

template <class Ta, class T, IFTEMP1(class V0) IFTEMP1(class CV0) class MM, class V> static void DoTestVM1R(
    IFTEMP1(V0& temp) IFTEMP1(CV0& ctemp) const MM& a, CONST V& b, std::string label)
{
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
#ifndef NOVIEWS
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b.Reverse(),label);
#endif

#ifdef XTEST
  tmv::Vector<T> b0 = b;

  b.Zero();
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);

  b = b0;
#ifndef NOVIEWS
  b.SubVector(0,b.size()/2).Zero();
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);

  b = b0;
  b.SubVector(b.size()/2,b.size()).Zero();
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);

  b = b0;
  b.SubVector(0,b.size()/4).Zero();
  b.SubVector(3*b.size()/4,b.size()).Zero();
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);

  if (b.size() > 1) {
    b = b0;
    b.SubVector(0,1).Zero();
    DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);

    b = b0;
    b.SubVector(b.size()-1,b.size()).Zero();
    DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);

    b = b0;
    b.SubVector(0,1).Zero();
    b.SubVector(b.size()-1,b.size()).Zero();
    DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
  }
  b = b0;
#endif
#endif
}

template <class Ta, class T, IFTEMP1(class V0) IFTEMP1(class CV0) class MM, class V> static void DoTestMV1C(
    IFTEMP1(V0& temp) IFTEMP1(CV0& ctemp) const MM& a, CONST V& b, std::string label)
{
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
#ifndef NOVIEWS
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b.Reverse(),label);

  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,label+" Conj");
#endif
#ifdef XTEST
#ifndef NOVIEWS
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b.Reverse(),
      label+" Conj");
#endif

  tmv::Vector<T> b0 = b;

  b.Zero();
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
#ifndef NOVIEWS
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,label+" Conj");
#endif

  b = b0; 
#ifndef NOVIEWS
  b.SubVector(0,b.size()/2).Zero();
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,label+" Conj");

  b = b0; 
  b.SubVector(b.size()/2,b.size()).Zero();
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,label+" Conj");

  b = b0; 
  b.SubVector(0,b.size()/4).Zero();
  b.SubVector(3*b.size()/4,b.size()).Zero();
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
  DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,label+" Conj");

  if (b.size() > 1) {
    b = b0; 
    b.SubVector(0,1).Zero();
    DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
    DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,
        label+" Conj");

    b = b0; 
    b.SubVector(b.size()-1,b.size()).Zero();
    DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
    DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,
        label+" Conj");

    b = b0; 
    b.SubVector(0,1).Zero();
    b.SubVector(b.size()-1,b.size()).Zero();
    DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
    DoTestMV1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,
        label+" Conj");
  }
  b = b0; 
#endif
#endif
}

template <class Ta, class T, IFTEMP1(class V0) IFTEMP1(class CV0) class MM, class V> static void DoTestVM1C(
    IFTEMP1(V0& temp) IFTEMP1(CV0& ctemp) const MM& a, CONST V& b, std::string label)
{
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
#ifndef NOVIEWS
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b.Reverse(),label);

  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,
      label+" Conj");
#endif
#ifdef XTEST
#ifndef NOVIEWS
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b.Reverse(),
      label+" Conj");
#endif

  tmv::Vector<T> b0 = b;

  b.Zero();
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
#ifndef NOVIEWS
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,label+" Conj");
#endif

  b = b0;
#ifndef NOVIEWS
  b.SubVector(0,b.size()/2).Zero();
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,label+" Conj");

  b = b0;
  b.SubVector(b.size()/2,b.size()).Zero();
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,label+" Conj");

  b = b0;
  b.SubVector(0,b.size()/4).Zero();
  b.SubVector(3*b.size()/4,b.size()).Zero();
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
  DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,label+" Conj");

  if (b.size() > 1) {
    b = b0;
    b.SubVector(0,1).Zero();
    DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
    DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,
        label+" Conj");

    b = b0;
    b.SubVector(b.size()-1,b.size()).Zero();
    DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
    DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,
        label+" Conj");

    b = b0;
    b.SubVector(0,1).Zero();
    b.SubVector(b.size()-1,b.size()).Zero();
    DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);
    DoTestVM1a<Ta,T>(IFTEMP1(temp) IFTEMP1(ctemp) Conjugate(a),b,
        label+" Conj");
  }
  b = b0;
#endif
#endif
}

template <class T> inline void SetZ(T& z)
{ z = T(5); }
template <class T> inline void SetZ(std::complex<T>& z)
{ z = std::complex<T>(3,4); }

template <class Ta, class T, class MM, class V> static void DoTestMV2a(
    const MM& a, CONST V& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV2a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
  }

  tmv::Matrix<Ta> m = a;
  tmv::Vector<T> v = b;

  double eps = EPS * Norm(v) * Norm(m);
  double eps2 = EPS * Norm(v) * (1.+Norm(m));

  if (CanMultMV(a,b)) {
    tmv::Vector<T> prod = m*v;
    tmv::Vector<T> c = prod;
    tmv::Vector<T> c0 = c;
    c0 /= T(2);
    if (XXDEBUG5) {
      std::cout<<"CanMult("<<tmv::TypeText(m)<<","<<tmv::TypeText(b)<<","<<tmv::TypeText(v)<<")\n";
      std::cout<<"v = "<<v<<std::endl;
      std::cout<<"m = "<<m<<std::endl;
      std::cout<<"prod = "<<prod<<std::endl;
    }
    c = a*b;
    Assert(Norm(c-prod) <= eps,label+" b=a*v");
    c = c0;
    c += a*b;
    Assert(Norm(c-(c0+prod)) <= eps2,label+" b+=a*v");
    c = c0;
    RealType(T) x(5);
    T z; SetZ(z);
    c = x*a*b;
    Assert(Norm(c-x*prod) <= x*eps,label+" b=x*a*v");
    c = z*a*b;
    Assert(Norm(c-z*prod) <= x*eps,label+" b=z*a*v");
    c = c0;
#ifdef XTEST
    c += x*a*b;
    Assert(Norm(c-(c0+x*prod)) <= x*eps2,label+" b+=x*a*v");
    c = c0;
    c += z*a*b;
    Assert(Norm(c-(c0+z*prod)) <= x*eps2,label+" b+=z*a*v");
    c = c0;
    c = -a*b;
    Assert(Norm(c-(-prod)) <= eps,label+" b=-a*v");
    c = -x*a*b;
    Assert(Norm(c-(-x*prod)) <= x*eps,label+" b=-x*a*v");
    c = -z*a*b;
    Assert(Norm(c-(-z*prod)) <= x*eps,label+" b=-z*a*v");
    c = c0;
    c += -a*b;
    Assert(Norm(c-(c0-prod)) <= eps2,label+" b+=-a*v");
    c = c0;
    c -= a*b;
    Assert(Norm(c-(c0-prod)) <= eps2,label+" b-=a*v");
    c = c0;
    c -= -a*b;
    Assert(Norm(c-(c0+prod)) <= eps2,label+" b-=-a*v");
    c = c0;
    c += -x*a*b;
    Assert(Norm(c-(c0-x*prod)) <= x*eps2,label+" b+=-x*a*v");
    c = c0;
    c -= x*a*b;
    Assert(Norm(c-(c0-x*prod)) <= x*eps2,label+" b-=x*a*v");
    c = c0;
    c -= -x*a*b;
    Assert(Norm(c-(c0+x*prod)) <= x*eps2,label+" b-=-x*a*v");
    c = c0;
    c += -z*a*b;
    Assert(Norm(c-(c0-z*prod)) <= x*eps2,label+" b+=-z*a*v");
    c = c0;
    c -= z*a*b;
    Assert(Norm(c-(c0-z*prod)) <= x*eps2,label+" b-=z*a*v");
    c = c0;
    c -= -z*a*b;
    Assert(Norm(c-(c0+z*prod)) <= x*eps2,label+" b-=-z*a*v");
    c = c0;
#endif
  }

#ifndef NOMIX
  if (CanMultMV(a,v,b)) {
    tmv::Vector<T> prod = m*v;
    if (XXDEBUG5) {
      std::cout<<"CanMult("<<tmv::TypeText(a)<<","<<tmv::TypeText(v)<<","<<tmv::TypeText(b)<<")\n";
      std::cout<<"v = "<<v<<std::endl;
      std::cout<<"m = "<<m<<std::endl;
      std::cout<<"prod = "<<prod<<std::endl;
    }
    b = a*v;
    Assert(Norm(VEC(T,b)-prod) <= eps,label+" b=a*v");
    b = v;
    b += a*v;
    Assert(Norm(VEC(T,b)-(v+prod)) <= eps2,label+" b+=a*v");
    b = v;
    RealType(T) x(5);
    T z; SetZ(z);
    b = x*a*v;
    Assert(Norm(VEC(T,b)-x*prod) <= x*eps,label+" b=x*a*v");
    b = z*a*v;
    Assert(Norm(VEC(T,b)-z*prod) <= x*eps,label+" b=z*a*v");
    b = v;
#ifdef XTEST
    b += x*a*v;
    Assert(Norm(VEC(T,b)-(v+x*prod)) <= x*eps2,label+" b+=x*a*v");
    b = v;
    b += z*a*v;
    Assert(Norm(VEC(T,b)-(v+z*prod)) <= x*eps2,label+" b+=z*a*v");
    b = v;
    b = -a*v;
    Assert(Norm(VEC(T,b)-(-prod)) <= eps,label+" b=-a*v");
    b = -x*a*v;
    Assert(Norm(VEC(T,b)-(-x*prod)) <= x*eps,label+" b=-x*a*v");
    b = -z*a*v;
    Assert(Norm(VEC(T,b)-(-z*prod)) <= x*eps,label+" b=-z*a*v");
    b = v;
    b += -a*v;
    Assert(Norm(VEC(T,b)-(v-prod)) <= eps2,label+" b+=-a*v");
    b = v;
    b -= a*v;
    Assert(Norm(VEC(T,b)-(v-prod)) <= eps2,label+" b-=a*v");
    b = v;
    b -= -a*v;
    Assert(Norm(VEC(T,b)-(v+prod)) <= eps2,label+" b-=-a*v");
    b = v;
    b += -x*a*v;
    Assert(Norm(VEC(T,b)-(v-x*prod)) <= x*eps2,label+" b+=-x*a*v");
    b = v;
    b -= x*a*v;
    Assert(Norm(VEC(T,b)-(v-x*prod)) <= x*eps2,label+" b-=x*a*v");
    b = v;
    b -= -x*a*v;
    Assert(Norm(VEC(T,b)-(v+x*prod)) <= x*eps2,label+" b-=-x*a*v");
    b = v;
    b += -z*a*v;
    Assert(Norm(VEC(T,b)-(v-z*prod)) <= x*eps2,label+" b+=-z*a*v");
    b = v;
    b -= z*a*v;
    Assert(Norm(VEC(T,b)-(v-z*prod)) <= x*eps2,label+" b-=z*a*v");
    b = v;
    b -= -z*a*v;
    Assert(Norm(VEC(T,b)-(v+z*prod)) <= x*eps2,label+" b-=-z*a*v");
    b = v;
#endif
  }
#endif

#ifdef ALIASOK
  if (CanMultMV(a,b,b)) {
    tmv::Vector<T> prod = m*v;
    b = a*b;
    if (XXDEBUG5) {
      std::cout<<"CanMult("<<tmv::TypeText(a)<<","<<tmv::TypeText(b)<<","<<tmv::TypeText(b)<<")\n";
      std::cout<<"b = a*b = "<<b<<std::endl;
      std::cout<<"b-prod = "<<(VEC(T,b)-prod)<<std::endl;
      std::cout<<"Norm(b-prod) = "<<Norm(VEC(T,b)-prod)<<std::endl;
      std::cout<<"eps = "<<eps<<std::endl;
    }
    Assert(Norm(VEC(T,b)-prod) <= eps,label+" b=a*b");
    b = v;
    b += a*b;
    Assert(Norm(VEC(T,b)-(v+prod)) <= eps2,label+" b+=a*b");
    b = v;
    b -= a*b;
    Assert(Norm(VEC(T,b)-(v-prod)) <= eps2,label+" b-=a*b");
    b = v;
    RealType(T) x(5);
    T z; SetZ(z);
    b = x*a*b;
    Assert(Norm(VEC(T,b)-x*prod) <= x*eps,label+" b=x*a*b");
    b = v;
    b = z*a*b;
    Assert(Norm(VEC(T,b)-z*prod) <= x*eps,label+" b=z*a*b");
    b = v;
#ifdef XTEST
    b += x*a*b;
    Assert(Norm(VEC(T,b)-(v+x*prod)) <= x*eps2,label+" b+=x*a*b");
    b = v;
    b += z*a*b;
    Assert(Norm(VEC(T,b)-(v+z*prod)) <= x*eps2,label+" b+=z*a*b");
    b = v;
    b = -a*b;
    Assert(Norm(VEC(T,b)-(-prod)) <= eps,label+" b=-a*b");
    b = v;
    b = -x*a*b;
    Assert(Norm(VEC(T,b)-(-x*prod)) <= x*eps,label+" b=-x*a*b");
    b = v;
    b = -z*a*b;
    Assert(Norm(VEC(T,b)-(-z*prod)) <= x*eps,label+" b=-z*a*b");
    b = v;
    b += -a*b;
    Assert(Norm(VEC(T,b)-(v-prod)) <= eps2,label+" b+=-a*b");
    b = v;
    b -= -a*b;
    Assert(Norm(VEC(T,b)-(v+prod)) <= eps2,label+" b-=-a*b");
    b = v;
    b += -x*a*b;
    Assert(Norm(VEC(T,b)-(v-x*prod)) <= x*eps2,label+" b+=-x*a*b");
    b = v;
    b -= x*a*b;
    Assert(Norm(VEC(T,b)-(v-x*prod)) <= x*eps2,label+" b-=x*a*b");
    b = v;
    b -= -x*a*b;
    Assert(Norm(VEC(T,b)-(v+x*prod)) <= x*eps2,label+" b-=-x*a*b");
    b = v;
    b += -z*a*b;
    Assert(Norm(VEC(T,b)-(v-z*prod)) <= x*eps2,label+" b+=-z*a*b");
    b = v;
    b -= z*a*b;
    Assert(Norm(VEC(T,b)-(v-z*prod)) <= x*eps2,label+" b-=z*a*b");
    b = v;
    b -= -z*a*b;
    Assert(Norm(VEC(T,b)-(v+z*prod)) <= x*eps2,label+" b-=-z*a*b");
    b = v;
#endif
  }
#endif // ALIAS
  if (showstartdone)
    std::cout<<"Done MV2a"<<std::endl;
}

template <class Ta, class T, class MM, class V> static void DoTestVM2a(
    const MM& a, CONST V& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start VM2a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
  }

  tmv::Matrix<Ta> m = a;
  tmv::Vector<T> v = b;

  double eps = EPS * Norm(v) * Norm(m);
  double eps2 = EPS * Norm(v) * (1.+Norm(m));

  if (CanMultVM(b,a)) {
    tmv::Vector<T> prod = v*m;
    tmv::Vector<T> c = prod;
    tmv::Vector<T> c0 = c;
    c0 /= T(2);

    if (XXDEBUG5) {
      std::cout<<"CanMult("<<tmv::TypeText(b)<<","<<tmv::TypeText(a)<<")\n";
      std::cout<<"v = "<<v<<std::endl;
      std::cout<<"m = "<<m<<std::endl;
      std::cout<<"prod = "<<prod<<std::endl;
    }
    c = b*a;
    Assert(Norm(c-prod) <= eps,label+" c=b*a");
    c = c0;
    c += b*a;
    Assert(Norm(c-(c0+prod)) <= eps2,label+" c+=b*a");
    c = c0;
    RealType(T) x(5);
    T z; SetZ(z);
    c = x*b*a;
    Assert(Norm(c-x*prod) <= x*eps,label+" c=x*b*a");
    c = z*b*a;
    Assert(Norm(c-z*prod) <= x*eps,label+" c=z*b*a");
    c = c0;
#ifdef XTEST
    c += -b*a;
    Assert(Norm(c-(c0-prod)) <= eps2,label+" c+=-b*a");
    c = c0;
    c += x*b*a;
    Assert(Norm(c-(c0+x*prod)) <= x*eps2,label+" c+=x*b*a");
    c = c0;
    c += z*b*a;
    Assert(Norm(c-(c0+z*prod)) <= x*eps2,label+" c+=z*b*a");
    c = c0;
    c = -b*a;
    Assert(Norm(c-(-prod)) <= eps,label+" c=-b*a");
    c = -x*b*a;
    Assert(Norm(c-(-x*prod)) <= x*eps,label+" c=-x*b*a");
    c = -z*b*a;
    Assert(Norm(c-(-z*prod)) <= x*eps,label+" c=-z*b*a");
    c = c0;
    c -= b*a;
    Assert(Norm(c-(c0-prod)) <= eps2,label+" c-=b*a");
    c = c0;
    c -= -b*a;
    Assert(Norm(c-(c0+prod)) <= eps2,label+" c-=-b*a");
    c = c0;
    c += -x*b*a;
    Assert(Norm(c-(c0-x*prod)) <= x*eps2,label+" c+=-x*b*a");
    c = c0;
    c -= x*b*a;
    Assert(Norm(c-(c0-x*prod)) <= x*eps2,label+" c-=x*b*a");
    c = c0;
    c -= -x*b*a;
    Assert(Norm(c-(c0+x*prod)) <= x*eps2,label+" c-=-x*b*a");
    c = c0;
    c += -z*b*a;
    Assert(Norm(c-(c0-z*prod)) <= x*eps2,label+" c+=-z*b*a");
    c = c0;
    c -= z*b*a;
    Assert(Norm(c-(c0-z*prod)) <= x*eps2,label+" c-=z*b*a");
    c = c0;
    c -= -z*b*a;
    Assert(Norm(c-(c0+z*prod)) <= x*eps2,label+" c-=-z*b*a");
    c = c0;
#endif
  }

#ifndef NOMIX
  if (CanMultVM(v,a,b)) {
    tmv::Vector<T> prod = v*m;
    if (XXDEBUG5) {
      std::cout<<"CanMult("<<tmv::TypeText(v)<<","<<tmv::TypeText(a)<<","<<tmv::TypeText(b)<<")\n";
      std::cout<<"v = "<<v<<std::endl;
      std::cout<<"m = "<<m<<std::endl;
      std::cout<<"prod = "<<prod<<std::endl;
    }
    b = v*a;
    Assert(Norm(b-prod) <= eps,label+" b=v*a");
    b = v;
    b += v*a;
    Assert(Norm(b-(v+prod)) <= eps2,label+" b+=v*a");
    b = v;
    RealType(T) x(5);
    T z; SetZ(z);
    b = x*v*a;
    Assert(Norm(b-x*prod) <= x*eps,label+" b=x*v*a");
    b = z*v*a;
    Assert(Norm(b-z*prod) <= x*eps,label+" b=z*v*a");
    b = v;
#ifdef XTEST
    b += -v*a;
    Assert(Norm(b-(v-prod)) <= eps2,label+" b+=-v*a");
    b = v;
    b += x*v*a;
    Assert(Norm(b-(v+x*prod)) <= x*eps2,label+" b+=x*v*a");
    b = v;
    b += z*v*a;
    Assert(Norm(b-(v+z*prod)) <= x*eps2,label+" b+=z*v*a");
    b = v;
    b = -v*a;
    Assert(Norm(b-(-prod)) <= eps,label+" b=-v*a");
    b = -x*v*a;
    Assert(Norm(b-(-x*prod)) <= x*eps,label+" b=-x*v*a");
    b = -z*v*a;
    Assert(Norm(b-(-z*prod)) <= x*eps,label+" b=-z*v*a");
    b = v;
    b -= v*a;
    Assert(Norm(b-(v-prod)) <= eps2,label+" b-=v*a");
    b = v;
    b -= -v*a;
    Assert(Norm(b-(v+prod)) <= eps2,label+" b-=-v*a");
    b = v;
    b += -x*v*a;
    Assert(Norm(b-(v-x*prod)) <= x*eps2,label+" b+=-x*v*a");
    b = v;
    b -= x*v*a;
    Assert(Norm(b-(v-x*prod)) <= x*eps2,label+" b-=x*v*a");
    b = v;
    b -= -x*v*a;
    Assert(Norm(b-(v+x*prod)) <= x*eps2,label+" b-=-x*v*a");
    b = v;
    b += -z*v*a;
    Assert(Norm(b-(v-z*prod)) <= x*eps2,label+" b+=-z*v*a");
    b = v;
    b -= z*v*a;
    Assert(Norm(b-(v-z*prod)) <= x*eps2,label+" b-=z*v*a");
    b = v;
    b -= -z*v*a;
    Assert(Norm(b-(v+z*prod)) <= x*eps2,label+" b-=-z*v*a");
    b = v;
#endif
  }
#endif

  if (CanMultVM(b,a,b)) {
    tmv::Vector<T> prod = v*m;
    if (XXDEBUG5) {
      std::cout<<"CanMult("<<tmv::TypeText(b)<<","<<tmv::TypeText(a)<<","<<tmv::TypeText(b)<<")\n";
      std::cout<<"v = "<<v<<std::endl;
      std::cout<<"m = "<<m<<std::endl;
      std::cout<<"prod = "<<prod<<std::endl;
    }
    b *= a;
    Assert(Norm(VEC(T,b)-prod) <= eps,label+" b*=a");
    b = v;
#ifdef ALIASOK
    b = b*a;
    if (XXDEBUG5) {
      std::cout<<"b = b*a = "<<b<<std::endl;
      std::cout<<"b-prod = "<<(VEC(T,b)-prod)<<std::endl;
      std::cout<<"Norm(b-prod) = "<<Norm(VEC(T,b)-prod)<<std::endl;
      std::cout<<"eps = "<<eps<<std::endl;
    }
    Assert(Norm(VEC(T,b)-prod) <= eps,label+" b=b*a");
    b = v;
    b += b*a;
    Assert(Norm(VEC(T,b)-(v+prod)) <= eps2,label+" b+=b*a");
    b = v;
#endif
#ifdef XTEST
    RealType(T) x(5);
    T z; SetZ(z);
    b *= x*a;
    Assert(Norm(VEC(T,b)-x*prod) <= x*eps,label+" b*=(x*a)");
    b = v;
    b *= z*a;
    Assert(Norm(VEC(T,b)-z*prod) <= x*eps,label+" b*=(z*a)");
    b = v;
    b *= -a;
    Assert(Norm(VEC(T,b)-(-prod)) <= eps,label+" b*=-a");
    b = v;
#ifdef ALIASOK
    b = x*b*a;
    Assert(Norm(VEC(T,b)-x*prod) <= x*eps,label+" b=x*b*a");
    b = v;
    b = z*b*a;
    Assert(Norm(VEC(T,b)-z*prod) <= x*eps,label+" b=z*b*a");
    b = v;
    b += x*b*a;
    Assert(Norm(VEC(T,b)-(v+x*prod)) <= x*eps2,label+" b+=x*b*a");
    b = v;
    b += z*b*a;
    Assert(Norm(VEC(T,b)-(v+z*prod)) <= x*eps2,label+" b+=z*b*a");
    b = v;
    b = -b*a;
    Assert(Norm(VEC(T,b)-(-prod)) <= eps,label+" b=-b*a");
    b = v;
    b = -x*b*a;
    Assert(Norm(VEC(T,b)-(-x*prod)) <= x*eps,label+" b=-x*b*a");
    b = v;
    b = -z*b*a;
    Assert(Norm(VEC(T,b)-(-z*prod)) <= x*eps,label+" b=-z*b*a");
    b = v;
    b += -b*a;
    Assert(Norm(VEC(T,b)-(v-prod)) <= eps2,label+" b+=-b*a");
    b = v;
    b -= b*a;
    Assert(Norm(VEC(T,b)-(v-prod)) <= eps2,label+" b-=b*a");
    b = v;
    b -= -b*a;
    Assert(Norm(VEC(T,b)-(v+prod)) <= eps2,label+" b-=-b*a");
    b = v;
    b += -x*b*a;
    Assert(Norm(VEC(T,b)-(v-x*prod)) <= x*eps2,label+" b+=-x*b*a");
    b = v;
    b -= x*b*a;
    Assert(Norm(VEC(T,b)-(v-x*prod)) <= x*eps2,label+" b-=x*b*a");
    b = v;
    b -= -x*b*a;
    Assert(Norm(VEC(T,b)-(v+x*prod)) <= x*eps2,label+" b-=-x*b*a");
    b = v;
    b += -z*b*a;
    Assert(Norm(VEC(T,b)-(v-z*prod)) <= x*eps2,label+" b+=-z*b*a");
    b = v;
    b -= z*b*a;
    Assert(Norm(VEC(T,b)-(v-z*prod)) <= x*eps2,label+" b-=z*b*a");
    b = v;
    b -= -z*b*a;
    Assert(Norm(VEC(T,b)-(v+z*prod)) <= x*eps2,label+" b-=-z*b*a");
    b = v;
#endif
#endif
  }
  if (showstartdone)
    std::cout<<"Done VM2a"<<std::endl;
}

template <class Ta, class T, class MM, class V> static void DoTestMV2R(
    const MM& a, CONST V& b, std::string label)
{
  DoTestMV2a<Ta,T>(a,b,label);
#ifndef NOVIEWS
  DoTestMV2a<Ta,T>(a,b.Reverse(),label+" Rev");
#endif

#ifdef XTEST
  tmv::Vector<T> b0 = b;
  b.Zero();
  DoTestMV2a<Ta,T>(a,b,label+" 1");
  b = b0;

#ifndef NOVIEWS
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
#endif
}

template <class Ta, class T, class MM, class V> static void DoTestVM2R(
    const MM& a, CONST V& b, std::string label)
{
  DoTestVM2a<Ta,T>(a,b,label);
#ifndef NOVIEWS
  DoTestMV2a<Ta,T>(a,b.Reverse(),label+" Rev");
#endif

#ifdef XTEST
  tmv::Vector<T> b0 = b;
  b.Zero();
  DoTestVM2a<Ta,T>(a,b,label+" 1");
  b = b0;

#ifndef NOVIEWS
  b.SubVector(0,b.size()/2).Zero();
  DoTestVM2a<Ta,T>(a,b,label+" 2");
  b = b0;

  b.SubVector(b.size()/2,b.size()).Zero();
  DoTestVM2a<Ta,T>(a,b,label+" 3");
  b = b0;

  b.SubVector(0,b.size()/4).Zero();
  b.SubVector(3*b.size()/4,b.size()).Zero();
  DoTestVM2a<Ta,T>(a,b,label+" 4");
  b = b0;
#endif
#endif
}

template <class Ta, class T, class MM, class V> static void DoTestMV2C(
    const MM& a, CONST V& b, std::string label)
{
  DoTestMV2a<Ta,T>(a,b,label);
#ifndef NOVIEWS
  DoTestMV2a<Ta,T>(a,b.Reverse(),label+" Rev");

  DoTestMV2a<Ta,T>(Conjugate(a),b,label+" Conj");
#endif
#ifdef XTEST
#ifndef NOVIEWS
  DoTestMV2a<Ta,T>(Conjugate(a),b.Reverse(),label+" Rev,Conj");
#endif

  tmv::Vector<T> b0 = b;
  b.Zero();
  DoTestMV2a<Ta,T>(a,b,label+" 1");
#ifndef NOVIEWS
  DoTestMV2a<Ta,T>(Conjugate(a),b,label+" Conj1");
#endif
  b = b0;

#ifndef NOVIEWS
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
#endif
}

template <class Ta, class T, class MM, class V> static void DoTestVM2C(
    const MM& a, CONST V& b, std::string label)
{
  DoTestVM2a<Ta,T>(a,b,label);
#ifndef NOVIEWS
  DoTestVM2a<Ta,T>(a,b.Reverse(),label+" Rev");

  DoTestVM2a<Ta,T>(Conjugate(a),b,label+" Conj");
#endif
#ifdef XTEST
#ifndef NOVIEWS
  DoTestVM2a<Ta,T>(Conjugate(a),b.Reverse(),label+" Rev,Conj");
#endif

  tmv::Vector<T> b0 = b;
  b.Zero();
  DoTestVM2a<Ta,T>(a,b,label+" 1");
#ifndef NOVIEWS
  DoTestVM2a<Ta,T>(Conjugate(a),b,label+" Conj1");
#endif
  b = b0;

#ifndef NOVIEWS
  b.SubVector(0,b.size()/2).Zero();
  DoTestVM2a<Ta,T>(a,b,label+" 2");
  DoTestVM2a<Ta,T>(Conjugate(a),b,label+" Conj2");
  b = b0;

  b.SubVector(b.size()/2,b.size()).Zero();
  DoTestVM2a<Ta,T>(a,b,label+" 3");
  DoTestVM2a<Ta,T>(Conjugate(a),b,label+" Conj3");
  b = b0;

  b.SubVector(0,b.size()/4).Zero();
  b.SubVector(3*b.size()/4,b.size()).Zero();
  DoTestVM2a<Ta,T>(a,b,label+" 4");
  DoTestVM2a<Ta,T>(Conjugate(a),b,label+" Conj4");
  b = b0;
#endif
#endif
}

template <class Ta, class Tb, class T, class MM, class V1, class V2> 
static void DoTestMV3a(
    const MM& a, const V1& b, CONST V2& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV3a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TypeText(c)<<" "<<c<<std::endl;
  }

  tmv::Matrix<Ta> m = a;
  tmv::Matrix<Ta> mt = Transpose(a);
  tmv::Vector<T> c0 = c;
  tmv::Vector<Tb> v1 = b;
  tmv::Vector<T> v2 = c;

  double eps = EPS * Norm(b) * Norm(a);
  double eps2 = EPS * (Norm(c0) + Norm(b) * Norm(a));

  if (XXDEBUG6) {
    std::cout<<"a = "<<TypeText(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<TypeText(b)<<"  "<<b.step()<<"  "<<b<<std::endl;
    std::cout<<"c = "<<TypeText(c)<<"  "<<c.step()<<"  "<<c<<std::endl;
    std::cout<<"a*b = "<<m*v1<<std::endl;
  }

  if (CanMultMV(a,b,c)) {
    c = a*b;
    v2 = m*v1;
    Assert(Norm(VEC(T,c)-v2) <= eps,label+" c=a*b");
    c = c0;
    c += a*b;
    v2 = c0 + m*v1;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c+=a*b");
    c = c0;
    RealType(T) x(5);
    T z; SetZ(z);
    c = x*a*b;
    v2 = x*m*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps,label+" c=x*a*b");
    c = c0;
    c = z*a*b;
    v2 = z*m*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps,label+" c=z*a*b");
    c = c0;
#ifdef XTEST
    c += x*a*b;
    v2 = c0 + x*m*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=x*a*b");
    c = c0;
    c += z*a*b;
    v2 = c0 + z*m*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=z*a*b");
    c = c0;
    c = -a*b;
    v2 = -m*v1;
    Assert(Norm(VEC(T,c)-v2) <= eps,label+" c=-a*b");
    c = c0;
    c += -a*b;
    v2 = c0 - m*v1;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c+=-a*b");
    c = c0;
    c -= a*b;
    v2 = c0 - m*v1;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c-=a*b");
    c = c0;
    c -= -a*b;
    v2 = c0 + m*v1;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c-=-a*b");
    c = c0;
    c += -x*a*b;
    v2 = c0 - x*m*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=-x*a*b");
    c = c0;
    c -= x*a*b;
    v2 = c0 - x*m*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=x*a*b");
    c = c0;
    c -= -x*a*b;
    v2 = c0 + x*m*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=-x*a*b");
    c = c0;
    c += -z*a*b;
    v2 = c0 - z*m*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=-z*a*b");
    c = c0;
    c -= z*a*b;
    v2 = c0 - z*m*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=z*a*b");
    c = c0;
    c -= -z*a*b;
    v2 = c0 + z*m*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=-z*a*b");
    c = c0;
#endif

#ifndef NOVIEWS
    c = b*Transpose(a);
    v2 = v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= eps,label+" c=b*at");
    c = c0;
    c += b*Transpose(a);
    v2 = c0 + v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c+=b*at");
    c = c0;
    c = x*b*Transpose(a);
    v2 = x*v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= x*eps,label+" c=x*b*at");
    c = c0;
    c = z*b*Transpose(a);
    v2 = z*v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= x*eps,label+" c=z*b*at");
    c = c0;
#ifdef XTEST
    c += x*b*Transpose(a);
    v2 = c0 + x*v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=x*b*at");
    c = c0;
    c += z*b*Transpose(a);
    v2 = c0 + z*v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=z*b*at");
    c = c0;
    c = -b*Transpose(a);
    v2 = -v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= eps,label+" c=-b*at");
    c = c0;
    c += -b*Transpose(a);
    v2 = c0 - v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c+=-b*at");
    c = c0;
    c -= b*Transpose(a);
    v2 = c0 - v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c-=b*at");
    c = c0;
    c -= -b*Transpose(a);
    v2 = c0 + v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c-=-b*at");
    c = c0;
    c += -x*b*Transpose(a);
    v2 = c0 - x*v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=-x*b*at");
    c = c0;
    c -= x*b*Transpose(a);
    v2 = c0 - x*v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=x*b*at");
    c = c0;
    c -= -x*b*Transpose(a);
    v2 = c0 + x*v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=-x*b*at");
    c = c0;
    c += -z*b*Transpose(a);
    v2 = c0 - z*v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=-z*b*at");
    c = c0;
    c -= z*b*Transpose(a);
    v2 = c0 - z*v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=z*b*at");
    c = c0;
    c -= -z*b*Transpose(a);
    v2 = c0 + z*v1*mt;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=-z*b*at");
    c = c0;
#endif
#endif
  }

  if (showstartdone)
    std::cout<<"Done MV3a"<<std::endl;
}

template <class Ta, class Tb, class T, class MM, class V1, class V2> 
static void DoTestVM3a(
    const MM& a, const V1& b, CONST V2& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start VM3a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TypeText(c)<<" "<<c<<std::endl;
  }

  tmv::Matrix<Ta> m = a;
  tmv::Matrix<Ta> mt = Transpose(a);
  tmv::Vector<T> c0 = c;
  tmv::Vector<Tb> v1 = b;
  tmv::Vector<T> v2 = c;

  double eps = EPS * Norm(b) * Norm(a);
  double eps2 = EPS * (Norm(c0) + Norm(b) * Norm(a));

  if (XXDEBUG6) {
    std::cout<<"a = "<<TypeText(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<TypeText(b)<<"  "<<b.step()<<"  "<<b<<std::endl;
    std::cout<<"c = "<<TypeText(c)<<"  "<<c.step()<<"  "<<c<<std::endl;
    std::cout<<"b*a = "<<v1*m<<std::endl;
  }

  if (CanMultVM(b,a,c)) {
    c = b*a;
    v2 = v1*m;
    Assert(Norm(VEC(T,c)-v2) <= eps,label+" c=b*a");
    c = c0;
    c += b*a;
    v2 = c0 + v1*m;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c+=b*a");
    c = c0;
    RealType(T) x(5);
    T z; SetZ(z);
    c = x*b*a;
    v2 = x*v1*m;
    Assert(Norm(VEC(T,c)-v2) <= x*eps,label+" c=x*b*a");
    c = c0;
    c = z*b*a;
    v2 = z*v1*m;
    Assert(Norm(VEC(T,c)-v2) <= x*eps,label+" c=z*b*a");
    c = c0;
#ifdef XTEST
    c += x*b*a;
    v2 = c0 + x*v1*m;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=x*b*a");
    c = c0;
    c += z*b*a;
    v2 = c0 + z*v1*m;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=z*b*a");
    c = c0;
    c = -b*a;
    v2 = -v1*m;
    Assert(Norm(VEC(T,c)-v2) <= eps,label+" c=-b*a");
    c = c0;
    c += -b*a;
    v2 = c0 - v1*m;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c+=-b*a");
    c = c0;
    c -= b*a;
    v2 = c0 - v1*m;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c-=b*a");
    c = c0;
    c -= -b*a;
    v2 = c0 + v1*m;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c-=-b*a");
    c = c0;
    c += -x*b*a;
    v2 = c0 - x*v1*m;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=-x*b*a");
    c = c0;
    c -= x*b*a;
    v2 = c0 - x*v1*m;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=x*b*a");
    c = c0;
    c -= -x*b*a;
    v2 = c0 + x*v1*m;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=-x*b*a");
    c = c0;
    c += -z*b*a;
    v2 = c0 - z*v1*m;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=-z*b*a");
    c = c0;
    c -= z*b*a;
    v2 = c0 - z*v1*m;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=z*b*a");
    c = c0;
    c -= -z*b*a;
    v2 = c0 + z*v1*m;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=-z*b*a");
    c = c0;
#endif

#ifndef NOVIEWS
    c = Transpose(a)*b;
    v2 = mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= eps,label+" c=at*b");
    c = c0;
    c += Transpose(a)*b;
    v2 = c0 + mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c+=at*b");
    c = c0;
    c = x*Transpose(a)*b;
    v2 = x*mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps,label+" c=x*at*b");
    c = c0;
    c = z*Transpose(a)*b;
    v2 = z*mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps,label+" c=z*at*b");
    c = c0;
#ifdef XTEST
    c += x*Transpose(a)*b;
    v2 = c0 + x*mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=x*at*b");
    c = c0;
    c += z*Transpose(a)*b;
    v2 = c0 + z*mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=z*at*b");
    c = c0;
    c = -Transpose(a)*b;
    v2 = -mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= eps,label+" c=-at*b");
    c = c0;
    c += -Transpose(a)*b;
    v2 = c0 - mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c+=-at*b");
    c = c0;
    c -= Transpose(a)*b;
    v2 = c0 - mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c-=at*b");
    c = c0;
    c -= -Transpose(a)*b;
    v2 = c0 + mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= eps2,label+" c-=-at*b");
    c = c0;
    c += -x*Transpose(a)*b;
    v2 = c0 - x*mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=-x*at*b");
    c = c0;
    c -= x*Transpose(a)*b;
    v2 = c0 - x*mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=x*at*b");
    c = c0;
    c -= -x*Transpose(a)*b;
    v2 = c0 + x*mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=-x*at*b");
    c = c0;
    c += -z*Transpose(a)*b;
    v2 = c0 - z*mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c+=-z*at*b");
    c = c0;
    c -= z*Transpose(a)*b;
    v2 = c0 - z*mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=z*at*b");
    c = c0;
    c -= -z*Transpose(a)*b;
    v2 = c0 + z*mt*v1;
    Assert(Norm(VEC(T,c)-v2) <= x*eps2,label+" c-=-z*at*b");
    c = c0;
#endif
#endif
  }

  if (showstartdone)
    std::cout<<"Done VM3a"<<std::endl;
}

template <class Ta, class Tb, class T, class MM, class V1, class V2> 
static void DoTestMV3R(
    const MM& a, const V1& b, CONST V2& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV3"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TypeText(c)<<" "<<c<<std::endl;
  }

  DoTestMV3a<Ta,Tb,T>(a,b,c,label);
#ifndef NOVIEWS
  DoTestMV3a<Ta,Tb,T>(a,b.Reverse(),c,label);
  DoTestMV3a<Ta,Tb,T>(a,b,c.Reverse(),label);
  DoTestMV3a<Ta,Tb,T>(a,b.Reverse(),c.Reverse(),label);
#endif

  if (showstartdone)
    std::cout<<"Done MV3"<<std::endl;
}

template <class Ta, class Tb, class T, class MM, class V1, class V2> 
static void DoTestVM3R(
    const MM& a, const V1& b, CONST V2& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start VM3"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TypeText(c)<<" "<<c<<std::endl;
  }

  DoTestVM3a<Ta,Tb,T>(a,b,c,label);
#ifndef NOVIEWS
  DoTestVM3a<Ta,Tb,T>(a,b.Reverse(),c,label);
  DoTestVM3a<Ta,Tb,T>(a,b,c.Reverse(),label);
  DoTestVM3a<Ta,Tb,T>(a,b.Reverse(),c.Reverse(),label);
#endif

  if (showstartdone)
    std::cout<<"Done VM3"<<std::endl;
}

template <class Ta, class Tb, class T, class MM, class V1, class V2> 
static void DoTestMV3C(
    const MM& a, const V1& b, CONST V2& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV3"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TypeText(c)<<" "<<c<<std::endl;
  }

  DoTestMV3a<Ta,Tb,T>(a,b,c,label);
#ifndef NOVIEWS
  DoTestMV3a<Ta,Tb,T>(a,b.Reverse(),c,label);
  DoTestMV3a<Ta,Tb,T>(a,b,c.Reverse(),label);
  DoTestMV3a<Ta,Tb,T>(a,b.Reverse(),c.Reverse(),label);

  DoTestMV3a<Ta,Tb,T>(Conjugate(a),b,c,label+" Conj");
#ifdef XTEST
  DoTestMV3a<Ta,Tb,T>(Conjugate(a),b.Reverse(),c,label+" Conj");
  DoTestMV3a<Ta,Tb,T>(Conjugate(a),b,c.Reverse(),label+" Conj");
  DoTestMV3a<Ta,Tb,T>(Conjugate(a),b.Reverse(),c.Reverse(),label+" Conj");
#endif
#endif

  if (showstartdone)
    std::cout<<"Done MV3"<<std::endl;
}

template <class Ta, class Tb, class T, class MM, class V1, class V2> 
static void DoTestVM3C(
    const MM& a, const V1& b, CONST V2& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start VM3"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TypeText(c)<<" "<<c<<std::endl;
  }

  DoTestVM3a<Ta,Tb,T>(a,b,c,label);
#ifndef NOVIEWS
  DoTestVM3a<Ta,Tb,T>(a,b.Reverse(),c,label);
  DoTestVM3a<Ta,Tb,T>(a,b,c.Reverse(),label);
  DoTestVM3a<Ta,Tb,T>(a,b.Reverse(),c.Reverse(),label);

  DoTestVM3a<Ta,Tb,T>(Conjugate(a),b,c,label+" Conj");
#ifdef XTEST
  DoTestVM3a<Ta,Tb,T>(Conjugate(a),b.Reverse(),c,label+" Conj");
  DoTestVM3a<Ta,Tb,T>(Conjugate(a),b,c.Reverse(),label+" Conj");
  DoTestVM3a<Ta,Tb,T>(Conjugate(a),b.Reverse(),c.Reverse(),label+" Conj");
#endif
#endif

  if (showstartdone)
    std::cout<<"Done VM3"<<std::endl;
}

template <class T, class Tb, class Tsum, IFTEMP1(class M0) IFTEMP1(class CM0) class M1, class M2> static void DoTestMM1a(
    IFTEMP1(M0& temp) IFTEMP1(CM0& ctemp) const M1& a, const M2& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM1a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
  }

  if (XXDEBUG7) {
    std::cout<<"a = "<<tmv::TypeText(a)<<" = "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" = "<<b<<std::endl;
  }

  if (CanAdd(a,b)) {
    tmv::Matrix<T> m1 = a;
    tmv::Matrix<Tb> m2 = b;
    double eps = EPS*(Norm(m1)+Norm(m2));
    {
      tmv::Matrix<Tsum> sum = m1+m2;
      tmv::Matrix<Tsum> diff = m1-m2;

      if (XXDEBUG7) {
        std::cout<<"CanAdd("<<tmv::TypeText(a)<<","<<tmv::TypeText(b)<<")\n";
        std::cout<<"m1-m2 = "<<m1-m2<<std::endl;
        std::cout<<"a-b = "<<(IFTEMP(temp=)a-b)<<std::endl;
        std::cout<<"m1+m2 = "<<m1+m2<<std::endl;
        std::cout<<"a+b = "<<(IFTEMP(temp=)a+b)<<std::endl;
      }
      Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)a-b)-diff) <= eps,label+" a-b");
      Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)a+b)-sum) <= eps,label+" a+b");
#ifdef XTEST
#ifndef NOMIX
      Assert(Norm((a-m2)-diff) <= eps,label+" a-m");
      Assert(Norm((m1-b)-diff) <= eps,label+" m-b");
      Assert(Norm((a+m2)-sum) <= eps,label+" a+m");
      Assert(Norm((m1+b)-sum) <= eps,label+" m+b");
#endif
      RealType(T) x(5);
      sum = m1+x*m2;
      diff = m1-x*m2;
      Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)a-x*b)-diff) <= x*eps,label+" a-x*b");
      Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)a+x*b)-sum) <= x*eps,label+" a+x*b");
#ifndef NOMIX
      Assert(Norm((a-x*m2)-diff) <= x*eps,label+" a-x*m");
      Assert(Norm((m1-x*b)-diff) <= x*eps,label+" m-x*b");
      Assert(Norm((a+x*m2)-sum) <= x*eps,label+" a+x*m");
      Assert(Norm((m1+x*b)-sum) <= x*eps,label+" m+x*b");
#endif
      sum = x*m1+m2;
      diff = x*m1-m2;
      Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)x*a-b)-diff) <= x*eps,label+" x*a-b");
      Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)x*a+b)-sum) <= x*eps,label+" x*a+b");
#ifndef NOMIX
      Assert(Norm((x*a-m2)-diff) <= x*eps,label+" x*a-m");
      Assert(Norm((x*m1-b)-diff) <= x*eps,label+" x*m-b");
      Assert(Norm((x*a+m2)-sum) <= x*eps,label+" x*a+m");
      Assert(Norm((x*m1+b)-sum) <= x*eps,label+" x*m+b");
#endif
      sum = x*m1+x*m2;
      diff = x*m1-x*m2;
      Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)x*a-x*b)-diff) <= x*eps,label+" x*a-x*b");
      Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)x*a+x*b)-sum) <= x*eps,label+" x*a+x*b");
#ifndef NOMIX
      Assert(Norm((x*a-x*m2)-diff) <= x*eps,label+" x*a-x*m");
      Assert(Norm((x*m1-x*b)-diff) <= x*eps,label+" x*m-x*b");
      Assert(Norm((x*a+x*m2)-sum) <= x*eps,label+" x*a+x*m");
      Assert(Norm((x*m1+x*b)-sum) <= x*eps,label+" x*m+x*b");
#endif
#endif
    }

    {
      RealType(T) x(5);
      ComplexType(T) z(3,4);
      tmv::Matrix<ComplexType(T)> sum = m1+z*m2;
      tmv::Matrix<ComplexType(T)> diff = m1-z*m2;
      Assert(Norm(MAT(ComplexType(T),IFTEMP(ctemp=)a+z*b)-sum) <= x*eps,label+" a+z*b");
      Assert(Norm(MAT(ComplexType(T),IFTEMP(ctemp=)a-z*b)-diff) <= x*eps,label+" a-z*b");
#ifdef XTEST
#ifndef NOMIX
      Assert(Norm((a-z*m2)-diff) <= x*eps,label+" a-z*m");
      Assert(Norm((m1-z*b)-diff) <= x*eps,label+" m-z*b");
      Assert(Norm((a+z*m2)-sum) <= x*eps,label+" a+z*m");
      Assert(Norm((m1+z*b)-sum) <= x*eps,label+" m+z*b");
#endif
      sum = x*m1+z*m2;
      diff = x*m1-z*m2;
      Assert(Norm(MAT(ComplexType(T),IFTEMP(ctemp=)x*a-z*b)-diff) <= x*eps,label+" x*a-z*b");
      Assert(Norm(MAT(ComplexType(T),IFTEMP(ctemp=)x*a+z*b)-sum) <= x*eps,label+" x*a+z*b");
#ifndef NOMIX
      Assert(Norm((x*a-z*m2)-diff) <= x*eps,label+" x*a-z*m");
      Assert(Norm((x*m1-z*b)-diff) <= x*eps,label+" x*m-z*b");
      Assert(Norm((x*a+z*m2)-sum) <= x*eps,label+" x*a+z*m");
      Assert(Norm((x*m1+z*b)-sum) <= x*eps,label+" x*m+z*b");
#endif
      sum = z*m1+m2;
      diff = z*m1-m2;
      Assert(Norm(MAT(ComplexType(T),IFTEMP(ctemp=)z*a-b)-diff) <= x*eps,label+" z*a-b");
      Assert(Norm(MAT(ComplexType(T),IFTEMP(ctemp=)z*a+b)-sum) <= x*eps,label+" z*a+b");
#ifndef NOMIX
      Assert(Norm((z*a-m2)-diff) <= x*eps,label+" z*a-m");
      Assert(Norm((z*m1-b)-diff) <= x*eps,label+" z*m-b");
      Assert(Norm((z*a+m2)-sum) <= x*eps,label+" z*a+m");
      Assert(Norm((z*m1+b)-sum) <= x*eps,label+" z*m+b");
#endif
      sum = z*m1+x*m2;
      diff = z*m1-x*m2;
      Assert(Norm(MAT(ComplexType(T),IFTEMP(ctemp=)z*a-x*b)-diff) <= x*eps,label+" z*a-x*b");
      Assert(Norm(MAT(ComplexType(T),IFTEMP(ctemp=)z*a+x*b)-sum) <= x*eps,label+" z*a+x*b");
#ifndef NOMIX
      Assert(Norm((z*a-x*m2)-diff) <= x*eps,label+" z*a-x*m");
      Assert(Norm((z*m1-x*b)-diff) <= x*eps,label+" z*m-x*b");
      Assert(Norm((z*a+x*m2)-sum) <= x*eps,label+" z*a+x*m");
      Assert(Norm((z*m1+x*b)-sum) <= x*eps,label+" z*m+x*b");
#endif
      sum = z*m1+z*m2;
      diff = z*m1-z*m2;
      Assert(Norm(MAT(ComplexType(T),IFTEMP(ctemp=)z*a-z*b)-diff) <= x*eps,label+" z*a-z*b");
      Assert(Norm(MAT(ComplexType(T),IFTEMP(ctemp=)z*a+z*b)-sum) <= x*eps,label+" z*a+bz*");
#ifndef NOMIX
      Assert(Norm((z*a-z*m2)-diff) <= x*eps,label+" z*a-z*m");
      Assert(Norm((z*m1-z*b)-diff) <= x*eps,label+" z*m-z*b");
      Assert(Norm((z*a+z*m2)-sum) <= x*eps,label+" z*a+z*m");
      Assert(Norm((z*m1+z*b)-sum) <= x*eps,label+" z*m+z*b");
#endif
#endif
    }
  }
  if (showstartdone)
    std::cout<<"Done MM1a"<<std::endl;
}

template <class T, IFTEMP1(class M0) IFTEMP1(class CM0) class M1, class M2> static void DoTestMM1RR(
    IFTEMP1(M0& temp) IFTEMP1(CM0& ctemp) const M1& a, const M2& b, std::string label)
{
  DoTestMM1a<T,T,T>(IFTEMP1(temp) IFTEMP1(ctemp) a,b,label);

#ifdef XTEST
#ifndef NOVIEWS
#ifdef USETEMP
  M0 temp2 = temp.Transpose();
  CM0 ctemp2 = ctemp.Transpose();
#endif
  DoTestMM1a<T,T,T>(IFTEMP1(temp2) IFTEMP1(ctemp2) Transpose(b),Transpose(a),
      label+" TransB TransA");
#endif
#endif
}

template <class T, IFTEMP1(class CM0) class M1, class M2> static void DoTestMM1RC(
    IFTEMP1(CM0& ctemp) const M1& a, const M2& b, std::string label)
{
  DoTestMM1a<T,CT,CT>(IFTEMP1(ctemp) IFTEMP1(ctemp) a,b,label);

#ifdef XTEST
#ifndef NOVIEWS
  DoTestMM1a<T,CT,CT>(IFTEMP1(ctemp) IFTEMP1(ctemp) a,Conjugate(b),
      label+" ConjB");
#ifdef USETEMP
  CM0 ctemp2 = ctemp.Transpose();
#endif
  DoTestMM1a<CT,T,CT>(IFTEMP1(ctemp2) IFTEMP1(ctemp2) Transpose(b),Transpose(a),
      label+" TransB TransA");
  DoTestMM1a<CT,T,CT>(IFTEMP1(ctemp2) IFTEMP1(ctemp2) Adjoint(b),Transpose(a),
      label+" AdjB TransA");
#endif
#endif
}

template <class T, IFTEMP1(class CM0) class M1, class M2> static void DoTestMM1CR(
    IFTEMP1(CM0& ctemp) const M1& a, const M2& b, std::string label)
{
  DoTestMM1a<CT,T,CT>(IFTEMP1(ctemp) IFTEMP1(ctemp) a,b,label);

#ifdef XTEST
#ifndef NOVIEWS
  DoTestMM1a<CT,T,CT>(IFTEMP1(ctemp) IFTEMP1(ctemp) Conjugate(a),b,
      label+" ConjA");
#ifdef USETEMP
  CM0 ctemp2 = ctemp.Transpose();
#endif
  DoTestMM1a<T,CT,CT>(IFTEMP1(ctemp2) IFTEMP1(ctemp2) Transpose(b),Transpose(a),
      label+" TransB TransA");
  DoTestMM1a<T,CT,CT>(IFTEMP1(ctemp2) IFTEMP1(ctemp2) Transpose(b),Adjoint(a),
      label+" TransB AdjA");
#endif
#endif
}

template <class T, IFTEMP1(class CM0) class M1, class M2> static void DoTestMM1CC(
    IFTEMP1(CM0& ctemp) const M1& a, const M2& b, std::string label)
{
  DoTestMM1a<CT,CT,CT>(IFTEMP1(ctemp) IFTEMP1(ctemp) a,b,label);

#ifdef XTEST
#ifndef NOVIEWS
  DoTestMM1a<CT,CT,CT>(IFTEMP1(ctemp) IFTEMP1(ctemp) Conjugate(a),b,
      label+" ConjA");
  DoTestMM1a<CT,CT,CT>(IFTEMP1(ctemp) IFTEMP1(ctemp) a,Conjugate(b),
      label+" ConjB");
  DoTestMM1a<CT,CT,CT>(IFTEMP1(ctemp) IFTEMP1(ctemp) Conjugate(a),Conjugate(b),
      label+" ConjA ConjB");
#ifdef USETEMP
  CM0 ctemp2 = ctemp.Transpose();
#endif
  DoTestMM1a<CT,CT,CT>(IFTEMP1(ctemp2) IFTEMP1(ctemp2) Transpose(b),
      Transpose(a),label+" TransB TransA");
  DoTestMM1a<CT,CT,CT>(IFTEMP1(ctemp2) IFTEMP1(ctemp2) Transpose(b),Adjoint(a),
      label+" TransB AdjA");
  DoTestMM1a<CT,CT,CT>(IFTEMP1(ctemp2) IFTEMP1(ctemp2) Adjoint(b),Transpose(a),
      label+" AdjB TransA");
  DoTestMM1a<CT,CT,CT>(IFTEMP1(ctemp2) IFTEMP1(ctemp2) Adjoint(b),Adjoint(a),
      label+" AdjA ConjB");
#endif
#endif
} 

template <class T, class Tb, class BaseM, class M1, class M2> static void DoTestMM2a(
    BaseM& a0, CONST M1& a, const M2& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM2a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"a0 = "<<tmv::TypeText(a0)<<" "<<a0<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
  }

  if (XXDEBUG8) {
    std::cout<<"a = "<<tmv::TypeText(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<"  "<<b<<std::endl;
  }
  a0 = a;
  const tmv::Matrix<T> m1 = a;
  const tmv::Matrix<Tb> m2 = b;
  if (XXDEBUG8) {
    std::cout<<"m1 = "<<tmv::TypeText(m1)<<"  "<<m1<<std::endl;
    std::cout<<"m2 = "<<tmv::TypeText(m2)<<"  "<<m2<<std::endl;
  }
  double eps = EPS*(Norm(m1)+Norm(m2));

#ifndef NOMIX
  if (CanAddEq(m1,b)) {
    if (XXDEBUG8) {
      std::cout<<"CanAddEq("<<tmv::TypeText(m1)<<","<<tmv::TypeText(b)<<")\n";
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
#endif
#ifndef NOADDEQ
  if (CanAddEq(a,b)) {
    if (XXDEBUG8) {
      std::cout<<"CanAddEq("<<tmv::TypeText(a)<<","<<tmv::TypeText(b)<<")\n";
    }
    tmv::Matrix<T> m4 = a = a0;
    a += b;
    m4 = m1+m2;
    if (XXDEBUG8) {
      std::cout<<"a += b = "<<a<<std::endl;
      std::cout<<"m4 = "<<m4<<std::endl;
    }
    Assert(Norm(MAT(T,a)-m4) <= eps,label+" a += b");
    a = a0;
#ifdef XTEST
    a += -b;
    m4 = m1-m2;
    Assert(Norm(MAT(T,a)-m4) <= eps,label+" a += -b");
    a = a0;
    a -= b;
    m4 = m1-m2;
    Assert(Norm(MAT(T,a)-m4) <= eps,label+" a -= b");
    a = a0;
    a -= -b;
    m4 = m1+m2;
    Assert(Norm(MAT(T,a)-m4) <= eps,label+" a -= -b");
    a = a0;
#endif
#ifdef ALIASOK
    a = a+b;
    m4 = m1+m2;
    Assert(Norm(MAT(T,a)-m4) <= eps,label+" a = a+b");
    a = a0;
#ifdef XTEST
    a = a-b;
    m4 = m1-m2;
    Assert(Norm(MAT(T,a)-m4) <= eps,label+" a = a-b");
    a = a0;
    a = b+a; 
    m4 = m2+m1;
    Assert(Norm(MAT(T,a)-m4) <= eps,label+" a = b+a");
    a = a0;
    a = b-a;
    m4 = m2-m1;
    Assert(Norm(MAT(T,a)-m4) <= eps,label+" a = b-a");
    a = a0;
#endif
#endif
  }
#endif // NOADDEQ

  if (showstartdone)
    std::cout<<"Done MM2a"<<std::endl;
}

template <class T, class BaseM, class M1, class M2> static void DoTestMM2RR(
    BaseM& a0, CONST M1& a, const M2& b, std::string label)
{
  DoTestMM2a<T,T>(a0,a,b,label);
}

template <class T, class BaseM, class M1, class M2> static void DoTestMM2RC(
    BaseM& a0, CONST M1& a, const M2& b, std::string label)
{
  DoTestMM2a<T,CT>(a0,a,b,label);
#ifdef XTEST
#ifndef NOVIEWS
  DoTestMM2a<T,CT>(a0,a,Conjugate(b),label+" ConjB");
#endif
#endif
}

template <class T, class BaseM, class M1, class M2> static void DoTestMM2CR(
    BaseM& a0, CONST M1& a, const M2& b, std::string label)
{
  DoTestMM2a<CT,T>(a0,a,b,label);
#ifndef NOVIEWS
  DoTestMM2a<CT,T>(a0,Conjugate(a),b,label+" ConjA");
#endif
}

template <class T, class BaseM, class M1, class M2> static void DoTestMM2CC(
    BaseM& a0, CONST M1& a, const M2& b, std::string label)
{
  DoTestMM2a<CT,CT>(a0,a,b,label);
#ifndef NOVIEWS
  DoTestMM2a<CT,CT>(a0,Conjugate(a),b,label+" ConjA");

#ifdef XTEST
  DoTestMM2a<CT,CT>(a0,a,Conjugate(b),label+" ConjB");
  DoTestMM2a<CT,CT>(a0,Conjugate(a),Conjugate(b),label+" ConjA ConjB");
#endif
#endif
}

template <class T, class Tb, IFTEMP1(class M0) class M1, class M2> static void DoTestMM3a(
    IFTEMP1(M0& temp) const M1& a, const M2& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM3a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
  }

  if (XXDEBUG7) {
    std::cout<<"a = "<<tmv::TypeText(a)<<" = "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" = "<<b<<std::endl;
  }

  if (CanMultMM(a,b)) {
    tmv::Matrix<T> m1 = a;
    tmv::Matrix<Tb> m2 = b;
    tmv::Matrix<ProductType(T,Tb)> mm = m1*m2;
    double eps = EPS*Norm(m1)*Norm(m2);
    if (XXDEBUG7) {
      std::cout<<"CanMult("<<tmv::TypeText(a)<<","<<tmv::TypeText(b)<<")\n";
      std::cout<<"m1*m2 = "<<mm<<std::endl;
      std::cout<<"m1*b = "<<(IFTEMP(temp=)m1*b)<<std::endl;
      std::cout<<"a*m2 = "<<(IFTEMP(temp=)a*m2)<<std::endl;
      std::cout<<"a*b = "<<(IFTEMP(temp=)a*b)<<std::endl;
    }
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)m1*b)-mm) <= eps,label+" m*b");
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)a*m2)-mm) <= eps,label+" a*m");
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)a*b)-mm) <= eps,label+" a*b");
  }

  if (showstartdone)
    std::cout<<"Done MM3a"<<std::endl;
}

template <class T, IFTEMP1(class M0) class M1, class M2> static void DoTestMM3RR(
    IFTEMP1(M0& temp) const M1& a, const M2& b, std::string label)
{
  DoTestMM3a<T,T>(IFTEMP1(temp) a,b,label);

#ifdef XTEST
#ifndef NOVIEWS
#ifdef USETEMP
  M0 temp2 = temp.Transpose();
#endif
  DoTestMM3a<T,T>(IFTEMP1(temp2) Transpose(b),Transpose(a),
      label+" TransB TransA");
#endif
#endif
}

template <class T, IFTEMP1(class M0) class M1, class M2> static void DoTestMM3RC(
    IFTEMP1(M0& temp) const M1& a, const M2& b, std::string label)
{
  DoTestMM3a<T,CT>(IFTEMP1(temp) a,b,label);

#ifdef XTEST
#ifndef NOVIEWS
  DoTestMM3a<T,CT>(IFTEMP1(temp) a,Conjugate(b),label+" ConjB");
#ifdef USETEMP
  M0 temp2 = temp.Transpose();
#endif
  DoTestMM3a<CT,T>(IFTEMP1(temp2) Transpose(b),Transpose(a),
      label+" TransB TransA");
  DoTestMM3a<CT,T>(IFTEMP1(temp2) Adjoint(b),Transpose(a),
      label+" AdjB TransA");
#endif
#endif
}

template <class T, IFTEMP1(class M0) class M1, class M2> static void DoTestMM3CR(
    IFTEMP1(M0& temp) const M1& a, const M2& b, std::string label)
{
  DoTestMM3a<CT,T>(IFTEMP1(temp) a,b,label);

#ifdef XTEST
#ifndef NOVIEWS
  DoTestMM3a<CT,T>(IFTEMP1(temp) Conjugate(a),b,label+" ConjA");
#ifdef USETEMP
  M0 temp2 = temp.Transpose();
#endif
  DoTestMM3a<T,CT>(IFTEMP1(temp2) Transpose(b),Transpose(a),
      label+" TransB TransA");
  DoTestMM3a<T,CT>(IFTEMP1(temp2) Transpose(b),Adjoint(a),
      label+" TransB AdjA");
#endif
#endif
}

template <class T, IFTEMP1(class M0) class M1, class M2> static void DoTestMM3CC(
    IFTEMP1(M0& temp) const M1& a, const M2& b, std::string label)
{
  DoTestMM3a<CT,CT>(IFTEMP1(temp) a,b,label);

#ifdef XTEST
#ifndef NOVIEWS
  DoTestMM3a<CT,CT>(IFTEMP1(temp) Conjugate(a),b,label+" ConjA");
  DoTestMM3a<CT,CT>(IFTEMP1(temp) a,Conjugate(b),label+" ConjB");
  DoTestMM3a<CT,CT>(IFTEMP1(temp) Conjugate(a),Conjugate(b),
      label+" ConjA ConjB");
#ifdef USETEMP
  M0 temp2 = temp.Transpose();
#endif
  DoTestMM3a<CT,CT>(IFTEMP1(temp2) Transpose(b),Transpose(a),
      label+" TransB TransA");
  DoTestMM3a<CT,CT>(IFTEMP1(temp2) Transpose(b),Adjoint(a),
      label+" TransB AdjA");
  DoTestMM3a<CT,CT>(IFTEMP1(temp2) Adjoint(b),Transpose(a),
      label+" AdjB TransA");
  DoTestMM3a<CT,CT>(IFTEMP1(temp2) Adjoint(b),Adjoint(a),
      label+" AdjA ConjB");
#endif
#endif
}

template <class T, class Tb, class BaseM, class M1, class M2> static void DoTestMM4a(
    BaseM& a0, CONST M1& a, const M2& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM4a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"a0 = "<<tmv::TypeText(a0)<<" "<<a0<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
  }

  if (XXDEBUG8) {
    std::cout<<"a = "<<tmv::TypeText(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<"  "<<b<<std::endl;
  }
  a0 = a;
  const tmv::Matrix<T> m1 = a;
  const tmv::Matrix<Tb> m2 = b;
  if (XXDEBUG8) {
    std::cout<<"m1 = "<<tmv::TypeText(m1)<<"  "<<m1<<std::endl;
    std::cout<<"m2 = "<<tmv::TypeText(m2)<<"  "<<m2<<std::endl;
  }

  double eps = EPS*Norm(m1)*Norm(m2);
  double eps2 = EPS*Norm(m1)*(1.+Norm(m2));

  if (m1.rowsize() == m2.colsize()) {
    tmv::Matrix<ProductType(T,Tb)> mm = m1*m2;
    if (CanMultMM(a,b,m1)) {
      if (XXDEBUG8) {
        std::cout<<"CanMult("<<tmv::TypeText(a)<<","<<tmv::TypeText(b)<<","<<tmv::TypeText(m1)<<")\n";
      }
      tmv::Matrix<T> m3 = m1;
      tmv::Matrix<T> m4 = m1;
      RealType(T) x(5);
      T z; SetZ(z);
      m3 = a*b;
      m4 = mm;
      Assert(Norm(m3-m4) <= eps,label+" m = a*b");
      m3 = m1;
      m3 += a*b;
      m4 = m1 + mm;
      Assert(Norm(m3-m4) <= eps2,label+" m += a*b");
      m3 = m1;
      m3 = x*a*b;
      m4 = x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m = x*a*b");
      m3 = m1;
      m3 = z*a*b;
      m4 = z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m = z*a*b");
      m3 = m1;
#ifdef XTEST
      m3 += x*a*b;
      m4 = m1 + x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m += x*a*b");
      m3 = m1;
      m3 += z*a*b;
      m4 = m1 + z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m += z*a*b");
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
      m3 = -x*a*b;
      m4 = -x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m = -x*a*b");
      m3 = m1;
      m3 += -x*a*b;
      m4 = m1 - x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m += -x*a*b");
      m3 = m1;
      m3 -= x*a*b;
      m4 = m1 - x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m -= x*a*b");
      m3 = m1;
      m3 -= -x*a*b;
      m4 = m1 + x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m -= -x*a*b");
      m3 = m1;
      m3 = -z*a*b;
      m4 = -z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m = -z*a*b");
      m3 = m1;
      m3 += -z*a*b;
      m4 = m1 - z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m += -z*a*b");
      m3 = m1;
      m3 -= z*a*b;
      m4 = m1 - z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m -= z*a*b");
      m3 = m1;
      m3 -= -z*a*b;
      m4 = m1 + z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m -= -z*a*b");
      m3 = m1;
#endif
    }
#ifndef NOMIX
    if (CanMultMM(m1,b,m1)) {
      if (XXDEBUG8) {
        std::cout<<"CanMult("<<tmv::TypeText(m1)<<","<<tmv::TypeText(b)<<","<<tmv::TypeText(m1)<<")\n";
      }
      tmv::Matrix<T> m3 = m1;
      tmv::Matrix<T> m4 = m1;
      m3 += m1*b;
      m4 = m1 + mm;
      Assert(Norm(m3-m4) <= eps2,label+" m += m1*b");
      m3 = m1;
      RealType(T) x(5);
      T z; SetZ(z);
      m3 += x*m1*b;
      m4 = m1 + x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m += x*m1*b");
      m3 = m1;
      m3 += z*m1*b;
      m4 = m1 + z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m += z*m1*b");
      m3 = m1;
#ifndef NOMULTEQ
      m3 *= b;
      m4 = mm;
      Assert(Norm(m3-m4) <= eps,label+" m *= b");
      m3 = m1;
      m3 *= x*b;
      m4 = x*mm;
      Assert(Norm(m3-m4) <= x*eps,label+" m *= x*b");
      m3 = m1;
      m3 *= z*b;
      m4 = z*mm;
      Assert(Norm(m3-m4) <= x*eps,label+" m *= z*b");
      m3 = m1;
#endif
#ifdef XTEST
      m3 += m3*b;
      m4 = m1 + mm;
      Assert(Norm(m3-m4) <= eps2,label+" m += m*b");
      m3 = m1;
      m3 += x*m3*b;
      m4 = m1 + x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m += x*m*b");
      m3 = m1;
      m3 += z*m3*b;
      m4 = m1 + z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m += z*m*b");
      m3 = m1;
#ifndef NOMULTEQ
      m3 *= -b;
      m4 = -mm;
      Assert(Norm(m3-m4) <= eps,label+" m *= -b");
      m3 = m1;
      m3 *= -x*b;
      m4 = -x*mm;
      Assert(Norm(m3-m4) <= x*eps,label+" m *= -x*b");
      m3 = m1;
      m3 *= -z*b;
      m4 = -z*mm;
      Assert(Norm(m3-m4) <= x*eps,label+" m *= -z*b");
      m3 = m1;
#endif
      m3 = m3*b;
      m4 = mm;
      Assert(Norm(m3-m4) <= eps,label+" m = m*b");
      m3 = m1;
      m3 = -m3*b;
      m4 = -mm;
      Assert(Norm(m3-m4) <= eps,label+" m = -m*b");
      m3 = m1;
      m3 = m3*(x*b);
      m4 = x*mm;
      Assert(Norm(m3-m4) <= x*eps,label+" m = m*(x*b)");
      m3 = m1;
      m3 = (x*m3)*b;
      m4 = x*mm;
      Assert(Norm(m3-m4) <= x*eps,label+" m = (x*m)*b");
      m3 = m1;
      m3 = m3*(z*b);
      m4 = z*mm;
      Assert(Norm(m3-m4) <= x*eps,label+" m = m*(z*b)");
      m3 = m1;
      m3 = (z*m3)*b;
      m4 = z*mm;
      Assert(Norm(m3-m4) <= x*eps,label+" m = (z*m)*b");
      m3 = m1;
      m3 = (x*m3)*(x*b);
      m4 = x*x*mm;
      Assert(Norm(m3-m4) <= x*x*eps,label+" m = (x*m)*(x*b)");
      m3 = m1;
      m3 = (z*m3)*(x*b);
      m4 = z*x*mm;
      Assert(Norm(m3-m4) <= x*x*eps,label+" m = (z*m)*(x*b)");
      m3 = m1;
      m3 = (x*m3)*(z*b);
      m4 = z*x*mm;
      Assert(Norm(m3-m4) <= x*x*eps,label+" m = (x*m)*(z*b)");
      m3 = m1;
      m3 = (z*m3)*(z*b);
      m4 = z*z*mm;
      Assert(Norm(m3-m4) <= x*x*eps,label+" m = (z*m)*(z*b)");
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
      m3 -= x*m1*b;
      m4 = m1 - x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m -= x*m1*b");
      m3 = m1;
      m3 -= z*m1*b;
      m4 = m1 - z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m -= z*m1*b");
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
      m3 -= x*m3*b;
      m4 = m1 - x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m -= x*m*b");
      m3 = m1;
      m3 -= z*m3*b;
      m4 = m1 - z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m -= z*m*b");
      m3 = m1;
#endif
    }
#endif
#ifndef NOMULTEQ
    if (CanMultMM(a,b,a)) {
      if (XXDEBUG8) {
        std::cout<<"CanMult("<<tmv::TypeText(a)<<","<<tmv::TypeText(b)
        <<","<<tmv::TypeText(a)<<")\n";
      }
      tmv::Matrix<T> m4 = m1;
      a *= b;
      m4 = mm;
      if (XXDEBUG8) {
        std::cout<<"a *= b = "<<a<<std::endl;
        std::cout<<"m4 = "<<m4<<std::endl;
      }
      Assert(Norm(MAT(T,a)-m4) <= eps,label+" a *= b");
      a = a0;
#ifdef ALIASOK
      a = a*b;
      m4 = mm;
      Assert(Norm(MAT(T,a)-m4) <= eps,label+" a = a*b");
      a = a0;
#endif
    }
#ifdef XTEST
    if (CanMultXMM(a,b,a)) {
      if (XXDEBUG8) {
        std::cout<<"CanMultXM("<<tmv::TypeText(a)<<","<<tmv::TypeText(b)
        <<","<<tmv::TypeText(a)<<")\n";
      }
      tmv::Matrix<T> m4 = m1;
      a *= -b;
      m4 = -mm;
      Assert(Norm(MAT(T,a)-m4) <= eps,label+" a *= -b");
      a = a0;
      RealType(T) x(5);
      T z; SetZ(z);
      a *= x*b;
      m4 = x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a *= x*b");
      a = a0;
      a *= z*b;
      m4 = z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a *= z*b");
      a = a0;
      a *= -x*b;
      m4 = -x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a *= -x*b");
      a = a0;
      a *= -z*b;
      m4 = -z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a *= -z*b");
      a = a0;
#ifdef ALIASOK
      a = -a*b;
      m4 = -mm;
      Assert(Norm(MAT(T,a)-m4) <= eps,label+" a = -a*b");
      a = a0;
      a += a*b;
      m4 = m1 + mm;
      Assert(Norm(MAT(T,a)-m4) <= eps2,label+" a += a*b");
      a = a0;
      a = (x*a)*b;
      m4 = x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a = (x*a)*b");
      a = a0;
      a += x*a*b;
      m4 = m1 + x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a += x*a*b");
      a = a0;
      a = (z*a)*b;
      m4 = z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a = (z*a)*b");
      a = a0;
      a += z*a*b;
      m4 = m1 + z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a += z*a*b");
      a = a0;
      a += -a*b;
      m4 = m1 - mm;
      Assert(Norm(MAT(T,a)-m4) <= eps2,label+" a += -a*b");
      a = a0;
      a -= a*b;
      m4 = m1 - mm;
      Assert(Norm(MAT(T,a)-m4) <= eps2,label+" a -= a*b");
      a = a0;
      a -= -a*b;
      m4 = m1 + mm;
      Assert(Norm(MAT(T,a)-m4) <= eps2,label+" a -= -a*b");
      a = a0;
      a = a*(x*b);
      m4 = x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a = a*(x*b)");
      a = a0;
      a = a*(z*b);
      m4 = z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a = a*(z*b)");
      a = a0;
      a -= x*a*b;
      m4 = m1 - x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a -= x*a*b");
      a = a0;
      a -= z*a*b;
      m4 = m1 - z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a -= z*a*b");
      a = a0;
      a = (x*a)*(x*b);
      m4 = x*x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*x*eps,label+" a = (x*a)*(x*b)");
      a = a0;
      a = (z*a)*(x*b);
      m4 = z*x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*x*eps,label+" a = (z*a)*(x*b)");
      a = a0;
      a = (x*a)*(z*b);
      m4 = x*z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*x*eps,label+" a = (x*a)*(z*b)");
      a = a0;
      a = (z*a)*(z*b);
      m4 = z*z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*x*eps,label+" a = (z*a)*(z*b)");
      a = a0;
#endif
    }
#endif
#endif // NOMULTEQ
  }
#ifndef INORDER
  if (m2.rowsize() == m1.colsize()) {
    tmv::Matrix<T> mm = m2*m1;
    if (CanMultMM(b,a,m1)) {
      if (XXDEBUG8) {
        std::cout<<"CanMult("<<tmv::TypeText(b)<<","<<tmv::TypeText(a)
        <<","<<tmv::TypeText(m1)<<")\n";
      }
      tmv::Matrix<T> m3 = m1;
      tmv::Matrix<T> m4 = m1;
      m3 = b * a;
      m4 = mm;
      Assert(Norm(m3-m4) <= eps,label+" m = b*a");
      m3 = m1;
      m3 += b * a;
      m4 = m1 + mm;
      Assert(Norm(m3-m4) <= eps2,label+" m += b*a");
      RealType(T) x(5);
      T z; SetZ(z);
      m3 = m1;
      m3 = x*b*a;
      m4 = x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m = x*b*a");
      m3 = m1;
      m3 = z*b*a;
      m4 = z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m = z*b*a");
#ifdef XTEST
      m3 = m1;
      m3 += x*b*a;
      m4 = m1 + x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m += x*b*a");
      m3 = m1;
      m3 += z*b*a;
      m4 = m1 + z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m += z*b*a");
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
      m3 = -x*b*a;
      m4 = -x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m = -x*b*a");
      m3 = m1;
      m3 += -x*b*a;
      m4 = m1 - x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m += -x*b*a");
      m3 = m1;
      m3 -= x*b*a;
      m4 = m1 - x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m -= x*b*a");
      m3 = m1;
      m3 -= -x*b*a;
      m4 = m1 + x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m -= -x*b*a");
      m3 = m1;
      m3 = -z*b*a;
      m4 = -z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m = -z*b*a");
      m3 = m1;
      m3 += -z*b*a;
      m4 = m1 - z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m += -z*b*a");
      m3 = m1;
      m3 -= z*b*a;
      m4 = m1 - z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m -= z*b*a");
      m3 = m1;
      m3 -= -z*b*a;
      m4 = m1 + z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m -= -z*b*a");
#endif
    }
#ifndef NOMIX
    if (CanMultMM(b,m1,m1)) {
      if (XXDEBUG8) {
        std::cout<<"CanMult("<<tmv::TypeText(b)<<","<<tmv::TypeText(m1)
        <<","<<tmv::TypeText(m1)<<")\n";
      }
      tmv::Matrix<T> m3 = m1;
      tmv::Matrix<T> m4 = m1;
      m3 = m1;
      m3 = b * m3;
      m4 = mm;
      Assert(Norm(m3-m4) <= eps,label+" m = b*m");
      m3 = m1;
      RealType(T) x(5);
      T z; SetZ(z);
      m3 = (x*b) * m3;
      m4 = x*mm;
      Assert(Norm(m3-m4) <= x*eps,label+" m = (x*b)*m");
      m3 = m1;
      m3 = (z*b) * m3;
      m4 = z*mm;
      Assert(Norm(m3-m4) <= x*eps,label+" m = (z*b)*m");
      m3 = m1;
#ifdef XTEST
      m3 += b * m3;
      m4 = m1 + mm;
      Assert(Norm(m3-m4) <= eps2,label+" m += b*m");
      m3 = m1;
      m3 += x*b * m3;
      m4 = m1 + x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m += x*b*m");
      m3 = m1;
      m3 += z*b * m3;
      m4 = m1 + z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m += z*b*m");
      m3 = m1;
      m3 = -b * m3;
      m4 = -mm;
      Assert(Norm(m3-m4) <= eps,label+" m = -b*m");
      m3 = m1;
      m3 = b * (x*m3);
      m4 = x*mm;
      Assert(Norm(m3-m4) <= x*eps,label+" m = b*(x*m)");
      m3 = m1;
      m3 = b * (z*m3);
      m4 = z*mm;
      Assert(Norm(m3-m4) <= x*eps,label+" m = b*(z*m)");
      m3 = m1;
      m3 = (x*b) * (x*m3);
      m4 = x*x*mm;
      Assert(Norm(m3-m4) <= x*x*eps,label+" m = (x*b)*(x*m)");
      m3 = m1;
      m3 = (z*b) * (x*m3);
      m4 = z*x*mm;
      Assert(Norm(m3-m4) <= x*x*eps,label+" m = (z*b)*(x*m)");
      m3 = m1;
      m3 = (x*b) * (z*m3);
      m4 = x*z*mm;
      Assert(Norm(m3-m4) <= x*x*eps,label+" m = (x*b)*(z*m)");
      m3 = m1;
      m3 = (z*b) * (z*m3);
      m4 = z*z*mm;
      Assert(Norm(m3-m4) <= x*x*eps,label+" m = (z*b)*(z*m)");
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
      m3 -= x*b * m3;
      m4 = m1 - x*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m -= x*b*m");
      m3 = m1;
      m3 -= z*b * m3;
      m4 = m1 - z*mm;
      Assert(Norm(m3-m4) <= x*eps2,label+" m -= z*b*m");
#endif
    }
#endif
#ifndef NOMULTEQ
#ifdef ALIASOK
    if (CanMultMM(b,a,a)) {
      if (XXDEBUG8) {
        std::cout<<"CanMult("<<tmv::TypeText(b)<<","<<tmv::TypeText(a)
        <<","<<tmv::TypeText(a)<<")\n";
      }
      a = b * a;
      Assert(Norm(MAT(T,a)-mm) <= eps,label+" a = b*a");
      a = a0;
    }
    if (CanMultXMM(b,a,a)) {
      if (XXDEBUG8) {
        std::cout<<"CanMultXM("<<tmv::TypeText(b)<<","<<tmv::TypeText(a)
        <<","<<tmv::TypeText(a)<<")\n";
      }
      tmv::Matrix<T> m4 = m1;
      a = -b * a;
      m4 = -mm;
      Assert(Norm(MAT(T,a)-m4) <= eps,label+" a = -b*a");
      a = a0;
      a += b*a;
      m4 = m1 + mm;
      Assert(Norm(MAT(T,a)-m4) <= eps2,label+" a += b*a");
      a = a0;
      RealType(T) x(5);
      T z; SetZ(z);
      a = a0;
      a = x*b*a;
      m4 = x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a = x*b*a");
      a = a0;
      a = z*b*a;
      m4 = z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a = z*b*a");
      a = a0;
#ifdef XTEST
      a = x*(b * a);
      m4 = x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a = x*(b*a)");
      a = a0;
      a += x*b*a;
      m4 = m1 + x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a += x*b*a");
      a = a0;
      a = z*(b * a);
      m4 = z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a = z*(b*a)");
      a = a0;
      a += z*b*a;
      m4 = m1 + z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a += z*b*a");
      a = a0;
      a -= b*a;
      m4 = m1 - mm;
      Assert(Norm(MAT(T,a)-m4) <= eps2,label+" a -= b*a");
      a = a0;
      a = b*(x*a);
      m4 = x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a = b*(x*a)");
      a = a0;
      a = b*(z*a);
      m4 = z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a = b*(z*a)");
      a = a0;
      a = (x*b)*(x*a);
      m4 = x*x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*x*eps,label+" a = (x*b)*(x*a)");
      a = a0;
      a = (z*b)*(x*a);
      m4 = z*x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*x*eps,label+" a = (z*b)*(x*a)");
      a = a0;
      a = (x*b)*(z*a);
      m4 = x*z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*x*eps,label+" a = (x*b)*(z*a)");
      a = a0;
      a = (z*b)*(z*a);
      m4 = z*z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*x*eps,label+" a = (z*b)*(z*a)");
      a = a0;
      a -= x*b*a;
      m4 = m1 - x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a -= x*b*a");
      a = a0;
      a -= z*b*a;
      m4 = m1 - z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a -= z*b*a");
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
    if (CanMultMM(b,b,a)) {
      if (XXDEBUG8) {
        std::cout<<"CanMult("<<tmv::TypeText(b)<<","<<tmv::TypeText(b)
        <<","<<tmv::TypeText(a)<<")\n";
      }
      tmv::Matrix<T> m4 = m1;
      a = b*b;
      m4 = mm;
      Assert(Norm(MAT(T,a)-m4) <= eps,label+" a = b*b");
      a = a0;
    }
#ifdef XTEST
    if (CanMultXMM(b,b,a)) {
      if (XXDEBUG8) {
        std::cout<<"CanMultXM("<<tmv::TypeText(b)<<","<<tmv::TypeText(b)
        <<","<<tmv::TypeText(a)<<")\n";
      }
      tmv::Matrix<T> m4 = m1;
      a = -b*b;
      m4 = -mm;
      Assert(Norm(MAT(T,a)-m4) <= eps,label+" a = -b*b");
      a = a0;
      a += b*b;
      m4 = m1 + mm;
      Assert(Norm(MAT(T,a)-m4) <= eps2,label+" a += b*b");
      a = a0;
      RealType(T) x(5);
      T z; SetZ(z);
      a = x*b*b;
      m4 = x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a = x*b*b");
      a = a0;
      a += x*b*b;
      m4 = m1 + x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a += x*b*b");
      a = a0;
      a = z*b*b;
      m4 = z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a = z*b*b");
      a = a0;
      a += z*b*b;
      m4 = m1 + z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a += z*b*b");
      a = a0;
      a += -b*b;
      m4 = m1 - mm;
      Assert(Norm(MAT(T,a)-m4) <= eps2,label+" a += -b*b");
      a = a0;
      a -= b*b;
      m4 = m1 - mm;
      Assert(Norm(MAT(T,a)-m4) <= eps2,label+" a -= b*b");
      a = a0;
      a -= -b*b;
      m4 = m1 + mm;
      Assert(Norm(MAT(T,a)-m4) <= eps2,label+" a -= -b*b");
      a = a0;
      a = -x*b*b;
      m4 = -x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a = -x*b*b");
      a = a0;
      a += -x*b*b;
      m4 = m1 - x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a += -x*b*b");
      a = a0;
      a -= x*b*b;
      m4 = m1 - x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a -= x*b*b");
      a = a0;
      a -= -x*b*b;
      m4 = m1 + x*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a -= -x*b*b");
      a = a0;
      a = -z*b*b;
      m4 = -z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps,label+" a = -z*b*b");
      a = a0;
      a += -z*b*b;
      m4 = m1 - z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a += -z*b*b");
      a = a0;
      a -= z*b*b;
      m4 = m1 - z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a -= z*b*b");
      a = a0;
      a -= -z*b*b;
      m4 = m1 + z*mm;
      Assert(Norm(MAT(T,a)-m4) <= x*eps2,label+" a -= -z*b*b");
      a = a0;
    }
#endif
  }
#endif // NOMULTEQ
#endif // INORDER

  if (showstartdone)
    std::cout<<"Done MM4a"<<std::endl;
}

template <class T, class BaseM, class M1, class M2> static void DoTestMM4RR(
    BaseM& a0, CONST M1& a, const M2& b, std::string label)
{
  DoTestMM4a<T,T>(a0,a,b,label);
}

template <class T, class BaseM, class M1, class M2> static void DoTestMM4RC(
    BaseM& a0, CONST M1& a, const M2& b, std::string label)
{
  DoTestMM4a<T,CT>(a0,a,b,label);
#ifdef XTEST
#ifndef NOVIEWS
  DoTestMM4a<T,CT>(a0,a,Conjugate(b),label+" ConjB");
#endif
#endif
}

template <class T, class BaseM, class M1, class M2> static void DoTestMM4CR(
    BaseM& a0, CONST M1& a, const M2& b, std::string label)
{
  DoTestMM4a<CT,T>(a0,a,b,label);
#ifndef NOVIEWS
  DoTestMM4a<CT,T>(a0,Conjugate(a),b,label+" ConjA");
#endif
}

template <class T, class BaseM, class M1, class M2> static void DoTestMM4CC(
    BaseM& a0, CONST M1& a, const M2& b, std::string label)
{
  DoTestMM4a<CT,CT>(a0,a,b,label);
#ifndef NOVIEWS
  DoTestMM4a<CT,CT>(a0,Conjugate(a),b,label+" ConjA");

#ifdef XTEST
  DoTestMM4a<CT,CT>(a0,a,Conjugate(b),label+" ConjB");
  DoTestMM4a<CT,CT>(a0,Conjugate(a),Conjugate(b),label+" ConjA ConjB");
#endif
#endif
}

template <class Ta, class Tb, class T, class M1, class M2, class M3> 
static void DoTestMM5a(
    const M1& a, const M2& b, CONST M3& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM5a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TypeText(c)<<" "<<c<<std::endl;
  }

  tmv::Matrix<Ta> m = a;
  tmv::Matrix<Ta> mt = Transpose(a);
  tmv::Matrix<T> c0 = c;
  tmv::Matrix<Tb> m2 = b;
  tmv::Matrix<T> m3 = c;

  double eps = EPS * Norm(b) * Norm(a);
  double eps2 = EPS * (Norm(c0) + Norm(b) * Norm(a));

  if (XXDEBUG9) {
    std::cout<<"a = "<<TypeText(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<TypeText(b)<<"  "<<b<<std::endl;
    std::cout<<"c = "<<TypeText(c)<<"  "<<c<<std::endl;
  }

  if (CanMultMM(a,b,c)) {
    c = a*b;
    m3 = m*m2;
    if (XXDEBUG9) {
      std::cout<<"c = a*b = "<<c<<std::endl;
      std::cout<<"m*m2 = "<<m3<<std::endl;
    }
    Assert(Norm(MAT(T,c)-m3) <= eps,label+" c=a*b");
    c = c0;
    c += a*b;
    m3 = c0 + m*m2;
    if (XXDEBUG9) {
      std::cout<<"c += a*b = "<<c<<std::endl;
      std::cout<<"c0 + m*m2 = "<<m3<<std::endl;
    }
    Assert(Norm(MAT(T,c)-m3) <= eps2,label+" c+=a*b");
    c = c0;
    RealType(T) x(5);
    T z; SetZ(z);
    c = x*a*b;
    m3 = x*m*m2;
    Assert(Norm(MAT(T,c)-m3) <= x*eps,label+" c=x*a*b");
    c = c0;
    c = z*a*b;
    m3 = z*m*m2;
    Assert(Norm(MAT(T,c)-m3) <= x*eps,label+" c=z*a*b");
    c = c0;
#ifdef XTEST
    c += x*a*b;
    m3 = c0 + x*m*m2;
    Assert(Norm(MAT(T,c)-m3) <= x*eps2,label+" c+=x*a*b");
    c = c0;
    c += z*a*b;
    m3 = c0 + z*m*m2;
    Assert(Norm(MAT(T,c)-m3) <= x*eps2,label+" c+=z*a*b");
    c = c0;
    c = -a*b;
    m3 = -m*m2;
    Assert(Norm(MAT(T,c)-m3) <= eps,label+" c=-a*b");
    c = c0;
    c += -a*b;
    m3 = c0 - m*m2;
    Assert(Norm(MAT(T,c)-m3) <= eps2,label+" c+=-a*b");
    c = c0;
    c -= a*b;
    m3 = c0 - m*m2;
    Assert(Norm(MAT(T,c)-m3) <= eps2,label+" c-=a*b");
    c = c0;
    c -= -a*b;
    m3 = c0 + m*m2;
    Assert(Norm(MAT(T,c)-m3) <= eps2,label+" c-=-a*b");
    c = c0;
    c += -x*a*b;
    m3 = c0 - x*m*m2;
    Assert(Norm(MAT(T,c)-m3) <= x*eps2,label+" c+=-x*a*b");
    c = c0;
    c -= x*a*b;
    m3 = c0 - x*m*m2;
    Assert(Norm(MAT(T,c)-m3) <= x*eps2,label+" c-=x*a*b");
    c = c0;
    c -= -x*a*b;
    m3 = c0 + x*m*m2;
    Assert(Norm(MAT(T,c)-m3) <= x*eps2,label+" c-=-x*a*b");
    c = c0;
    c += -z*a*b;
    m3 = c0 - z*m*m2;
    Assert(Norm(MAT(T,c)-m3) <= x*eps2,label+" c+=-z*a*b");
    c = c0;
    c -= z*a*b;
    m3 = c0 - z*m*m2;
    Assert(Norm(MAT(T,c)-m3) <= x*eps2,label+" c-=z*a*b");
    c = c0;
    c -= -z*a*b;
    m3 = c0 + z*m*m2;
    Assert(Norm(MAT(T,c)-m3) <= x*eps2,label+" c-=-z*a*b");
    c = c0;
#endif
  }

  if (showstartdone)
    std::cout<<"Done MM5a"<<std::endl;
}

template <class Ta, class Tb, class T, class M1, class M2, class M3> 
static void DoTestMM5R(
    const M1& a, const M2& b, CONST M3& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM5"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TypeText(c)<<" "<<c<<std::endl;
  }

  DoTestMM5a<Ta,Tb,T>(a,b,c,label);

  if (showstartdone)
    std::cout<<"Done MM5"<<std::endl;
}

template <class Ta, class Tb, class T, class M1, class M2, class M3> 
static void DoTestMM5C(
    const M1& a, const M2& b, CONST M3& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM5"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TypeText(c)<<" "<<c<<std::endl;
  }

  DoTestMM5a<Ta,Tb,T>(a,b,c,label);
#ifndef NOVIEWS
  DoTestMM5a<Ta,Tb,T>(Conjugate(a),b,c,label+" Conj");
#ifdef XTEST
  DoTestMM5a<Ta,Tb,T>(a,Conjugate(b),c,label+" Conj");
  DoTestMM5a<Ta,Tb,T>(a,b,Conjugate(c),label+" Conj");
  DoTestMM5a<Ta,Tb,T>(a,Conjugate(b),Conjugate(c),label+" Conj");
  DoTestMM5a<Ta,Tb,T>(Conjugate(a),Conjugate(b),c,label+" Conj");
  DoTestMM5a<Ta,Tb,T>(Conjugate(a),b,Conjugate(c),label+" Conj");
  DoTestMM5a<Ta,Tb,T>(Conjugate(a),Conjugate(b),Conjugate(c),label+" Conj");
#endif
#endif

  if (showstartdone)
    std::cout<<"Done MM5"<<std::endl;
}

template <class T, class BaseM, class M, class V1, class V2> static void DoTestOProda(
    BaseM& a0, CONST M& a, const V1& v1, const V2& v2, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start OProd"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"v1 = "<<tmv::TypeText(v1)<<" "<<v1<<std::endl;
    std::cout<<"v2 = "<<tmv::TypeText(v2)<<" "<<v2<<std::endl;
  }

  a0 = a;
  tmv::Matrix<T> vv = v1^v2;

  double eps = EPS*(Norm(a0)+Norm(v1)*Norm(v2));

  a = v1^v2;
  Assert(Norm(MAT(T,a)-vv) <= eps,label+" a = v1^v2");
  a = a0;
  a += v1^v2;
  Assert(Norm(MAT(T,a)-(a0+vv)) <= eps,label+" a += v1^v2");
  a = a0;
  a -= v1^v2;
  Assert(Norm(MAT(T,a)-(a0-vv)) <= eps,label+" a -= v1^v2");
  a = a0;
  RealType(T) x(5);
  T z; SetZ(z);
  a = x * (v1^v2);
  Assert(Norm(MAT(T,a)-x*vv) <= x*eps,label+" a = x * (v1^v2)");
  a = a0;
#ifdef SYMOPROD
  if (a.issym()) {
#endif
    a = z * (v1^v2);
    Assert(Norm(MAT(T,a)-z*vv) <= x*eps,label+" a = z * (v1^v2)");
    a = a0;
#ifdef SYMOPROD
  }
#endif
#ifdef XTEST
  a = (x * v1)^v2;
  Assert(Norm(MAT(T,a)-x*vv) <= x*eps,label+" a = (x*v1) ^ v2)");
  a = a0;
  a = v1 ^ (x * v2);
  Assert(Norm(MAT(T,a)-x*vv) <= x*eps,label+" a = v1 ^ (x*v2)");
  a = a0;
  a += x * (v1^v2);
  Assert(Norm(MAT(T,a)-(a0+x*vv)) <= x*eps,label+" a += x * (v1^v2)");
  a = a0;
  a += (x * v1)^v2;
  Assert(Norm(MAT(T,a)-(a0+x*vv)) <= x*eps,label+" a += (x*v1) ^ v2)");
  a = a0;
  a += v1 ^ (x * v2);
  Assert(Norm(MAT(T,a)-(a0+x*vv)) <= x*eps,label+" a += v1 ^ (x*v2)");
  a = a0;
  a -= x * (v1^v2);
  Assert(Norm(MAT(T,a)-(a0-x*vv)) <= x*eps,label+" a -= x * (v1^v2)");
  a = a0;
  a -= (x * v1)^v2;
  Assert(Norm(MAT(T,a)-(a0-x*vv)) <= x*eps,label+" a -= (x*v1) ^ v2)");
  a = a0;
  a -= v1 ^ (x * v2);
  Assert(Norm(MAT(T,a)-(a0-x*vv)) <= x*eps,label+" a -= v1 ^ (x*v2)");
  a = a0;
#ifdef SYMOPROD
  if (a.issym()) {
#endif
    a = (z * v1)^v2;
    Assert(Norm(MAT(T,a)-z*vv) <= x*eps,label+" a = (z*v1) ^ v2)");
    a = a0;
    a = v1 ^ (z * v2);
    Assert(Norm(MAT(T,a)-z*vv) <= x*eps,label+" a = v1 ^ (z*v2)");
    a = a0;
    a += z * (v1^v2);
    Assert(Norm(MAT(T,a)-(a0+z*vv)) <= x*eps,label+" a += z * (v1^v2)");
    a = a0;
    a += (z * v1)^v2;
    Assert(Norm(MAT(T,a)-(a0+z*vv)) <= x*eps,label+" a += (z*v1) ^ v2)");
    a = a0;
    a += v1 ^ (z * v2);
    Assert(Norm(MAT(T,a)-(a0+z*vv)) <= x*eps,label+" a += v1 ^ (z*v2)");
    a = a0;
    a -= z * (v1^v2);
    Assert(Norm(MAT(T,a)-(a0-z*vv)) <= x*eps,label+" a -= z * (v1^v2)");
    a = a0;
    a -= (z * v1)^v2;
    Assert(Norm(MAT(T,a)-(a0-z*vv)) <= x*eps,label+" a -= (z*v1) ^ v2)");
    a = a0;
    a -= v1 ^ (z * v2);
    Assert(Norm(MAT(T,a)-(a0-z*vv)) <= x*eps,label+" a -= v1 ^ (z*v2)");
    a = a0;
    a = (x * v1)^(x * v2);
    Assert(Norm(MAT(T,a)-x*x*vv) <= x*x*eps,label+" a = (x*v1) ^ (x*v2))");
    a = a0;
    a += (x * v1)^(x * v2);
    Assert(Norm(MAT(T,a)-(a0+x*x*vv)) <= x*x*eps,label+" a += (x*v1) ^ (x*v2)");
    a = a0;
    a -= (x * v1)^(x * v2);
    Assert(Norm(MAT(T,a)-(a0-x*x*vv)) <= x*x*eps,label+" a -= (x*v1) ^ (x*v2)");
    a = a0;
    a = (x * v1)^(z * v2);
    Assert(Norm(MAT(T,a)-x*z*vv) <= x*x*eps,label+" a = (x*v1) ^ (z*v2))");
    a = a0;
    a += (x * v1)^(z * v2);
    Assert(Norm(MAT(T,a)-(a0+x*z*vv)) <= x*x*eps,label+" a += (x*v1) ^ (z*v2)");
    a = a0;
    a -= (x * v1)^(z * v2);
    Assert(Norm(MAT(T,a)-(a0-x*z*vv)) <= x*x*eps,label+" a -= (x*v1) ^ (z*v2)");
    a = a0;
    a = (z * v1)^(x * v2);
    Assert(Norm(MAT(T,a)-z*x*vv) <= x*x*eps,label+" a = (z*v1) ^ (x*v2))");
    a = a0;
    a += (z * v1)^(x * v2);
    Assert(Norm(MAT(T,a)-(a0+z*x*vv)) <= x*x*eps,label+" a += (z*v1) ^ (x*v2)");
    a = a0;
    a -= (z * v1)^(x * v2);
    Assert(Norm(MAT(T,a)-(a0-z*x*vv)) <= x*x*eps,label+" a -= (z*v1) ^ (x*v2)");
    a = a0;
    a = (z * v1)^(z * v2);
    Assert(Norm(MAT(T,a)-z*z*vv) <= x*x*eps,label+" a = (z*v1) ^ (z*v2))");
    a = a0;
    a += (z * v1)^(z * v2);
    Assert(Norm(MAT(T,a)-(a0+z*z*vv)) <= x*x*eps,label+" a += (z*v1) ^ (z*v2)");
    a = a0;
    a -= (z * v1)^(z * v2);
    Assert(Norm(MAT(T,a)-(a0-z*z*vv)) <= x*x*eps,label+" a -= (z*v1) ^ (z*v2)");
    a = a0;
#ifdef SYMOPROD
  }
#endif
#endif

  if (showstartdone)
    std::cout<<"Done OProd"<<std::endl;
}

template <class T, class BaseM, class M, class V1, class V2> static void DoTestOProdR(
    BaseM& a0, CONST M& a, const V1& v1, const V2& v2, std::string label)
{
  DoTestOProda<T>(a0,a,v1,v2,label);
#ifndef NOVIEWS
  DoTestOProda<T>(a0,a,v1.Reverse(),v2.Reverse(),label+" RevBC");
#ifndef SYMOPROD
  DoTestOProda<T>(a0,a,v1.Reverse(),v2,label+" RevB");
  DoTestOProda<T>(a0,a,v1,v2.Reverse(),label+" RevC");
#endif
#endif
}

template <class T, class BaseM, class M, class V1, class V2> static void DoTestOProdC(
    BaseM& a0, CONST M& a, const V1& v1, const V2& v2, std::string label)
{
  DoTestOProda<T>(a0,a,v1,v2,label);
#ifndef NOVIEWS
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
#endif
}

template <class T, class BaseM, class BaseCM, class M, class CM> 
static void TestMatrixArith1(
    BaseM& a0, BaseCM& ca0, CONST M& a, CONST CM& ca,
    std::string label)
{
  if (showstartdone) {
    std::cout<<"Start TestMatrixArith1 "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"a0 = "<<tmv::TypeText(a0)<<" "<<a0<<std::endl;
    std::cout<<"ca = "<<tmv::TypeText(ca)<<" "<<ca<<std::endl;
    std::cout<<"ca0 = "<<tmv::TypeText(ca0)<<" "<<ca0<<std::endl;
  }

  CT z(9,-2);
  T x = 12;
  DoTestMR<T>(a,label+" R");
  DoTestMC<CT>(ca,label+" C");

  DoTestMX1R<T>(IFTEMP1(a0) a,x,label+" R,R");
  DoTestMX2R<T>(a0,a,x,label+" R,R");
  DoTestMX1C<CT>(IFTEMP1(ca0) ca,z,label+" C,C");
  DoTestMX2C<CT>(ca0,ca,z,label+" C,C");
#ifdef XTEST
  DoTestMX1R<T>(IFTEMP1(ca0) a,z,label+" R,C");
  DoTestMX1C<CT>(IFTEMP1(ca0) ca,x,label+" C,R");
  DoTestMX2C<CT>(ca0,ca,x,label+" C,R");
#endif

#ifdef ALIASOK
  DoTestMM2RR<T>(a0,a,a,label+" self_arith");
  DoTestMM4RR<T>(a0,a,a,label+" self_arith");
  DoTestMM2CC<T>(ca0,ca,ca,label+" self_arith");
  DoTestMM4CC<T>(ca0,ca,ca,label+" self_arith");
#endif

  if (showstartdone)
    std::cout<<"Done Test1"<<std::endl;
}

template <class T, class M, class CM, class V1, class CV1, class V2, class CV2> 
static void TestMatrixArith2a(
    const M& a, const CM& ca, CONST V1& b, CONST CV1& cb, 
    V2& c, CV2& cc, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start TestMatrixArith2a "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"ca = "<<tmv::TypeText(ca)<<" "<<ca<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
    std::cout<<"cb = "<<tmv::TypeText(cb)<<" "<<cb<<std::endl;
    std::cout<<"c = "<<tmv::TypeText(c)<<" "<<c<<std::endl;
    std::cout<<"cc = "<<tmv::TypeText(cc)<<" "<<cc<<std::endl;
  }

  DoTestMV1R<T,T>(IFTEMP1(c) IFTEMP1(cc) a,b,label+" R,R");
  DoTestMV1C<CT,CT>(IFTEMP1(cc) IFTEMP1(cc) ca,cb,label+" C,C");
#ifdef XTEST
  DoTestMV1R<T,CT>(IFTEMP1(cc) IFTEMP1(cc) a,cb,label+" R,C");
  DoTestMV1C<CT,T>(IFTEMP1(cc) IFTEMP1(cc) ca,b,label+" C,R");
#endif

#ifndef NONSQUARE
  DoTestMV2R<T,T>(a,b,label+" R,R");
  DoTestMV2C<CT,CT>(ca,cb,label+" C,C");
#ifdef XTEST
  DoTestMV2R<T,CT>(a,cb,label+" R,C");
#endif
#endif

  if (showstartdone)
    std::cout<<"Done Test2a"<<std::endl;
}

template <class T, class M, class CM, class V1, class CV1, class V2, class CV2> 
static void TestMatrixArith2b(
    const M& a, const CM& ca, V1& b, CV1& cb, 
    CONST V2& c, CONST CV2& cc, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start TestMatrixArith2b "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"ca = "<<tmv::TypeText(ca)<<" "<<ca<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
    std::cout<<"cb = "<<tmv::TypeText(cb)<<" "<<cb<<std::endl;
    std::cout<<"c = "<<tmv::TypeText(c)<<" "<<c<<std::endl;
    std::cout<<"cc = "<<tmv::TypeText(cc)<<" "<<cc<<std::endl;
  }

  DoTestVM1R<T,T>(IFTEMP1(b) IFTEMP1(cb) a,c,label+" R,R");
  DoTestVM1C<CT,CT>(IFTEMP1(cb) IFTEMP1(cb) ca,cc,label+" C,C");
#ifdef XTEST
  DoTestVM1R<T,CT>(IFTEMP1(cb) IFTEMP1(cb) a,cc,label+" R,C");
  DoTestVM1C<CT,T>(IFTEMP1(cb) IFTEMP1(cb) ca,c,label+" C,R");
#endif

#ifndef NONSQUARE
  DoTestVM2R<T,T>(a,c,label+" R,R");
  DoTestVM2C<CT,CT>(ca,cc,label+" C,C");
#ifdef XTEST
  DoTestVM2R<T,CT>(a,cc,label+" R,C");
#endif
#endif

  if (showstartdone)
    std::cout<<"Done Test2b"<<std::endl;
}

template <class T, class M, class CM, class V1, class CV1, class V2, class CV2> 
static void TestMatrixArith3(
    const M& a, const CM& ca, CONST V1& b, CONST CV1& cb,
    CONST V2& c, CONST CV2& cc, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start TestMatrixArith3 "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"ca = "<<tmv::TypeText(ca)<<" "<<ca<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
    std::cout<<"cb = "<<tmv::TypeText(cb)<<" "<<cb<<std::endl;
    std::cout<<"c = "<<tmv::TypeText(c)<<" "<<c<<std::endl;
    std::cout<<"cc = "<<tmv::TypeText(cc)<<" "<<cc<<std::endl;
  }

  DoTestMV3R<T,T,T>(a,b,c,label+" R,R,R");
  DoTestVM3R<T,T,T>(a,c,b,label+" R,R,R");
  DoTestMV3C<CT,CT,CT>(ca,cb,cc,label+" C,C,C");
  DoTestVM3C<CT,CT,CT>(ca,cc,cb,label+" C,C,C");
#ifdef XTEST
  DoTestMV3R<T,T,CT>(a,b,cc,label+" C,R,R");
  DoTestVM3R<T,T,CT>(a,c,cb,label+" C,R,R");
  DoTestMV3R<T,CT,CT>(a,cb,cc,label+" C,R,C");
  DoTestVM3R<T,CT,CT>(a,cc,cb,label+" C,R,C");
  DoTestMV3C<CT,T,CT>(ca,b,cc,label+" C,C,R");
  DoTestVM3C<CT,T,CT>(ca,c,cb,label+" C,C,R");
#endif

  if (showstartdone)
    std::cout<<"Done Test3"<<std::endl;
}

template <class T, class BaseM, class BaseCM, class M, class CM> 
static void TestMatrixArith123(
    BaseM& a0, BaseCM& ca0, CONST M& a, CONST CM& ca,
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

  TestMatrixArith2a<T>(a,ca,vv,cvv,w,cw,label);
  TestMatrixArith2a<T>(a,ca,vs,cvs,w,cw,label);
  TestMatrixArith2b<T>(a,ca,v,cv,wv,cwv,label);
  TestMatrixArith2b<T>(a,ca,v,cv,ws,cws,label);

  TestMatrixArith3<T>(a,ca,vv,cvv,wv,cwv,label);
  TestMatrixArith3<T>(a,ca,vv,cvv,ws,cws,label);
  TestMatrixArith3<T>(a,ca,vs,cvs,wv,cwv,label);
  TestMatrixArith3<T>(a,ca,vs,cvs,ws,cws,label); 
}

template <class T, class BaseM, class BaseCM, class M1, class CM1, class M2, class CM2> 
static void TestMatrixArith4(
    BaseM& a0, BaseCM& ca0,
    CONST M1& a, CONST CM1& ca, const M2& b, const CM2& cb, 
    std::string label)
{
  if (showstartdone) {
    std::cout<<"Start TestMatrixArith4 "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
    std::cout<<"ca = "<<tmv::TypeText(ca)<<" "<<ca<<std::endl;
    std::cout<<"cb = "<<tmv::TypeText(cb)<<" "<<cb<<std::endl;
    std::cout<<"a0 = "<<tmv::TypeText(a0)<<std::endl;
    std::cout<<"ca0 = "<<tmv::TypeText(ca0)<<std::endl;
  }

  DoTestMM1RR<T>(IFTEMP1(a0) IFTEMP1(ca0) a,b,label+" R,R");
  DoTestMM2RR<T>(a0,a,b,label+" R,R");
  DoTestMM1CC<T>(IFTEMP1(ca0) ca,cb,label+" C,C");
  DoTestMM2CC<T>(ca0,ca,cb,label+" C,C");
#ifdef XTEST
  DoTestMM1RC<T>(IFTEMP1(ca0) a,cb,label+" R,C");
  DoTestMM1CR<T>(IFTEMP1(ca0) ca,b,label+" C,R");
  DoTestMM2CR<T>(ca0,ca,b,label+" C,R");
#endif

  if (showstartdone)
    std::cout<<"Done Test4"<<std::endl;
}

template <class T, class BaseM, class BaseCM, IFTEMP1(class M3) IFTEMP1(class CM3) class M1, class CM1, class M2, class CM2>
static void TestMatrixArith5(
    BaseM& a0, BaseCM& ca0, IFTEMP1(M3& c) IFTEMP1(CM3& cc) 
    CONST M1& a, CONST CM1& ca, const M2& b, const CM2& cb, 
    std::string label)
{
  if (showstartdone) {
    std::cout<<"Start TestMatrixArith5 "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
    std::cout<<"ca = "<<tmv::TypeText(ca)<<" "<<ca<<std::endl;
    std::cout<<"cb = "<<tmv::TypeText(cb)<<" "<<cb<<std::endl;
    std::cout<<"a0 = "<<tmv::TypeText(a0)<<std::endl;
    std::cout<<"ca0 = "<<tmv::TypeText(ca0)<<std::endl;
  }

  DoTestMM3RR<T>(IFTEMP1(c) a,b,label+" R,R");
  DoTestMM4RR<T>(a0,a,b,label+" R,R");
  DoTestMM3CC<T>(IFTEMP1(cc) ca,cb,label+" C,C");
  DoTestMM4CC<T>(ca0,ca,cb,label+" C,C");
#ifdef XTEST
  DoTestMM3RC<T>(IFTEMP1(cc) a,cb,label+" R,C");
  DoTestMM3CR<T>(IFTEMP1(cc) ca,b,label+" C,R");
  DoTestMM4CR<T>(ca0,ca,b,label+" C,R");
#endif

  if (showstartdone)
    std::cout<<"Done Test5"<<std::endl;
}

template <class T, class M1, class CM1, class M2, class CM2, class M3, class CM3> 
static void TestMatrixArith6(
    const M1& a, const CM1& ca, CONST M2& b, CONST CM2& cb,
    CONST M3& c, CONST CM3& cc, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start TestMatrixArith6 "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"ca = "<<tmv::TypeText(ca)<<" "<<ca<<std::endl;
    std::cout<<"b = "<<tmv::TypeText(b)<<" "<<b<<std::endl;
    std::cout<<"cb = "<<tmv::TypeText(cb)<<" "<<cb<<std::endl;
    std::cout<<"c = "<<tmv::TypeText(c)<<" "<<c<<std::endl;
    std::cout<<"cc = "<<tmv::TypeText(cc)<<" "<<cc<<std::endl;
  }

  DoTestMM5R<T,T,T>(a,b,c,label+" R,R,R");
  DoTestMM5C<CT,CT,CT>(ca,cb,cc,label+" C,C,C");
#ifdef XTEST
  DoTestMM5R<T,T,CT>(a,b,cc,label+" C,R,R");
  DoTestMM5R<T,CT,CT>(a,cb,cc,label+" C,R,C");
  DoTestMM5C<CT,T,CT>(ca,b,cc,label+" C,C,R");
#endif

  if (showstartdone)
    std::cout<<"Done Test6"<<std::endl;
}

template <class T, class BaseM, class BaseCM, class M1, class CM1, class M2, class CM2> 
static void TestMatrixArith456(
    BaseM& a0, BaseCM& ca0,
    CONST M1& a, CONST CM1& ca, const M2& b, const CM2& cb, 
    std::string label)
{
  if (CanAdd(a,b)) {
    TestMatrixArith4<T>(a0,ca0,a,ca,b,cb,label);
  }
  if (CanMultMM(a,b)) {
    tmv::Matrix<T,tmv::ColMajor> c1(a*b);
    tmv::Matrix<CT,tmv::ColMajor> cc1(ca*cb);
    TestMatrixArith5<T>(a0,ca0,IFTEMP1(c1) IFTEMP1(cc1) a,ca,b,cb,label);

    TestMatrixArith6<T>(a,ca,b,cb,c1.View(),cc1.View(),label);
#ifdef XTEST
    tmv::Matrix<T,tmv::RowMajor> c2(c1);
    tmv::Matrix<CT,tmv::RowMajor> cc2(cc1);
    TestMatrixArith6<T>(a,ca,b,cb,c2.View(),cc2.View(),label);

    tmv::Matrix<T> c3(4*c1.colsize(),4*c1.rowsize());
    tmv::Matrix<CT> cc3(4*c1.colsize(),4*c1.rowsize());
    tmv::MatrixView<T> c3v = c3.SubMatrix(0,c3.colsize(),0,c3.rowsize(),4,4);
    tmv::MatrixView<CT> cc3v = cc3.SubMatrix(0,c3.colsize(),0,c3.rowsize(),4,4);
    c3v = c1;
    cc3v = cc1;
    TestMatrixArith6<T>(a,ca,b,cb,c3v,cc3v,label);
#endif
  }
}

template <class T, class BaseM, class BaseCM, class M, class CM, class V1, class CV1, class V2, class CV2> 
static void TestMatrixArith7(
    BaseM& a0, BaseCM& ca0, CONST M& a, CONST CM& ca,
    const V1& v1, const CV1& cv1, const V2& v2, const CV2& cv2, 
    std::string label)
{
  if (showstartdone) {
    std::cout<<"Start TestMatrixArith7 "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TypeText(a)<<" "<<a<<std::endl;
    std::cout<<"ca = "<<tmv::TypeText(ca)<<" "<<ca<<std::endl;
    std::cout<<"v1 = "<<tmv::TypeText(v1)<<" "<<v1<<std::endl;
    std::cout<<"cv1 = "<<tmv::TypeText(cv1)<<" "<<cv1<<std::endl;
    std::cout<<"v2 = "<<tmv::TypeText(v2)<<" "<<v2<<std::endl;
    std::cout<<"cv2 = "<<tmv::TypeText(cv2)<<" "<<cv2<<std::endl;
    std::cout<<"a0 = "<<tmv::TypeText(a0)<<std::endl;
    std::cout<<"ca0 = "<<tmv::TypeText(ca0)<<std::endl;
  }

  DoTestOProdR<T>(a0,a,v1,v2,label+" R,R,R");
  DoTestOProdC<CT>(ca0,ca,cv1,cv2,label+" C,C,C");
#ifdef XTEST
  DoTestOProdC<CT>(ca0,ca,v1,v2,label+" C,R,R");
#ifndef SYMOPROD
  DoTestOProdC<CT>(ca0,ca,cv1,v2,label+" C,C,R");
  DoTestOProdC<CT>(ca0,ca,v1,cv2,label+" C,R,C");
#endif
#endif

  if (showstartdone)
    std::cout<<"Done Test7"<<std::endl;
}

#undef CT
#undef CONST
