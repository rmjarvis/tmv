// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#define CT std::complex<T>

template <class SM1, class SM2> inline bool CanLDiv(
    const SM1& a, const SM2& b)
{ return a.colsize() == b.colsize(); }

template <class SM1, class SM2, class SM3> inline bool CanLDiv(
    const SM1& a, const SM2& b, const SM3& c)
{ 
  return a.colsize() == b.colsize() && a.rowsize() == c.colsize() 
  && b.rowsize() == c.rowsize(); 
}

template <class SM1, class SM2> inline bool CanLDivEq(
    const SM1& a, const SM2& b)
{ return CanLDiv(a,b) && b.IsSquare(); }

template <class V1, class SM2> inline bool CanLDivVM(
    const V1& a, const SM2& b)
{ return a.size() == b.colsize(); }

template <class V1, class SM2, class V3> inline bool CanLDivVM(
    const V1& a, const SM2& b, const V3& c)
{ return a.size() == b.colsize() && c.size() == b.rowsize(); }

template <class V1, class SM2> inline bool CanLDivEqVM(
    const V1& a, const SM2& b)
{ return CanLDivVM(a,b) && b.IsSquare(); }

template <class SM1, class SM2> inline bool CanRDiv(
    const SM1& a, const SM2& b)
{ return a.rowsize() == b.rowsize(); }

template <class SM1, class SM2, class SM3> inline bool CanRDiv(
    const SM1& a, const SM2& b, const SM3& c)
{ 
  return a.rowsize() == b.rowsize() && a.colsize() == c.colsize() 
  && b.colsize() == c.rowsize(); 
}

template <class SM1, class SM2> inline bool CanRDivEq(
    const SM1& a, const SM2& b)
{ return CanRDiv(a,b) && b.IsSquare(); }

template <class V1, class SM2> inline bool CanRDivVM(
    const V1& a, const SM2& b)
{ return a.size() == b.rowsize(); }

template <class V1, class SM2, class V3> inline bool CanRDivVM(
    const V1& a, const SM2& b, const V3& c)
{ return a.size() == b.rowsize() && c.size() == b.colsize(); }

template <class V1, class SM2> inline bool CanRDivEqVM(
    const V1& a, const SM2& b)
{ return CanRDivVM(a,b) && b.IsSquare(); }

template <class T1, class T2> struct ProdType
{ typedef T1 Tprod; };

template <class T> struct ProdType<T,std::complex<T> >
{ typedef std::complex<T> Tprod; };

#define ProductType(T1,T2) typename ProdType<T1,T2>::Tprod

#ifdef NOVIEWS
#define CONST
#else
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

template <class T, class Tb, IFTEMP1(class V0) IFTEMP1(class CV0) class V, class MM> static void DoTestLDivVM1a(
    tmv::DivType dt, IFTEMP1(V0& temp) IFTEMP1(CV0& ctemp) 
    const V& a, const MM& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start LDiv VM1a: "<<label<<std::endl;
  }

  tmv::Vector<T> v = a;
  tmv::Matrix<Tb> m = b;
  m.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  m.SaveDiv();

  double eps = EPS*Norm(m)*Norm(tmv::Matrix<Tb>(m.Inverse()))*v.size();

  if (CanLDivVM(a,b)) {
    tmv::Vector<ProductType(T,Tb)> frac = v/m;
    eps *= Norm(frac);
    if (XXDEBUG1) {
      std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
      std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
      std::cout<<"v = "<<tmv::TMV_Text(v)<<"  "<<v<<std::endl;
      std::cout<<"m = "<<tmv::TMV_Text(m)<<"  "<<m<<std::endl;
      std::cout<<"a/b = "<<(IFTEMP(temp=)a/b)<<std::endl;
      std::cout<<"v/m = "<<frac<<std::endl;
      std::cout<<"Norm(diff) = "<<Norm(VEC2(T,Tb,IFTEMP(temp=)a/b)-frac)<<std::endl;
      std::cout<<"eps = "<<eps<<std::endl;
    }
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)a/b)-frac) <= eps,label+" a/b");
    TMV_RealType(T) x(5);
    TMV_ComplexType(T) z(3,4);
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)x*a/b)-x*frac) <= x*eps,label+" x*a/b");
    Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)z*a/b)-z*frac) <= x*eps,label+" z*a/b");
#ifdef XTEST
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)b.Inverse()*a)-frac) <= eps,label+" b^-1*a");
#ifndef NOMIX
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)a/m)-frac) <= eps,label+" a/m");
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)m.Inverse()*a)-frac) <= eps,label+" m^-1*a");
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)v/b)-frac) <= eps,label+" v/b");
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)b.Inverse()*v)-frac) <= eps,label+" b^-1*v");
#endif

    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)(x*a)/b)-x*frac) <= x*eps,label+" (x*a)/b");
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)x*(a/b))-x*frac) <= x*eps,label+" x*(a/b)");
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)a/(x*b))-frac/x) <= eps/x,label+" a/(x*b)");

    Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)(z*a)/b)-z*frac) <= x*eps,label+" (z*a)/b");
    Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)z*(a/b))-z*frac) <= x*eps,label+" z*(a/b)");
    Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)a/(z*b))-frac/z) <= eps/x,label+" a/(z*b)");

    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)(x*a)/(x*b))-frac) <= eps,
        label+" (x*a)/(x*b)");
    Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)(z*a)/(x*b))-(z/x)*frac) <= eps,
        label+" (z*a)/(x*b)");
    Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)(x*a)/(z*b))-(x/z)*frac) <= eps,
        label+" (x*a)/(z*b)");
    Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)(z*a)/(z*b))-frac) <= eps,
        label+" (z*a)/(z*b)");
#endif
  }

  if (showstartdone) 
    std::cout<<"Done LDiv VM1a"<<std::endl;
}

template <class T, class Tb, IFTEMP1(class V0) IFTEMP1(class CV0) class V, class MM> static void DoTestRDivVM1a(
    tmv::DivType dt, IFTEMP1(V0& temp) IFTEMP1(CV0& ctemp) 
    const V& a, const MM& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start RDiv VM1a: "<<label<<std::endl;
  }

  tmv::Vector<T> v = a;
  tmv::Matrix<Tb> m = b;
  m.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  m.SaveDiv();

  double eps = EPS*Norm(m)*Norm(tmv::Matrix<Tb>(m.Inverse()))*v.size();

  if (CanRDivVM(a,b)) {
    tmv::Vector<ProductType(T,Tb)> frac = v%m;
    eps *= Norm(frac);
    if (XXDEBUG2) {
      std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
      std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
      std::cout<<"v = "<<tmv::TMV_Text(v)<<"  "<<v<<std::endl;
      std::cout<<"m = "<<tmv::TMV_Text(m)<<"  "<<m<<std::endl;
      std::cout<<"a%b = "<<(IFTEMP(temp=)a%b)<<std::endl;
      std::cout<<"v%m = "<<frac<<std::endl;
      std::cout<<"Norm(diff) = "<<Norm(VEC2(T,Tb,IFTEMP(temp=)a%b)-frac)<<std::endl;
      std::cout<<"eps = "<<eps<<std::endl;
    }
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)a%b)-frac) <= eps,label+" a%b");
    TMV_RealType(T) x(5);
    TMV_ComplexType(T) z(3,4);
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)x*a%b)-x*frac) <= x*eps,label+" x*a%b");
    Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)z*a%b)-z*frac) <= x*eps,label+" z*a%b");
#ifdef XTEST
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)a*b.Inverse())-frac) <= eps,label+" a*b^-1");
#ifndef NOMIX
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)a%m)-frac) <= eps,label+" a%m");
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)a*m.Inverse())-frac) <= eps,label+" a*m^-1");
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)v%b)-frac) <= eps,label+" v%b");
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)v*b.Inverse())-frac) <= eps,label+" v*b^-1");
#endif

    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)(x*a)%b)-x*frac) <= x*eps,label+" (x*a)%b");
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)x*(a%b))-x*frac) <= x*eps,label+" x*(a%b)");
    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)a%(x*b))-frac/x) <= eps/x,label+" a%(x*b)");

    Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)(z*a)%b)-z*frac) <= x*eps,label+" (z*a)%b");
    Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)z*(a%b))-z*frac) <= x*eps,label+" z*(a%b)");
    Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)a%(z*b))-frac/z) <= eps/x,label+" a%(z*b)");

    Assert(Norm(VEC2(T,Tb,IFTEMP(temp=)(x*a)%(x*b))-frac) <= eps,
        label+" (x*a)%(x*b)");
    Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)(z*a)%(x*b))-(z/x)*frac) <= eps,
        label+" (z*a)%(x*b)");
    Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)(x*a)%(z*b))-(x/z)*frac) <= eps,
        label+" (x*a)%(z*b)");
    Assert(Norm(VEC(TMV_ComplexType(T),IFTEMP(ctemp=)(z*a)%(z*b))-frac) <= eps,
        label+" (z*a)%(z*b)");
#endif
  }

  if (showstartdone) 
    std::cout<<"Done RDiv VM1a"<<std::endl;
}

template <class T, class Tb, IFTEMP1(class V0) IFTEMP1(class CV0) class V, class MM> static void DoTestLDivVM1(
    tmv::DivType dt, IFTEMP1(V0& v0) IFTEMP1(CV0& cv0) 
    const V& a, const MM& b, std::string label)
{
  DoTestLDivVM1a<T,Tb>(dt,IFTEMP1(v0) IFTEMP1(cv0) a,b,label);

#ifndef NOVIEWS
  DoTestLDivVM1a<T,Tb>(dt,IFTEMP1(v0) IFTEMP1(cv0) a.Reverse(),b,
      label+" RevA");
  if (tmv::isComplex(T()))
    DoTestLDivVM1a<T,Tb>(dt,IFTEMP1(v0) IFTEMP1(cv0) a.Conjugate(),b,
        label+" ConjA");
#endif
}

template <class T, class Tb, IFTEMP1(class V0) IFTEMP1(class CV0) class V, class MM> static void DoTestRDivVM1(
    tmv::DivType dt, IFTEMP1(V0& v0) IFTEMP1(CV0& cv0) 
    const V& a, const MM& b, std::string label)
{
  DoTestRDivVM1a<T,Tb>(dt,IFTEMP1(v0) IFTEMP1(cv0) a,b,label);

#ifndef NOVIEWS
  DoTestRDivVM1a<T,Tb>(dt,IFTEMP1(v0) IFTEMP1(cv0) a.Reverse(),b,
      label+" RevA");
  if (tmv::isComplex(T()))
    DoTestRDivVM1a<T,Tb>(dt,IFTEMP1(v0) IFTEMP1(cv0) a.Conjugate(),b,
        label+" ConjA");
#endif
}

template <class T> inline void SetZ(T& z)
{ z = T(5); }
template <class T> inline void SetZ(std::complex<T>& z)
{ z = std::complex<T>(3,4); }

#ifndef NOLDIVEQ
template <class T, class Tb, class V, class MM> static void DoTestLDivVM2a(
    tmv::DivType dt, CONST V& a, const MM& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start LDiv VM2b: "<<label<<std::endl;
  }

  tmv::Vector<T> v = a;
  tmv::Matrix<Tb> m = b;
  m.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  m.SaveDiv();

  double eps = EPS*Norm(m)*Norm(tmv::Matrix<Tb>(m.Inverse()))*v.size();

  if (CanLDivEqVM(a,b)) {
    tmv::Vector<T> a0 = v;
    tmv::Vector<T> frac = v/m;
    eps *= Norm(frac);
    a /= b;
    Assert(Norm(VEC(T,a)-frac) <= eps,label+" a/=b");
    a = a0;
#ifdef ALIASOK
    a = a / b;
    Assert(Norm(VEC(T,a)-frac) <= eps,label+" a=a/b");
    a = a0;
    a = b.Inverse() * a;
    Assert(Norm(VEC(T,a)-frac) <= eps,label+" a=b^-1*a");
    a = a0;
#ifdef XTEST
    TMV_RealType(T) x(5);
    T z;  SetZ(z);
    a = x * a / b;
    Assert(Norm(VEC(T,a)-x*frac) <= x*eps,label+" a=x*a/b");
    a = a0;
    a = x * b.Inverse() * a;
    Assert(Norm(VEC(T,a)-x*frac) <= x*eps,label+" a=x*b^-1*a");
    a = a0;
    a = z * a / b;
    Assert(Norm(VEC(T,a)-z*frac) <= x*eps,label+" a=z*a/b");
    a = a0;
    a = z * b.Inverse() * a;
    Assert(Norm(VEC(T,a)-z*frac) <= x*eps,label+" a=z*b^-1*a");
    a = a0;
#endif
#endif
  }
  if (showstartdone) 
    std::cout<<"Done LDiv VM2a"<<std::endl;
}

template <class T, class Tb, class V, class MM> static void DoTestLDivVM2(
    tmv::DivType dt, CONST V& a, const MM& b, std::string label)
{
  DoTestLDivVM2a<T,Tb>(dt,a,b,label);
#ifndef NOVIEWS
  DoTestLDivVM2a<T,Tb>(dt,a.Reverse(),b,label+" RevA");
  if (tmv::isComplex(T())) 
    DoTestLDivVM2a<T,Tb>(dt,a.Conjugate(),b,label+" ConjA");
#endif

#ifdef XTEST
  tmv::Vector<T> a0 = a;
  a.Zero();
  DoTestLDivVM2a<T,Tb>(dt,a,b,label+" 1");
  a = a0;

#ifndef NOVIEWS
  a.SubVector(0,a.size()/2).Zero();
  DoTestLDivVM2a<T,Tb>(dt,a,b,label+" 2");
  a = a0;

  a.SubVector(a.size()/2,a.size()).Zero();
  DoTestLDivVM2a<T,Tb>(dt,a,b,label+" 3");
  a = a0;

  a.SubVector(0,a.size()/4).Zero();
  a.SubVector(3*a.size()/4,a.size()).Zero();
  DoTestLDivVM2a<T,Tb>(dt,a,b,label+" 4");
  a = a0;
#endif
#endif
}
#endif

#ifndef NORDIVEQ
template <class T, class Tb, class V, class MM> static void DoTestRDivVM2a(
    tmv::DivType dt, CONST V& a, const MM& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start RDiv VM2a: "<<label<<std::endl;
  }

  tmv::Vector<T> v = a;
  tmv::Matrix<Tb> m = b;
  m.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  m.SaveDiv();

  double eps = EPS*Norm(m)*Norm(tmv::Matrix<Tb>(m.Inverse()))*v.size();

  if (CanRDivEqVM(a,b)) {
    tmv::Vector<T> a0 = a;
    tmv::Vector<T> frac = v%m;
    eps *= Norm(frac);
    a %= b;
    Assert(Norm(VEC(T,a)-frac) <= eps,label+" a%=b");
    a = a0;
    a *= b.Inverse();
    Assert(Norm(VEC(T,a)-frac) <= eps,label+" a*=b^-1");
    a = a0;
#ifdef ALIASOK
    a = a % b;
    Assert(Norm(VEC(T,a)-frac) <= eps,label+" a=a%b");
    a = a0;
    a = a * b.Inverse();
    Assert(Norm(VEC(T,a)-frac) <= eps,label+" a=a*b^-1");
    a = a0;
#ifdef XTEST
    TMV_RealType(T) x(5);
    T z;  SetZ(z);
    a = x * a % b;
    Assert(Norm(VEC(T,a)-x*frac) <= x*eps,label+" a=x*a%b");
    a = a0;
    a = x * a * b.Inverse();
    Assert(Norm(VEC(T,a)-x*frac) <= x*eps,label+" a=x*a*b^-1");
    a = a0;
    a = z * a % b;
    Assert(Norm(VEC(T,a)-z*frac) <= x*eps,label+" a=z*a%b");
    a = a0;
    a = z * a * b.Inverse();
    Assert(Norm(VEC(T,a)-z*frac) <= x*eps,label+" a=z*a*b^-1");
    a = a0;
#endif
#endif
  }
  if (showstartdone) 
    std::cout<<"Done RDiv VM2a"<<std::endl;
}

template <class T, class Tb, class V, class MM> static void DoTestRDivVM2(
    tmv::DivType dt, CONST V& a, const MM& b, std::string label)
{
  DoTestRDivVM2a<T,Tb>(dt,a,b,label);
#ifndef NOVIEWS
  DoTestRDivVM2a<T,Tb>(dt,a.Reverse(),b,label+" RevA");
  if (tmv::isComplex(T()))
    DoTestRDivVM2a<T,Tb>(dt,a.Conjugate(),b,label+" ConjA");
#endif

#ifdef XTEST
  tmv::Vector<T> a0 = a;
  a.Zero();
  DoTestRDivVM2a<T,Tb>(dt,a,b,label+" 1");
  a = a0;

#ifndef NOVIEWS
  a.SubVector(0,a.size()/2).Zero();
  DoTestRDivVM2a<T,Tb>(dt,a,b,label+" 2");
  a = a0;

  a.SubVector(a.size()/2,a.size()).Zero();
  DoTestRDivVM2a<T,Tb>(dt,a,b,label+" 3");
  a = a0;

  a.SubVector(0,a.size()/4).Zero();
  a.SubVector(3*a.size()/4,a.size()).Zero();
  DoTestRDivVM2a<T,Tb>(dt,a,b,label+" 4");
  a = a0;
#endif
#endif
}
#endif

template <class Ta, class Tb, class T, class V1, class MM, class V2>
static void DoTestLDivVM3a(
    tmv::DivType dt, const V1& a, const MM& b, CONST V2& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start LDiv VM3a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TMV_Text(c)<<" "<<c<<std::endl;
  }

  tmv::Vector<Ta> v1 = a;
  tmv::Matrix<Tb> m = b;
  tmv::Vector<T> v2 = c;
  m.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  m.SaveDiv();

  double eps = EPS*Norm(m)*Norm(tmv::Matrix<Tb>(m.Inverse()))*v1.size();

  if (XXDEBUG3) {
    std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a.step()<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TMV_Text(c)<<"  "<<c.step()<<"  "<<c<<std::endl;
  }

  if (CanLDivVM(a,b,c)) {
    v2 = v1/m;
    eps *= Norm(v2);
    c = a/b;
    if (XXDEBUG3) {
      std::cout<<"v/m = "<<v2<<std::endl;
      std::cout<<"a/b = "<<c<<std::endl;
    }
    Assert(Norm(VEC(T,c)-v2) <= eps,label+" c=a/b");
#ifdef XTEST
    c = -a/b;
    Assert(Norm(VEC(T,c)-(-v2)) <= eps,label+" c=-a/b");
    TMV_RealType(T) x(5);
    T z; SetZ(z);
    c = x*a/b;
    Assert(Norm(VEC(T,c)-(x*v2)) <= x*eps,label+" c=x*a/b");
    c = a/b*x;
    Assert(Norm(VEC(T,c)-(x*v2)) <= x*eps,label+" c=a/b*x");
    c = (x*a)/b;
    Assert(Norm(VEC(T,c)-(x*v2)) <= x*eps,label+" c=(x*a)/b");
    c = z*a/b;
    Assert(Norm(VEC(T,c)-(z*v2)) <= x*eps,label+" c=z*a/b");
    c = a/b*z;
    Assert(Norm(VEC(T,c)-(z*v2)) <= x*eps,label+" c=a/b*z");
    c = (z*a)/b;
    Assert(Norm(VEC(T,c)-(z*v2)) <= x*eps,label+" c=(z*a)/b");
#endif
  }

  if (showstartdone)
    std::cout<<"Done LDiv VM3a"<<std::endl;
}

template <class Ta, class Tb, class T, class V1, class MM, class V2>
static void DoTestRDivVM3a(
    tmv::DivType dt, const V1& a, const MM& b, CONST V2& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start RDiv VM3a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TMV_Text(c)<<" "<<c<<std::endl;
  }

  tmv::Vector<Ta> v1 = a;
  tmv::Matrix<Tb> m = b;
  tmv::Vector<T> v2 = c;
  m.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  m.SaveDiv();

  double eps = EPS*Norm(m)*Norm(tmv::Matrix<Tb>(m.Inverse()))*v1.size();

  if (XXDEBUG4) {
    std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a.step()<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TMV_Text(c)<<"  "<<c.step()<<"  "<<c<<std::endl;
  }

  if (CanRDivVM(a,b,c)) {
    v2 = v1%m;
    eps *= Norm(v2);
    c = a%b;
    if (XXDEBUG4) {
      std::cout<<"v%m = "<<v2<<std::endl;
      std::cout<<"a%b = "<<c<<std::endl;
    }
    Assert(Norm(VEC(T,c)-v2) <= eps,label+" c=a%b");
#ifdef XTEST
    c = -a%b;
    Assert(Norm(VEC(T,c)-(-v2)) <= eps,label+" c=-a%b");
    TMV_RealType(T) x(5);
    T z; SetZ(z);
    c = x*a%b;
    Assert(Norm(VEC(T,c)-(x*v2)) <= x*eps,label+" c=x*a%b");
    c = a%b*x;
    Assert(Norm(VEC(T,c)-(x*v2)) <= x*eps,label+" c=x*a%b");
    c = (x*a)%b;
    Assert(Norm(VEC(T,c)-(x*v2)) <= x*eps,label+" c=(x*a)%b");
    c = z*a%b;
    Assert(Norm(VEC(T,c)-(z*v2)) <= x*eps,label+" c=z*a%b");
    c = a%b*z;
    Assert(Norm(VEC(T,c)-(z*v2)) <= x*eps,label+" c=a%b*z");
    c = (z*a)%b;
    Assert(Norm(VEC(T,c)-(z*v2)) <= x*eps,label+" c=(z*a)%b");
#endif
  }

  if (showstartdone)
    std::cout<<"Done RDiv VM3a"<<std::endl;
}

template <class Ta, class Tb, class T, class V1, class MM, class V2>
static void DoTestLDivVM3(
    tmv::DivType dt, const V1& a, const MM& b, CONST V2& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start LDiv VM3"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TMV_Text(c)<<" "<<c<<std::endl;
  }

  DoTestLDivVM3a<Ta,Tb,T>(dt,a,b,c,label);
#ifndef NOVIEWS
  DoTestLDivVM3a<Ta,Tb,T>(dt,a.Reverse(),b,c,label+" RevA");
  DoTestLDivVM3a<Ta,Tb,T>(dt,a,b,c.Reverse(),label+" RevC");
  DoTestLDivVM3a<Ta,Tb,T>(dt,a.Reverse(),b,c.Reverse(),label+" RevAC");

#ifdef XTEST
  if (tmv::isComplex(Ta())) 
    DoTestLDivVM3a<Ta,Tb,T>(dt,Conjugate(a),b,c,label+" ConjA");
  if (tmv::isComplex(T())) 
    DoTestLDivVM3a<Ta,Tb,T>(dt,a,b,Conjugate(c),label+" ConjC");
  if (tmv::isComplex(Ta()) && tmv::isComplex(T())) 
    DoTestLDivVM3a<Ta,Tb,T>(dt,Conjugate(a),b,Conjugate(c),label+" ConjAC");
#endif
#endif

  if (showstartdone)
    std::cout<<"Done LDiv VM3"<<std::endl;
}

template <class Ta, class Tb, class T, class V1, class MM, class V2>
static void DoTestRDivVM3(
    tmv::DivType dt, const V1& a, const MM& b, CONST V2& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start RDiv VM3"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TMV_Text(c)<<" "<<c<<std::endl;
  }

  DoTestRDivVM3a<Ta,Tb,T>(dt,a,b,c,label);
#ifndef NOVIEWS
  DoTestRDivVM3a<Ta,Tb,T>(dt,a.Reverse(),b,c,label+" RevA");
  DoTestRDivVM3a<Ta,Tb,T>(dt,a,b,c.Reverse(),label+" RevC");
  DoTestRDivVM3a<Ta,Tb,T>(dt,a.Reverse(),b,c.Reverse(),label+" RevAC");

#ifdef XTEST
  if (tmv::isComplex(Ta())) 
    DoTestRDivVM3a<Ta,Tb,T>(dt,Conjugate(a),b,c,label+" ConjA");
  if (tmv::isComplex(T())) 
    DoTestRDivVM3a<Ta,Tb,T>(dt,a,b,Conjugate(c),label+" ConjC");
  if (tmv::isComplex(Ta()) && tmv::isComplex(T())) 
    DoTestRDivVM3a<Ta,Tb,T>(dt,Conjugate(a),b,Conjugate(c),label+" ConjAC");
#endif
#endif

  if (showstartdone)
    std::cout<<"Done RDiv VM3"<<std::endl;
}

template <class T, class Tb, IFTEMP1(class M0) IFTEMP1(class CM0) class M1, class M2> static void DoTestLDivMM1a(
    tmv::DivType dt, IFTEMP1(M0& temp) IFTEMP1(CM0& ctemp) 
    const M1& a, const M2& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start LDiv MM1b: "<<label<<std::endl;
  }

  tmv::Matrix<T> m1 = a;
  tmv::Matrix<Tb> m2 = b;
  m2.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  m2.SaveDiv();

  double eps = EPS*Norm(m2)*Norm(tmv::Matrix<Tb>(m2.Inverse()))*m1.colsize();

  if (CanLDiv(a,b)) {
    tmv::Matrix<ProductType(T,Tb)> frac = m1/m2;
    eps *= Norm(frac);
    if (XXDEBUG5) {
      std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
      std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
      std::cout<<"v = "<<tmv::TMV_Text(m1)<<"  "<<m1<<std::endl;
      std::cout<<"m = "<<tmv::TMV_Text(m2)<<"  "<<m2<<std::endl;
      std::cout<<"a/b = "<<(IFTEMP(temp=)a/b)<<std::endl;
      std::cout<<"m1/m2 = "<<frac<<std::endl;
      std::cout<<"Norm(diff) = "<<Norm(MAT2(T,Tb,IFTEMP(temp=)a/b)-frac)<<std::endl;
      std::cout<<"eps = "<<eps<<std::endl;
    }
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)a/b)-frac) <= eps,label+" a/b");
    TMV_RealType(T) x(5);
    TMV_ComplexType(T) z(3,4);
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)x*a/b)-x*frac) <= x*eps,label+" x*a/b");
    Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)z*a/b)-z*frac) <= x*eps,label+" z*a/b");
#ifdef XTEST
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)b.Inverse()*a)-frac) <= eps,label+" b^-1*a");
#ifndef NOMIX
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)a/m2)-frac) <= eps,label+" a/m2");
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)m2.Inverse()*a)-frac) <= eps,label+" m2^-1*a");
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)m1/b)-frac) <= eps,label+" m1/b");
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)b.Inverse()*m1)-frac) <= eps,label+" b^-1*m1");
#endif

    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)(x*a)/b)-x*frac) <= x*eps,label+" (x*a)/b");
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)x*(a/b))-x*frac) <= x*eps,label+" x*(a/b)");
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)a/(x*b))-frac/x) <= eps/x,label+" a/(x*b)");

    Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)(z*a)/b)-z*frac) <= x*eps,label+" (z*a)/b");
    Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)z*(a/b))-z*frac) <= x*eps,label+" z*(a/b)");
    Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)a/(z*b))-frac/z) <= eps/x,label+" a/(z*b)");

    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)(x*a)/(x*b))-frac) <= eps,label+" (x*a)/(x*b)");
    Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)(z*a)/(x*b))-(z/x)*frac) <= eps,
        label+" (z*a)/(x*b)");
    Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)(x*a)/(z*b))-(x/z)*frac) <= eps,
        label+" (x*a)/(z*b)");
    Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)(z*a)/(z*b))-frac) <= eps,
        label+" (z*a)/(z*b)");
#endif
  }

  if (showstartdone) 
    std::cout<<"Done LDiv MM1b"<<std::endl;
}

template <class T, class Tb, IFTEMP1(class M0) IFTEMP1(class CM0) class M1, class M2> static void DoTestRDivMM1a(
    tmv::DivType dt, IFTEMP1(M0& temp) IFTEMP1(CM0& ctemp) const M1& a, const M2& b,
    std::string label)
{
  if (showstartdone) {
    std::cout<<"Start RDiv MM1b: "<<label<<std::endl;
  }

  tmv::Matrix<T> m1 = a;
  tmv::Matrix<Tb> m2 = b;
  m2.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  m2.SaveDiv();

  double eps = EPS*Norm(m2)*Norm(tmv::Matrix<Tb>(m2.Inverse()))*m1.colsize();

  if (CanRDiv(a,b)) {
    tmv::Matrix<ProductType(T,Tb)> frac = m1%m2;
    eps *= Norm(frac);
    if (XXDEBUG6) {
      std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
      std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
      std::cout<<"m1 = "<<tmv::TMV_Text(m1)<<"  "<<m1<<std::endl;
      std::cout<<"m2 = "<<tmv::TMV_Text(m2)<<"  "<<m2<<std::endl;
      std::cout<<"a%b = "<<(IFTEMP(temp=)a%b)<<std::endl;
      std::cout<<"m1%m2 = "<<frac<<std::endl;
      std::cout<<"Norm(diff) = "<<Norm(MAT2(T,Tb,IFTEMP(temp=)a%b)-frac)<<std::endl;
      std::cout<<"eps = "<<eps<<std::endl;
    }
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)a%b)-frac) <= eps,label+" a%b");
    TMV_RealType(T) x(5);
    TMV_ComplexType(T) z(3,4);
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)x*a%b)-x*frac) <= x*eps,label+" x*a%b");
    Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)z*a%b)-z*frac) <= x*eps,label+" z*a%b");
#ifdef XTEST
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)a*b.Inverse())-frac) <= eps,label+" a*b^-1");
#ifndef NOMIX
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)a%m2)-frac) <= eps,label+" a%m2");
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)a*m2.Inverse())-frac) <= eps,label+" a*m2^-1");
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)m1%b)-frac) <= eps,label+" m1%b");
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)m1*b.Inverse())-frac) <= eps,label+" m1*b^-1");
#endif

    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)(x*a)%b)-x*frac) <= x*eps,label+" (x*a)%b");
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)x*(a%b))-x*frac) <= x*eps,label+" x*(a%b)");
    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)a%(x*b))-frac/x) <= eps/x,label+" a%(x*b)");

    Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)(z*a)%b)-z*frac) <= x*eps,label+" (z*a)%b");
    Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)z*(a%b))-z*frac) <= x*eps,label+" z*(a%b)");
    Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)a%(z*b))-frac/z) <= eps/x,label+" a%(z*b)");

    Assert(Norm(MAT2(T,Tb,IFTEMP(temp=)(x*a)%(x*b))-frac) <= eps,label+" (x*a)%(x*b)");
    Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)(z*a)%(x*b))-(z/x)*frac) <= eps,
        label+" (z*a)%(x*b)");
    Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)(x*a)%(z*b))-(x/z)*frac) <= eps,
        label+" (x*a)%(z*b)");
    Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)(z*a)%(z*b))-frac) <= eps,
        label+" (z*a)%(z*b)");
#endif
  }

  if (showstartdone) 
    std::cout<<"Done RDiv MM1b"<<std::endl;
}

template <class T, class Tb, IFTEMP1(class M0) IFTEMP1(class CM0) class M1, class M2> static void DoTestLDivMM1(
    tmv::DivType dt, IFTEMP1(M0& v0) IFTEMP1(CM0& cv0) 
    const M1& a, const M2& b, std::string label)
{
  DoTestLDivMM1a<T,Tb>(dt,IFTEMP1(v0) IFTEMP1(cv0) a,b,label);

#ifndef NOVIEWS
  if (a.IsSquare())
    DoTestLDivMM1a<T,Tb>(dt,IFTEMP1(v0) IFTEMP1(cv0) a.Transpose(),b,
        label+" TranA");
  if (tmv::isComplex(T())) {
    DoTestLDivMM1a<T,Tb>(dt,IFTEMP1(v0) IFTEMP1(cv0) a.Conjugate(),b,
        label+" ConjA");
    if (a.IsSquare())
      DoTestLDivMM1a<T,Tb>(dt,IFTEMP1(v0) IFTEMP1(cv0) a.Adjoint(),b,
          label+" AdjA");
  }
#endif
}

template <class T, class Tb, IFTEMP1(class M0)  IFTEMP1(class CM0)  class M1, class M2> static void DoTestRDivMM1(
    tmv::DivType dt, IFTEMP1(M0& v0)  IFTEMP1(CM0& cv0)  const M1& a, const M2& b,
    std::string label)
{
  DoTestRDivMM1a<T,Tb>(dt,IFTEMP1(v0) IFTEMP1(cv0) a,b,label);

#ifndef NOVIEWS
  if (a.IsSquare())
    DoTestRDivMM1a<T,Tb>(dt,IFTEMP1(v0) IFTEMP1(cv0) a.Transpose(),b,
        label+" TranA");
  if (tmv::isComplex(T())) {
    DoTestRDivMM1a<T,Tb>(dt,IFTEMP1(v0) IFTEMP1(cv0) a.Conjugate(),b,
        label+" ConjA");
    if (a.IsSquare())
      DoTestRDivMM1a<T,Tb>(dt,IFTEMP1(v0) IFTEMP1(cv0) a.Adjoint(),b,
          label+" AdjA");
  }
#endif
}

#ifndef NOLDIVEQ
template <class T, class Tb, class M0, class M1, class M2> static void DoTestLDivMM2a(
    tmv::DivType dt, M0& a0, CONST M1& a, const M2& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start LDiv MM2b: "<<label<<std::endl;
    std::cout<<"a0 = "<<a0<<std::endl;
    std::cout<<"a = "<<a<<std::endl;
    std::cout<<"b = "<<b<<std::endl;
  }

  tmv::Matrix<T> v = a;
  tmv::Matrix<Tb> m = b;
  m.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  m.SaveDiv();

  double eps = EPS*Norm(m)*Norm(tmv::Matrix<Tb>(m.Inverse()))*v.colsize();

  if (CanLDivEq(a,b)) {
    a0 = a;
    tmv::Matrix<T> frac = v/m;
    eps *= Norm(frac);
    a /= b;
    Assert(Norm(MAT(T,a)-frac) <= eps,label+" a/=b");
    a = a0;
#ifdef ALIASOK
    a = a / b;
    Assert(Norm(MAT(T,a)-frac) <= eps,label+" a=a/b");
    a = a0;
    a = b.Inverse() * a;
    Assert(Norm(MAT(T,a)-frac) <= eps,label+" a=b^-1*a");
    a = a0;
#ifdef XTEST
    TMV_RealType(T) x(5);
    T z;  SetZ(z);
    a = x * a / b;
    Assert(Norm(MAT(T,a)-x*frac) <= x*eps,label+" a=x*a/b");
    a = a0;
    a = x * b.Inverse() * a;
    Assert(Norm(MAT(T,a)-x*frac) <= x*eps,label+" a=x*b^-1*a");
    a = a0;
    a = z * a / b;
    Assert(Norm(MAT(T,a)-z*frac) <= x*eps,label+" a=z*a/b");
    a = a0;
    a = z * b.Inverse() * a;
    Assert(Norm(MAT(T,a)-z*frac) <= x*eps,label+" a=z*b^-1*a");
    a = a0;
#endif
#endif
  }
  if (showstartdone) {
    std::cout<<"Done LDiv MM2a"<<std::endl;
    std::cout<<"a0 = "<<a0<<std::endl;
    std::cout<<"a = "<<a<<std::endl;
    std::cout<<"b = "<<b<<std::endl;
  }
}

template <class T, class Tb, class M0, class M1, class M2> static void DoTestLDivMM2(
    tmv::DivType dt, M0& a0, CONST M1& a, const M2& b, std::string label)
{
  DoTestLDivMM2a<T,Tb>(dt,a0,a,b,label);
#ifndef NOVIEWS
  if (tmv::isComplex(T())) {
    DoTestLDivMM2a<T,Tb>(dt,a0,a.Conjugate(),b,label+" ConjA");
  }
#endif
}
#endif

#ifndef NORDIVEQ
template <class T, class Tb, class M0, class M1, class M2> static void DoTestRDivMM2a(
    tmv::DivType dt, M0& a0, CONST M1& a, const M2& b, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start RDiv MM2a: "<<label<<std::endl;
    std::cout<<"a0 = "<<a0<<std::endl;
    std::cout<<"a = "<<a<<std::endl;
    std::cout<<"b = "<<b<<std::endl;
  }

  tmv::Matrix<T> v = a;
  tmv::Matrix<Tb> m = b;
  m.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  m.SaveDiv();

  double eps = EPS*Norm(m)*Norm(tmv::Matrix<Tb>(m.Inverse()));

  if (CanRDivEq(a,b)) {
    a0 = a;
    tmv::Matrix<T> frac = v%m;
    eps *= Norm(frac);
    a %= b;
    Assert(Norm(MAT(T,a)-frac) <= eps,label+" a%=b");
    a = a0;
    a *= b.Inverse();
    Assert(Norm(MAT(T,a)-frac) <= eps,label+" a*=b^-1");
    a = a0;
#ifdef ALIASOK
    a = a % b;
    Assert(Norm(MAT(T,a)-frac) <= eps,label+" a=a%b");
    a = a0;
    a = a * b.Inverse();
    Assert(Norm(MAT(T,a)-frac) <= eps,label+" a=a*b^-1");
    a = a0;
#ifdef XTEST
    TMV_RealType(T) x(5);
    T z;  SetZ(z);
    a = x * a % b;
    Assert(Norm(MAT(T,a)-x*frac) <= x*eps,label+" a=x*a%b");
    a = a0;
    a = x * a * b.Inverse();
    Assert(Norm(MAT(T,a)-x*frac) <= x*eps,label+" a=x*a*b^-1");
    a = a0;
    a = z * a % b;
    Assert(Norm(MAT(T,a)-z*frac) <= x*eps,label+" a=z*a%b");
    a = a0;
    a = z * a * b.Inverse();
    Assert(Norm(MAT(T,a)-z*frac) <= x*eps,label+" a=z*a*b^-1");
    a = a0;
#endif
#endif
  }
  if (showstartdone) {
    std::cout<<"Done RDiv MM2a"<<std::endl;
    std::cout<<"a0 = "<<a0<<std::endl;
    std::cout<<"a = "<<a<<std::endl;
    std::cout<<"b = "<<b<<std::endl;
  }
}

template <class T, class Tb, class M0, class M1, class M2> static void DoTestRDivMM2(
    tmv::DivType dt, M0& a0, CONST M1& a, const M2& b, std::string label)
{
  DoTestRDivMM2a<T,Tb>(dt,a0,a,b,label);
#ifndef NOVIEWS
  if (tmv::isComplex(T())) {
    DoTestRDivMM2a<T,Tb>(dt,a0,a.Conjugate(),b,label+" ConjA");
  }
#endif
}
#endif

template <class Ta, class Tb, class T, class M1, class M2, class M3>
static void DoTestLDivMM3a(
    tmv::DivType dt, const M1& a, const M2& b, CONST M3& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start LDiv MM3a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TMV_Text(c)<<" "<<c<<std::endl;
  }

  tmv::Matrix<Ta> v1 = a;
  tmv::Matrix<Tb> m = b;
  tmv::Matrix<T> v2 = c;
  m.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  m.SaveDiv();

  double eps = EPS*Norm(m)*Norm(tmv::Matrix<Tb>(m.Inverse()))*v1.colsize();

  if (XXDEBUG7) {
    std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TMV_Text(c)<<"  "<<c<<std::endl;
  }

  if (CanLDiv(a,b,c)) {
    v2 = v1/m;
    eps *= Norm(v2);
    c = a/b;
    if (XXDEBUG7) {
      std::cout<<"v/m = "<<v2<<std::endl;
      std::cout<<"a/b = "<<c<<std::endl;
    }
    Assert(Norm(MAT(T,c)-v2) <= eps,label+" c=a/b");
#ifdef XTEST
    c = -a/b;
    Assert(Norm(MAT(T,c)-(-v2)) <= eps,label+" c=-a/b");
    TMV_RealType(T) x(5);
    T z; SetZ(z);
    c = x*a/b;
    Assert(Norm(MAT(T,c)-(x*v2)) <= x*eps,label+" c=x*a/b");
    c = a/b*x;
    Assert(Norm(MAT(T,c)-(x*v2)) <= x*eps,label+" c=a/b*x");
    c = (x*a)/b;
    Assert(Norm(MAT(T,c)-(x*v2)) <= x*eps,label+" c=(x*a)/b");
    c = z*a/b;
    Assert(Norm(MAT(T,c)-(z*v2)) <= x*eps,label+" c=z*a/b");
    c = a/b*z;
    Assert(Norm(MAT(T,c)-(z*v2)) <= x*eps,label+" c=a/b*z");
    c = (z*a)/b;
    Assert(Norm(MAT(T,c)-(z*v2)) <= x*eps,label+" c=(z*a)/b");
#endif
  }

  if (showstartdone)
    std::cout<<"Done LDiv MM3a"<<std::endl;
}

template <class Ta, class Tb, class T, class M1, class M2, class M3>
static void DoTestRDivMM3a(
    tmv::DivType dt, const M1& a, const M2& b, CONST M3& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start RDiv MM3a"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TMV_Text(c)<<" "<<c<<std::endl;
  }

  tmv::Matrix<Ta> v1 = a;
  tmv::Matrix<Tb> m = b;
  tmv::Matrix<T> v2 = c;
  m.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  m.SaveDiv();

  double eps = EPS*Norm(m)*Norm(tmv::Matrix<Tb>(m.Inverse()));

  if (XXDEBUG8) {
    std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TMV_Text(c)<<"  "<<c<<std::endl;
  }

  if (CanRDiv(a,b,c)) {
    v2 = v1%m;
    eps *= Norm(v2);
    c = a%b;
    if (XXDEBUG8) {
      std::cout<<"v%m = "<<v2<<std::endl;
      std::cout<<"a%b = "<<c<<std::endl;
      std::cout<<"c-v2 = "<<c-v2<<std::endl;
      std::cout<<"Norm(c-v2) = "<<Norm(MAT(T,c)-v2)<<std::endl;
      std::cout<<"eps = "<<EPS<<" * "<<Norm(m)<<" * "<<
      Norm(tmv::Matrix<Tb>(m.Inverse()))<<" * "<<
      v1.colsize()<<" = "<<eps<<std::endl;
    }
    Assert(Norm(MAT(T,c)-v2) <= eps,label+" c=a%b");
#ifdef XTEST
    c = -a%b;
    Assert(Norm(MAT(T,c)-(-v2)) <= eps,label+" c=-a%b");
    TMV_RealType(T) x(5);
    T z; SetZ(z);
    c = x*a%b;
    Assert(Norm(MAT(T,c)-(x*v2)) <= x*eps,label+" c=x*a%b");
    c = a%b*x;
    Assert(Norm(MAT(T,c)-(x*v2)) <= x*eps,label+" c=x*a%b");
    c = (x*a)%b;
    Assert(Norm(MAT(T,c)-(x*v2)) <= x*eps,label+" c=(x*a)%b");
    c = z*a%b;
    Assert(Norm(MAT(T,c)-(z*v2)) <= x*eps,label+" c=z*a%b");
    c = a%b*z;
    Assert(Norm(MAT(T,c)-(z*v2)) <= x*eps,label+" c=a%b*z");
    c = (z*a)%b;
    Assert(Norm(MAT(T,c)-(z*v2)) <= x*eps,label+" c=(z*a)%b");
#endif
  }

  if (showstartdone)
    std::cout<<"Done RDiv MM3a"<<std::endl;
}

template <class Ta, class Tb, class T, class M1, class M2, class M3>
static void DoTestLDivMM3(
    tmv::DivType dt, const M1& a, const M2& b, CONST M3& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start LDiv MM3"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TMV_Text(c)<<" "<<c<<std::endl;
  }

  DoTestLDivMM3a<Ta,Tb,T>(dt,a,b,c,label);
#ifndef NOVIEWS
  if (a.IsSquare())
    DoTestLDivMM3a<Ta,Tb,T>(dt,Transpose(a),b,c,label+" TranA");
  if (c.IsSquare())
    DoTestLDivMM3a<Ta,Tb,T>(dt,a,b,Transpose(c),label+" TranC");
  if (a.IsSquare() && c.IsSquare())
    DoTestLDivMM3a<Ta,Tb,T>(dt,Transpose(a),b,Transpose(c),label+" TranAC");

#ifdef XTEST
  if (tmv::isComplex(Ta())) {
    DoTestLDivMM3a<Ta,Tb,T>(dt,Conjugate(a),b,c,label+" ConjA");
    if (c.IsSquare())
      DoTestLDivMM3a<Ta,Tb,T>(dt,Conjugate(a),b,Transpose(c),
          label+" ConjA TranC");
    if (a.IsSquare())
      DoTestLDivMM3a<Ta,Tb,T>(dt,Adjoint(a),b,c,label+" AdjA");
    if (a.IsSquare() && c.IsSquare())
      DoTestLDivMM3a<Ta,Tb,T>(dt,Adjoint(a),b,Transpose(c),
          label+" AdjA TranC");
    if (tmv::isComplex(T())) {
      DoTestLDivMM3a<Ta,Tb,T>(dt,Conjugate(a),b,Conjugate(c),
          label+" ConjA ConjC");
      if (c.IsSquare())
        DoTestLDivMM3a<Ta,Tb,T>(dt,Conjugate(a),b,Adjoint(c),
            label+" ConjA AdjC");
      if (a.IsSquare())
        DoTestLDivMM3a<Ta,Tb,T>(dt,Adjoint(a),b,Conjugate(c),
            label+" AdjA ConjC");
      if (a.IsSquare() && c.IsSquare())
        DoTestLDivMM3a<Ta,Tb,T>(dt,Adjoint(a),b,Adjoint(c),
            label+" AdjA AdjC");
    }
  } else if (tmv::isComplex(T())) {
    DoTestLDivMM3a<Ta,Tb,T>(dt,a,b,Conjugate(c),label+" ConjC");
    if (c.IsSquare())
      DoTestLDivMM3a<Ta,Tb,T>(dt,a,b,Adjoint(c),label+" AdjC");
    if (a.IsSquare())
      DoTestLDivMM3a<Ta,Tb,T>(dt,Transpose(a),b,Conjugate(c),
          label+" TranA ConjC");
    if (a.IsSquare() && c.IsSquare())
      DoTestLDivMM3a<Ta,Tb,T>(dt,Transpose(a),b,Adjoint(c),
          label+" TranA AdjC");
  }
#endif
#endif

  if (showstartdone)
    std::cout<<"Done LDiv MM3"<<std::endl;
}

template <class Ta, class Tb, class T, class M1, class M2, class M3>
static void DoTestRDivMM3(
    tmv::DivType dt, const M1& a, const M2& b, CONST M3& c, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start RDiv MM3"<<label<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<" "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<" "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TMV_Text(c)<<" "<<c<<std::endl;
  }

  DoTestRDivMM3a<Ta,Tb,T>(dt,a,b,c,label);
#ifndef NOVIEWS
  if (a.IsSquare())
    DoTestRDivMM3a<Ta,Tb,T>(dt,Transpose(a),b,c,label+" TranA");
  if (c.IsSquare())
    DoTestRDivMM3a<Ta,Tb,T>(dt,a,b,Transpose(c),label+" TranC");
  if (a.IsSquare() && c.IsSquare())
    DoTestRDivMM3a<Ta,Tb,T>(dt,Transpose(a),b,Transpose(c),label+" TranAC");

#ifdef XTEST
  if (tmv::isComplex(Ta())) {
    DoTestRDivMM3a<Ta,Tb,T>(dt,Conjugate(a),b,c,label+" ConjA");
    if (c.IsSquare())
      DoTestRDivMM3a<Ta,Tb,T>(dt,Conjugate(a),b,Transpose(c),
          label+" ConjA TranC");
    if (a.IsSquare())
      DoTestRDivMM3a<Ta,Tb,T>(dt,Adjoint(a),b,c,label+" AdjA");
    if (a.IsSquare() && c.IsSquare())
      DoTestRDivMM3a<Ta,Tb,T>(dt,Adjoint(a),b,Transpose(c),
          label+" AdjA TranC");
    if (tmv::isComplex(T())) {
      DoTestRDivMM3a<Ta,Tb,T>(dt,Conjugate(a),b,Conjugate(c),
          label+" ConjA ConjC");
      if (c.IsSquare())
        DoTestRDivMM3a<Ta,Tb,T>(dt,Conjugate(a),b,Adjoint(c),
            label+" ConjA AdjC");
      if (a.IsSquare())
        DoTestRDivMM3a<Ta,Tb,T>(dt,Adjoint(a),b,Conjugate(c),
            label+" AdjA ConjC");
      if (a.IsSquare() && c.IsSquare())
        DoTestRDivMM3a<Ta,Tb,T>(dt,Adjoint(a),b,Adjoint(c),
            label+" AdjA AdjC");
    }
  } else if (tmv::isComplex(T())) {
    DoTestRDivMM3a<Ta,Tb,T>(dt,a,b,Conjugate(c),label+" ConjC");
    if (c.IsSquare())
      DoTestRDivMM3a<Ta,Tb,T>(dt,a,b,Adjoint(c),label+" AdjC");
    if (a.IsSquare())
      DoTestRDivMM3a<Ta,Tb,T>(dt,Transpose(a),b,Conjugate(c),
          label+" TranA ConjC");
    if (a.IsSquare() && c.IsSquare())
      DoTestRDivMM3a<Ta,Tb,T>(dt,Transpose(a),b,Adjoint(c),
          label+" TranA AdjC");
  }
#endif
#endif

  if (showstartdone)
    std::cout<<"Done RDiv MM3"<<std::endl;
}

template <class T, IFTEMP1(class M0) IFTEMP1(class CM0) class MM> static void DoTestDivMX(
    tmv::DivType dt, IFTEMP1(M0& temp) IFTEMP1(CM0& ctemp) 
    const MM& a, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start Div MX: "<<label<<std::endl;
  }

  tmv::Matrix<T> m = a;
  m.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  m.SaveDiv();

  double eps = EPS*Norm(m)*Norm(m.Inverse())*std::max(m.colsize(),m.rowsize());

  TMV_RealType(T) x(5);
  TMV_ComplexType(T) z(3,4);
  tmv::Matrix<T> xfrac = x/m;
  tmv::Matrix<TMV_ComplexType(T)> zfrac = z/m;
  double normfrac = Norm(xfrac);
  double normm = Norm(m);

  if (XXDEBUG9) {
    std::cout<<"eps = "<<eps<<std::endl;
    std::cout<<"x = "<<x<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
    std::cout<<"m = "<<tmv::TMV_Text(m)<<"  "<<m<<std::endl;
    std::cout<<"x/a = "<<(IFTEMP(temp=)x/a)<<std::endl;
    std::cout<<"x/m = "<<xfrac<<std::endl;
    std::cout<<"a*(x/a) = "<<tmv::Matrix<T>(a*(IFTEMP(temp=)x/a))<<std::endl;
    std::cout<<"m*(x/m) = "<<m*xfrac<<std::endl;
    std::cout<<"Norm(diff) = "<<Norm(MAT(T,IFTEMP(temp=)x/a)-xfrac)<<std::endl;
    std::cout<<"eps*Norm(diff) = "<<eps*normfrac<<std::endl;
  }
  Assert(Norm(MAT(T,IFTEMP(temp=)x/a)-xfrac) <= x*eps*normfrac,label+" x/a");
  Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)z/a)-zfrac) <= x*eps*normfrac,label+" z/a");
#ifdef XTEST
  Assert(Norm(MAT(T,IFTEMP(temp=)a.Inverse()*x)-xfrac) <= x*eps*normfrac,
      label+" a^-1*x");
  Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)a.Inverse()*z)-zfrac) <= x*eps*normfrac,
      label+" a^-1*z");
#endif

  if (XXDEBUG9) {
    std::cout<<"x%a = "<<(IFTEMP(temp=)x%a)<<std::endl;
    std::cout<<"x%m = "<<xfrac<<std::endl;
    std::cout<<"x%a*a = "<<tmv::Matrix<T>((IFTEMP(temp=)x%a)*a)<<std::endl;
    std::cout<<"x%m*m = "<<xfrac*m<<std::endl;
    std::cout<<"Norm(diff) = "<<Norm(MAT(T,IFTEMP(temp=)x%a)-xfrac)<<std::endl;
  }
  Assert(Norm(MAT(T,IFTEMP(temp=)x%a)-xfrac) <= x*eps*normfrac,label+" x%a");
#ifdef XTEST
  Assert(Norm(MAT(T,IFTEMP(temp=)x*a.Inverse())-xfrac) <= x*eps*normfrac,
      label+" x*a^-1");
  Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)z%a)-zfrac) <= x*eps*normfrac,label+" z%a");
  Assert(Norm(MAT(TMV_ComplexType(T),IFTEMP(ctemp=)z*a.Inverse())-zfrac) <= x*eps*normfrac,
      label+" z*a^-1");
#endif

#ifdef USETEMP
  temp = a.Inverse();
#else
  tmv::Matrix<T> temp = a.Inverse();
#endif
  if (XXDEBUG9) {
    std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
    std::cout<<"1/a = "<<temp<<std::endl;
    std::cout<<"a*(1/a) = "<<tmv::Matrix<T>(a*temp)<<std::endl;
    std::cout<<"(1/a)*a = "<<tmv::Matrix<T>(temp*a)<<std::endl;
    std::cout<<"(1/a)*a-((1/a)*a)T = "<<
    tmv::Matrix<T>(temp*a)-Adjoint(tmv::Matrix<T>(temp*a))<<std::endl;
    std::cout<<"(a*(1/a))*a = "<<tmv::Matrix<T>(a*temp)*a<<std::endl;
    std::cout<<"(a*(1/a))*a-a = "<<tmv::Matrix<T>(a*temp)*a-a<<std::endl;
    std::cout<<Norm(tmv::Matrix<T>(a*temp)*a-a)<<"  "<<eps*normm<<std::endl;
    std::cout<<Norm(tmv::Matrix<T>(temp*a)-Adjoint(tmv::Matrix<T>(temp*a)))<<
    "  "<<eps<<std::endl;
  }
  Assert(Norm(tmv::Matrix<T>(a*temp)*a-MAT(T,a)) <= eps*normm,label+" a*(1/a)*a");
  Assert(Norm(tmv::Matrix<T>(temp*a)-Adjoint(tmv::Matrix<T>(temp*a))) <= eps,
      label+" (1/a)*a-((1/a)*a)T");

  if (showstartdone) 
    std::cout<<"Done MX"<<std::endl;
}

template <class T, class SM0, class CSM0, class SM1, class SM2, class CSM1, class CSM2> 
static void TestMatrixDivArith1(
    tmv::DivType dt, SM0& b0, CSM0& cb0, const SM1& a, const SM2& b,
    const CSM1& ca, const CSM2& cb, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start Test Div 1: "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
    std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
    std::cout<<"cb = "<<tmv::TMV_Text(cb)<<"  "<<cb<<std::endl;
    std::cout<<"b0 = "<<tmv::TMV_Text(b0)<<"  "<<b0<<std::endl;
    std::cout<<"cb0 = "<<tmv::TMV_Text(cb0)<<"  "<<cb0<<std::endl;
  }

  a.DivideUsing(dt);
  ca.DivideUsing(dt);
  a.SaveDiv();
  ca.SaveDiv();

  tmv::Matrix<CT> cc(a.rowsize(),b.rowsize());
  tmv::Matrix<CT> cd(b.colsize(),a.colsize());
#ifdef XTEST
  tmv::Matrix<T> c(a.rowsize(),b.rowsize());
  tmv::Matrix<T> d(b.colsize(),a.colsize());
  DoTestLDivMM1<T,T>(dt,IFTEMP1(c) IFTEMP1(cc) b,a,label+" R,R");
  DoTestRDivMM1<T,T>(dt,IFTEMP1(d) IFTEMP1(cd) b,a,label+" R,R");
  DoTestLDivMM1<T,CT>(dt,IFTEMP1(cc) IFTEMP1(cc) b,ca,label+" R,C");
  DoTestRDivMM1<T,CT>(dt,IFTEMP1(cd) IFTEMP1(cd) b,ca,label+" R,C");
  DoTestLDivMM1<CT,T>(dt,IFTEMP1(cc) IFTEMP1(cc) cb,a,label+" C,R");
  DoTestRDivMM1<CT,T>(dt,IFTEMP1(cd) IFTEMP1(cd) cb,a,label+" C,R");
#endif
  DoTestLDivMM1<CT,CT>(dt,IFTEMP1(cc) IFTEMP1(cc) cb,ca,label+" C,C");
  DoTestRDivMM1<CT,CT>(dt,IFTEMP1(cd) IFTEMP1(cd) cb,ca,label+" C,C");

#ifndef NOLDIVEQ
#ifdef XTEST
  DoTestLDivMM2<T,T>(dt,b0,b.View(),a,label+" R,R");
  DoTestLDivMM2<CT,T>(dt,cb0,cb.View(),a,label+" C,R");
#endif
  DoTestLDivMM2<CT,CT>(dt,cb0,cb.View(),ca,label+" C,C");
#endif

#ifndef NORDIVEQ
#ifdef XTEST
  DoTestRDivMM2<T,T>(dt,b0,b.View(),a,label+" R,R");
  DoTestRDivMM2<CT,T>(dt,cb0,cb.View(),a,label+" C,R");
#endif
  DoTestRDivMM2<CT,CT>(dt,cb0,cb.View(),ca,label+" C,C");
#endif

#ifdef XTEST
  DoTestLDivMM3<T,T,T>(dt,b,a,c.View(),label+" R,R,R");
  DoTestLDivMM3<T,T,CT>(dt,b,a,cc.View(),label+" R,R,C");
  DoTestRDivMM3<T,T,T>(dt,b,a,d.View(),label+" R,R,R");
  DoTestRDivMM3<T,T,CT>(dt,b,a,cd.View(),label+" R,R,C");
  DoTestLDivMM3<T,CT,CT>(dt,b,ca,cc.View(),label+" R,C,C");
  DoTestRDivMM3<T,CT,CT>(dt,b,ca,cd.View(),label+" R,C,C");
  DoTestLDivMM3<CT,T,CT>(dt,cb,a,cc.View(),label+" C,R,C");
  DoTestRDivMM3<CT,T,CT>(dt,cb,a,cd.View(),label+" C,R,C");
#endif
  DoTestLDivMM3<CT,CT,CT>(dt,cb,ca,cc.View(),label+" C,C,C");
  DoTestRDivMM3<CT,CT,CT>(dt,cb,ca,cd.View(),label+" C,C,C");
}

template <class T, class SM0, class CSM0, class SM1, class SM2, class CSM1, class CSM2> 
static void TestMatrixDivArith2(
    tmv::DivType dt, SM0& b0, CSM0& cb0,
    const SM1& a, const SM2& b,
    const CSM1& ca, const CSM2& cb, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start Test Div 2: "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
    std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
    std::cout<<"cb = "<<tmv::TMV_Text(cb)<<"  "<<cb<<std::endl;
    std::cout<<"b0 = "<<tmv::TMV_Text(b0)<<"  "<<b0<<std::endl;
    std::cout<<"cb0 = "<<tmv::TMV_Text(cb0)<<"  "<<cb0<<std::endl;
  }

  a.DivideUsing(dt);
  ca.DivideUsing(dt);
  a.SaveDiv();
  ca.SaveDiv();

#ifdef USETEMP
  tmv::Matrix<T> at = a.Transpose();
  tmv::Matrix<CT> cat = ca.Transpose();
#endif
#ifdef XTEST
  DoTestDivMX<T>(dt,IFTEMP1(at) IFTEMP1(cat) a,label+" R");
#endif
  DoTestDivMX<CT>(dt,IFTEMP1(cat) IFTEMP1(cat) ca,label+" C");

  tmv::Matrix<T> m2(b);
  tmv::Matrix<CT> cm2(cb);
  tmv::Vector<T> v(b.colsize());
  tmv::Vector<CT> cv(b.colsize());
  if (b.rowsize() > 0) {
    v = m2.col(0);
    cv = cm2.col(0);
  }
  else {
    v.SetAllTo(1);
    cv.SetAllTo(CT(1,2));
  }

  tmv::Vector<T> w(b.rowsize());
  tmv::Vector<CT> cw(b.rowsize());
  if (b.colsize() > 0) {
    w = m2.row(0);
    cw = cm2.row(0);
  }
  else {
    w.SetAllTo(1);
    cw.SetAllTo(CT(1,2));
  }
  tmv::Vector<CT> cx(a.rowsize());
  tmv::Vector<CT> cy(a.colsize());
#ifdef XTEST
  tmv::Vector<T> x(a.rowsize());
  tmv::Vector<T> y(a.colsize());
  DoTestLDivVM1<T,T>(dt,IFTEMP1(x) IFTEMP1(cx) v,a,label+" R,R");
  DoTestRDivVM1<T,T>(dt,IFTEMP1(y) IFTEMP1(cy) w,a,label+" R,R");
  DoTestLDivVM1<T,CT>(dt,IFTEMP1(cx) IFTEMP1(cx) v,ca,label+" R,C");
  DoTestRDivVM1<T,CT>(dt,IFTEMP1(cy) IFTEMP1(cy) w,ca,label+" R,C");
  DoTestLDivVM1<CT,T>(dt,IFTEMP1(cx) IFTEMP1(cx) cv,a,label+" C,R");
  DoTestRDivVM1<CT,T>(dt,IFTEMP1(cy) IFTEMP1(cy) cw,a,label+" C,R");
#endif
  DoTestLDivVM1<CT,CT>(dt,IFTEMP1(cx) IFTEMP1(cx) cv,ca,label+" C,C");
  DoTestRDivVM1<CT,CT>(dt,IFTEMP1(cy) IFTEMP1(cy) cw,ca,label+" C,C");

#ifndef NOLDIVEQ
#ifdef XTEST
  DoTestLDivVM2<T,T>(dt,v.View(),a,label+" R,R");
  DoTestLDivVM2<CT,T>(dt,cv.View(),a,label+" C,R");
#endif
  DoTestLDivVM2<CT,CT>(dt,cv.View(),ca,label+" C,C");
#endif

#ifndef NORDIVEQ
#ifdef XTEST
  DoTestRDivVM2<T,T>(dt,w.View(),a,label+" R,R");
  DoTestRDivVM2<CT,T>(dt,cw.View(),a,label+" C,R");
#endif
  DoTestRDivVM2<CT,CT>(dt,cw.View(),ca,label+" C,C");
#endif

#ifdef XTEST
  DoTestLDivVM3<T,T,T>(dt,v,a,x.View(),label+" R,R,R");
  DoTestLDivVM3<T,T,CT>(dt,v,a,cx.View(),label+" R,R,C");
  DoTestRDivVM3<T,T,T>(dt,w,a,y.View(),label+" R,R,R");
  DoTestRDivVM3<T,T,CT>(dt,w,a,cy.View(),label+" R,R,C");
  DoTestLDivVM3<T,CT,CT>(dt,v,ca,cx.View(),label+" R,C,C");
  DoTestRDivVM3<T,CT,CT>(dt,w,ca,cy.View(),label+" R,C,C");
  DoTestLDivVM3<CT,T,CT>(dt,cv,a,cx.View(),label+" C,R,C");
  DoTestRDivVM3<CT,T,CT>(dt,cw,a,cy.View(),label+" C,R,C");
#endif
  DoTestLDivVM3<CT,CT,CT>(dt,cv,ca,cx.View(),label+" C,C,C");
  DoTestRDivVM3<CT,CT,CT>(dt,cw,ca,cy.View(),label+" C,C,C");

  tmv::Matrix<CT> cc(a.rowsize(),b.rowsize());
  tmv::Matrix<CT> cd(b.colsize(),a.colsize());
#ifdef XTEST
  tmv::Matrix<T> c(a.rowsize(),b.rowsize());
  tmv::Matrix<T> d(b.colsize(),a.colsize());
  DoTestLDivMM1<T,T>(dt,IFTEMP1(c) IFTEMP1(cc) b,a,label+" R,R");
  DoTestRDivMM1<T,T>(dt,IFTEMP1(d) IFTEMP1(cd) b,a,label+" R,R");
  DoTestLDivMM1<T,CT>(dt,IFTEMP1(cc) IFTEMP1(cc) b,ca,label+" R,C");
  DoTestRDivMM1<T,CT>(dt,IFTEMP1(cd) IFTEMP1(cd) b,ca,label+" R,C");
  DoTestLDivMM1<CT,T>(dt,IFTEMP1(cc) IFTEMP1(cc) cb,a,label+" C,R");
  DoTestRDivMM1<CT,T>(dt,IFTEMP1(cd) IFTEMP1(cd) cb,a,label+" C,R");
#endif
  DoTestLDivMM1<CT,CT>(dt,IFTEMP1(cc) IFTEMP1(cc) cb,ca,label+" C,C");
  DoTestRDivMM1<CT,CT>(dt,IFTEMP1(cd) IFTEMP1(cd) cb,ca,label+" C,C");

#ifndef NOLDIVEQ
#ifdef XTEST
  DoTestLDivMM2<T,T>(dt,b0,b.View(),a,label+" R,R");
  DoTestLDivMM2<CT,T>(dt,cb0,cb.View(),a,label+" C,R");
#endif
  DoTestLDivMM2<CT,CT>(dt,cb0,cb.View(),ca,label+" C,C");
#endif

#ifndef NORDIVEQ
#ifdef XTEST
  DoTestRDivMM2<T,T>(dt,b0,b.View(),a,label+" R,R");
  DoTestRDivMM2<CT,T>(dt,cb0,cb.View(),a,label+" C,R");
#endif
  DoTestRDivMM2<CT,CT>(dt,cb0,cb.View(),ca,label+" C,C");
#endif

#ifdef XTEST
  DoTestLDivMM3<T,T,T>(dt,b,a,c.View(),label+" R,R,R");
  DoTestLDivMM3<T,T,CT>(dt,b,a,cc.View(),label+" R,R,C");
  DoTestRDivMM3<T,T,T>(dt,b,a,d.View(),label+" R,R,R");
  DoTestRDivMM3<T,T,CT>(dt,b,a,cd.View(),label+" R,R,C");
  DoTestLDivMM3<T,CT,CT>(dt,b,ca,cc.View(),label+" R,C,C");
  DoTestRDivMM3<T,CT,CT>(dt,b,ca,cd.View(),label+" R,C,C");
  DoTestLDivMM3<CT,T,CT>(dt,cb,a,cc.View(),label+" C,R,C");
  DoTestRDivMM3<CT,T,CT>(dt,cb,a,cd.View(),label+" C,R,C");
#endif
  DoTestLDivMM3<CT,CT,CT>(dt,cb,ca,cc.View(),label+" C,C,C");
  DoTestRDivMM3<CT,CT,CT>(dt,cb,ca,cd.View(),label+" C,C,C");
}

template <class T, IFTEMP1(class SM0) IFTEMP1(class CSM0) class SM1, class CSM1> 
static void TestMatrixDivArith3a(
    tmv::DivType dt, IFTEMP1(SM0& temp) IFTEMP1(CSM0& ctemp) 
    const SM1& a, const CSM1& ca, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start Test Div 3a: "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
    std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
  }

  DoTestDivMX<T>(dt,IFTEMP1(temp) IFTEMP1(ctemp) a,label+" R");
  DoTestDivMX<CT>(dt,IFTEMP1(ctemp) IFTEMP1(ctemp) ca,label+" C");
}

template <class T, class SM0, class CSM0, class SM1, class SM2, class SM3, class CSM1, class CSM2, class CSM3> 
static void TestMatrixDivArith3b(
    tmv::DivType dt, SM0& b0, CSM0& cb0,
    const SM1& a, SM2& b, SM3& c, 
    const CSM1& ca, CSM2& cb, CSM3& cc, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start Test Div 3b: "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TMV_Text(c)<<"  "<<c<<std::endl;
    std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
    std::cout<<"cb = "<<tmv::TMV_Text(cb)<<"  "<<cb<<std::endl;
    std::cout<<"cc = "<<tmv::TMV_Text(cc)<<"  "<<cc<<std::endl;
    std::cout<<"b0 = "<<tmv::TMV_Text(b0)<<"  "<<b0<<std::endl;
    std::cout<<"cb0 = "<<tmv::TMV_Text(cb0)<<"  "<<cb0<<std::endl;
  }

  DoTestLDivMM1<T,T>(dt,IFTEMP1(c) IFTEMP1(cc) b,a,label+" R,R");
  DoTestLDivMM1<CT,CT>(dt,IFTEMP1(cc) IFTEMP1(cc) cb,ca,label+" C,C");
#ifdef XTEST
  DoTestLDivMM1<T,CT>(dt,IFTEMP1(cc) IFTEMP1(cc) b,ca,label+" R,C");
  DoTestLDivMM1<CT,T>(dt,IFTEMP1(cc) IFTEMP1(cc) cb,a,label+" C,R");
#endif

#ifndef NOLDIVEQ
  DoTestLDivMM2<T,T>(dt,b0,b,a,label+" R,R");
  DoTestLDivMM2<CT,CT>(dt,cb0,cb,ca,label+" C,C");
#ifdef XTEST
  DoTestLDivMM2<CT,T>(dt,cb0,cb,a,label+" C,R");
#endif
#endif

  DoTestLDivMM3<T,T,T>(dt,b,a,c,label+" R,R,R");
  DoTestLDivMM3<CT,CT,CT>(dt,cb,ca,cc,label+" C,C,C");
#ifdef XTEST
  DoTestLDivMM3<T,T,CT>(dt,b,a,cc,label+" R,R,C");
  DoTestLDivMM3<T,CT,CT>(dt,b,ca,cc,label+" R,C,C");
  DoTestLDivMM3<CT,T,CT>(dt,cb,a,cc,label+" C,R,C");
#endif
}

template <class T, class SM0, class CSM0, class SM1, class SM2, class SM3, class CSM1, class CSM2, class CSM3> 
static void TestMatrixDivArith3c(
    tmv::DivType dt, SM0& b0, CSM0& cb0,
    const SM1& a, SM2& b, SM3& c, 
    const CSM1& ca, CSM2& cb, CSM3& cc, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start Test Div 3c: "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
    std::cout<<"c = "<<tmv::TMV_Text(c)<<"  "<<c<<std::endl;
    std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
    std::cout<<"cb = "<<tmv::TMV_Text(cb)<<"  "<<cb<<std::endl;
    std::cout<<"cc = "<<tmv::TMV_Text(cc)<<"  "<<cc<<std::endl;
    std::cout<<"b0 = "<<tmv::TMV_Text(b0)<<"  "<<b0<<std::endl;
    std::cout<<"cb0 = "<<tmv::TMV_Text(cb0)<<"  "<<cb0<<std::endl;
  }

  DoTestRDivMM1<T,T>(dt,IFTEMP1(c) IFTEMP1(cc) b,a,label+" R,R");
  DoTestRDivMM1<CT,CT>(dt,IFTEMP1(cc) IFTEMP1(cc) cb,ca,label+" C,C");
#ifdef XTEST
  DoTestRDivMM1<T,CT>(dt,IFTEMP1(cc) IFTEMP1(cc) b,ca,label+" R,C");
  DoTestRDivMM1<CT,T>(dt,IFTEMP1(cc) IFTEMP1(cc) cb,a,label+" C,R");
#endif

#ifndef NORDIVEQ
  DoTestRDivMM2<T,T>(dt,b0,b,a,label+" R,R");
  DoTestRDivMM2<CT,CT>(dt,cb0,cb,ca,label+" C,C");
#ifdef XTEST
  DoTestRDivMM2<CT,T>(dt,cb0,cb,a,label+" C,R");
#endif
#endif

  DoTestRDivMM3<T,T,T>(dt,b,a,c,label+" R,R,R");
  DoTestRDivMM3<CT,CT,CT>(dt,cb,ca,cc,label+" C,C,C");
#ifdef XTEST
  DoTestRDivMM3<T,T,CT>(dt,b,a,cc,label+" R,R,C");
  DoTestRDivMM3<T,CT,CT>(dt,b,ca,cc,label+" R,C,C");
  DoTestRDivMM3<CT,T,CT>(dt,cb,a,cc,label+" C,R,C");
#endif
}

template <class T, class SM1, class V1, class V2, class CSM1, class CV1, class CV2> 
static void TestMatrixDivArith3d(
    tmv::DivType dt, const SM1& a, const V1& v, V2& x,
    const CSM1& ca, const CV1& cv, CV2& cx, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start Test Div 3d: "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
    std::cout<<"v = "<<tmv::TMV_Text(v)<<"  "<<v<<std::endl;
    std::cout<<"x = "<<tmv::TMV_Text(x)<<"  "<<x<<std::endl;
    std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
    std::cout<<"cv = "<<tmv::TMV_Text(cv)<<"  "<<cv<<std::endl;
    std::cout<<"cx = "<<tmv::TMV_Text(cx)<<"  "<<cx<<std::endl;
  }

  DoTestLDivVM1<T,T>(dt,IFTEMP1(x) IFTEMP1(cx) v,a,label+" R,R");
  DoTestLDivVM1<CT,CT>(dt,IFTEMP1(cx) IFTEMP1(cx) cv,ca,label+" C,C");
#ifdef XTEST
  DoTestLDivVM1<T,CT>(dt,IFTEMP1(cx) IFTEMP1(cx) v,ca,label+" R,C");
  DoTestLDivVM1<CT,T>(dt,IFTEMP1(cx) IFTEMP1(cx) cv,a,label+" C,R");
#endif

#ifndef NOLDIVEQ
  DoTestLDivVM2<T,T>(dt,x,a,label+" R,R");
  DoTestLDivVM2<CT,CT>(dt,cx,ca,label+" C,C");
#ifdef XTEST
  DoTestLDivVM2<CT,T>(dt,cx,a,label+" C,R");
#endif
#endif

  DoTestLDivVM3<T,T,T>(dt,v,a,x,label+" R,R,R");
  DoTestLDivVM3<CT,CT,CT>(dt,cv,ca,cx,label+" C,C,C");
#ifdef XTEST
  DoTestLDivVM3<T,T,CT>(dt,v,a,cx,label+" R,R,C");
  DoTestLDivVM3<T,CT,CT>(dt,v,ca,cx,label+" R,C,C");
  DoTestLDivVM3<CT,T,CT>(dt,cv,a,cx,label+" C,R,C");
#endif
}

template <class T, class SM1, class V1, class V2, class CSM1, class CV1, class CV2> 
static void TestMatrixDivArith3e(
    tmv::DivType dt, const SM1& a, const V1& w, V2& y,
    const CSM1& ca, const CV1& cw, CV2& cy, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start Test Div 3e: "<<label<<std::endl;
    std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
    std::cout<<"w = "<<tmv::TMV_Text(w)<<"  "<<w<<std::endl;
    std::cout<<"y = "<<tmv::TMV_Text(y)<<"  "<<y<<std::endl;
    std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
    std::cout<<"cw = "<<tmv::TMV_Text(cw)<<"  "<<cw<<std::endl;
    std::cout<<"cy = "<<tmv::TMV_Text(cy)<<"  "<<cy<<std::endl;
  }

  DoTestRDivVM1<T,T>(dt,IFTEMP1(y) IFTEMP1(cy) w,a,label+" R,R");
  DoTestRDivVM1<CT,CT>(dt,IFTEMP1(cy) IFTEMP1(cy) cw,ca,label+" C,C");
#ifdef XTEST
  DoTestRDivVM1<T,CT>(dt,IFTEMP1(cy) IFTEMP1(cy) w,ca,label+" R,C");
  DoTestRDivVM1<CT,T>(dt,IFTEMP1(cy) IFTEMP1(cy) cw,a,label+" C,R");
#endif

#ifndef NORDIVEQ
  DoTestRDivVM2<T,T>(dt,y,a,label+" R,R");
  DoTestRDivVM2<CT,CT>(dt,cy,ca,label+" C,C");
#ifdef XTEST
  DoTestRDivVM2<CT,T>(dt,cy,a,label+" C,R");
#endif
#endif

  DoTestRDivVM3<T,T,T>(dt,w,a,y,label+" R,R,R");
  DoTestRDivVM3<CT,CT,CT>(dt,cw,ca,cy,label+" C,C,C");
#ifdef XTEST
  DoTestRDivVM3<T,T,CT>(dt,w,a,cy,label+" R,R,C");
  DoTestRDivVM3<T,CT,CT>(dt,w,ca,cy,label+" R,C,C");
  DoTestRDivVM3<CT,T,CT>(dt,cw,a,cy,label+" C,R,C");
#endif
}
