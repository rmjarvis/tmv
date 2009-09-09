#define CT std::complex<T>

template <class SM1, class SM2> inline bool CanLDiv(
    const SM1& a, const SM2& b)
{ return a.colsize() == b.colsize(); }

template <class T, class SM2> inline bool CanLDiv(
    const tmv::Vector<T>& a, const SM2& b)
{ return a.size() == b.colsize(); }

template <class SM1, class SM2> inline bool CanLDivEq(
    const SM1& a, const SM2& b)
{ return CanLDiv(a,b) && b.IsSquare(); }

template <class SM1, class SM2> inline bool CanRDiv(
    const SM1& a, const SM2& b)
{ return a.rowsize() == b.rowsize(); }

template <class T, class SM2> inline bool CanRDiv(
    const tmv::Vector<T>& a, const SM2& b)
{ return a.size() == b.rowsize(); }

template <class SM1, class SM2> inline bool CanRDivEq(
    const SM1& a, const SM2& b)
{ return CanRDiv(a,b) && b.IsSquare(); }

template <class SM1, class SM2, class T, class T2> 
inline void DoTestMatrixDivArithMM1b(
    const SM1& a, const SM2& b, const tmv::Matrix<T>& m1,
    const tmv::Matrix<T2>& m2, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM1b: "<<label<<std::endl;
  }
  b.SaveDiv();
  m2.SaveDiv();
  b.SetDiv();
  m2.SetDiv();

  if (XXDEBUG1) {
    std::cout<<"a = "<<tmv::Type(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<"  "<<b<<std::endl;
    std::cout<<"m1 = "<<tmv::Type(m1)<<"  "<<m1<<std::endl;
    std::cout<<"m2 = "<<tmv::Type(m2)<<"  "<<m2<<std::endl;
  }
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a!=m1");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b!=m2");

  double eps = EPS*Norm(m2)*Norm(m2.Inverse())*std::max(m2.colsize(),m2.rowsize());

  if (CanLDiv(a,b)) {
    tmv::Matrix<ComplexType(T)> frac = m1/m2;
    double normfrac = Norm(frac);
    if (XXDEBUG1) {
      std::cout<<"LDiv:\n";
      std::cout<<"a = "<<tmv::Type(a)<<a<<std::endl;
      std::cout<<"m1 = "<<tmv::Type(m1)<<m1<<std::endl;
      std::cout<<"b = "<<tmv::Type(b)<<b<<std::endl;
      std::cout<<"m2 = "<<tmv::Type(m2)<<m2<<std::endl;
      std::cout<<"m1/m2 = "<<frac<<std::endl;
      std::cout<<"a/m2 = "<<a/m2<<std::endl;
      std::cout<<"m1/b = "<<m1/b<<std::endl;
      std::cout<<"a/b = "<<a/b<<std::endl;
      std::cout<<"Norm(diff) = "<<Norm(a/b-frac)<<" <?  "<<eps*normfrac<<std::endl;
      std::cout<<"Norm(diff) = "<<Norm(a/m2-frac)<<" <?  "<<eps*normfrac<<std::endl;
      std::cout<<"Norm(diff) = "<<Norm(m1/b-frac)<<" <?  "<<eps*normfrac<<std::endl;
    }
    Assert(Norm((a/b)-frac) <= eps*normfrac,label+" a/b");
    Assert(Norm((a/m2)-frac) <= eps*normfrac,label+" a/m");
    Assert(Norm((m1/b)-frac) <= eps*normfrac,label+" m/b");
    Assert(Norm((b.Inverse()*a)-frac) <= eps*normfrac,label+" b^-1*a");
    Assert(Norm((m2.Inverse()*a)-frac) <= eps*normfrac,label+" m^-1*a");
    Assert(Norm((b.Inverse()*m1)-frac) <= eps*normfrac,label+" b^-1*m");
  }
  if (CanRDiv(a,b)) {
    tmv::Matrix<ComplexType(T)> frac = m1%m2;
    double normfrac = Norm(frac);
    if (XXDEBUG1) {
      std::cout<<"RDiv:\n";
      std::cout<<"a = "<<tmv::Type(a)<<a<<std::endl;
      std::cout<<"m1 = "<<tmv::Type(m1)<<m1<<std::endl;
      std::cout<<"b = "<<tmv::Type(b)<<b<<std::endl;
      std::cout<<"m2 = "<<tmv::Type(m2)<<m2<<std::endl;
      std::cout<<"m1%m2 = "<<m1%m2<<std::endl;
      std::cout<<"a%m2 = "<<a%m2<<std::endl;
      std::cout<<"m1%b = "<<m1%b<<std::endl;
      std::cout<<"a%b = "<<a%b<<std::endl;
      std::cout<<"Norm(diff) = "<<Norm(a%b-frac)<<" <?  "<<eps<<std::endl;
      std::cout<<"Norm(diff) = "<<Norm(a%m2-frac)<<" <?  "<<eps<<std::endl;
      std::cout<<"Norm(diff) = "<<Norm(m1%b-frac)<<" <?  "<<eps<<std::endl;
    }
    Assert(Norm((a%b)-frac) <= eps*normfrac,label+" a%b");
    Assert(Norm((a%m2)-frac) <= eps*normfrac,label+" a%m");
    Assert(Norm((m1%b)-frac) <= eps*normfrac,label+" m%b");
    Assert(Norm((a*b.Inverse())-frac) <= eps*normfrac,label+" a*b^-1");
    Assert(Norm((a*m2.Inverse())-frac) <= eps*normfrac,label+" a*m^-1");
    Assert(Norm((m1*b.Inverse())-frac) <= eps*normfrac,label+" m*b^-1");
  }

  if (showstartdone) 
    std::cout<<"Done MM1b"<<std::endl;
}

template <class SM1, class SM2, class T, class T2> 
inline void DoTestMatrixDivArithMM1a(
    tmv::DivType dt, const SM1& a, const SM2& b, const tmv::Matrix<T>& m1,
    const tmv::Matrix<T2>& m2, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM1a: "<<label<<std::endl;
  }
  b.DivideUsing(dt);
  m2.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  DoTestMatrixDivArithMM1b(a,b,m1,m2,label+" "+tmv::Text(dt));
  if (showstartdone) 
    std::cout<<"Done MM1a"<<std::endl;
}

template <class SM1, class SM2, class T, class T2> 
inline void DoTestMatrixDivArithMM1(
    tmv::DivType dt, const SM1& a, const SM2& b,
    const tmv::Matrix<T>& m1, const tmv::Matrix<T2>& m2, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM1: "<<label<<std::endl;
  }
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a!=m1");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b!=m2");

  tmv::Matrix<T> m1t = m1.Transpose();
  tmv::Matrix<T2> m2t = m2.Transpose();
  DoTestMatrixDivArithMM1a(dt,a.View(),b.View(),m1,m2,label);
  DoTestMatrixDivArithMM1a(dt,Transpose(a),b.View(),m1t,m2,label+" TransA");
  DoTestMatrixDivArithMM1a(dt,a.View(),Transpose(b),m1,m2t,label+" TransB");
  DoTestMatrixDivArithMM1a(dt,Transpose(a),Transpose(b),m1t,m2t,
      label+" TransA TransB");

  if (tmv::IsComplex(T())) {
    tmv::Matrix<T> m1c = m1.Conjugate();
    tmv::Matrix<T> m1a = m1.Adjoint();
    DoTestMatrixDivArithMM1a(dt,Conjugate(a),b.View(),m1c,m2,label+" ConjA");
    DoTestMatrixDivArithMM1a(dt,Adjoint(a),b.View(),m1a,m2,label+" AdjA");
    DoTestMatrixDivArithMM1a(dt,Conjugate(a),Transpose(b),m1c,m2t,
	label+" ConjA TransB");
    DoTestMatrixDivArithMM1a(dt,Adjoint(a),Transpose(b),m1a,m2t,
	label+" AdjA TransB");
  }

  if (tmv::IsComplex(T2())) {
    tmv::Matrix<T2> m2c = m2.Conjugate();
    tmv::Matrix<T2> m2a = m2.Adjoint();
    DoTestMatrixDivArithMM1a(dt,a.View(),Conjugate(b),m1,m2c,label+" ConjB");
    DoTestMatrixDivArithMM1a(dt,Transpose(a),Conjugate(b),m1t,m2c,
	label+" TransA ConjB");
    DoTestMatrixDivArithMM1a(dt,a.View(),Adjoint(b),m1,m2a,label+" AdjB");
    DoTestMatrixDivArithMM1a(dt,Transpose(a),Adjoint(b),m1t,m2a,
	label+" TransA AdjB");
  }

  if (tmv::IsComplex(T()) && tmv::IsComplex(T2())) {
    tmv::Matrix<T> m1c = m1.Conjugate();
    tmv::Matrix<T> m1a = m1.Adjoint();
    tmv::Matrix<T2> m2c = m2.Conjugate();
    tmv::Matrix<T2> m2a = m2.Adjoint();
    DoTestMatrixDivArithMM1a(dt,Conjugate(a),Conjugate(b),m1c,m2c,
	label+" ConjA ConjB");
    DoTestMatrixDivArithMM1a(dt,Adjoint(a),Conjugate(b),m1a,m2c,
	label+" AdjA ConjB");
    DoTestMatrixDivArithMM1a(dt,Conjugate(a),Adjoint(b),m1c,m2a,
	label+" ConjA AdjB");
    DoTestMatrixDivArithMM1a(dt,Adjoint(a),Adjoint(b),m1a,m2a,
	label+" AdjA AdjB");
  }

  if (showstartdone) 
    std::cout<<"Done MM1"<<std::endl;
}

template <class SM1, class SM2, class T, class T2> 
inline void DoTestMatrixDivArithMM2b(
    const SM1& a, const SM2& b, const tmv::Matrix<T>& m1,
    const tmv::Matrix<T2>& m2, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM2b: "<<label<<std::endl;
  }

  if (XXDEBUG2) {
    std::cout<<"a = "<<tmv::Type(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<"  "<<b<<std::endl;
    std::cout<<"m1 = "<<tmv::Type(m1)<<"  "<<m1<<std::endl;
    std::cout<<"m2 = "<<tmv::Type(m2)<<"  "<<m2<<std::endl;
  }

  m2.SaveDiv();
  b.SaveDiv();

  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a!=m1");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b!=m2");

  tmv::Matrix<T2,tmv::ColMajor> m2inv = Inverse(m2);
  double eps = EPS*std::max(m2.colsize(),m2.rowsize())*Norm(m2)*Norm(m2inv);

  if (CanLDivEq(m1,b)) {
    tmv::Matrix<T> m3 = m1;
    tmv::Matrix<T> m4 = m1;
    m3 /= b;
    m4 = tmv::Matrix<T>(m4/m2);
    if (XXDEBUG2) {
      std::cout<<"m3 = "<<m3<<std::endl;
      std::cout<<"m4 = "<<m4<<std::endl;
      std::cout<<"m3*b = "<<b*m3<<std::endl;
      std::cout<<"m4*m2 = "<<m2*m4<<std::endl;
      std::cout<<"Norm(m3-m4) = "<<Norm(m3-m4)<<std::endl;
      std::cout<<"eps*Norm(m4) = "<<eps*Norm(m4)<<std::endl;
    }
    Assert(Norm(m3-m4) <= eps*Norm(m4),label+" m/=b");
    m3 = m1;
    m3 = b.Inverse()*m3;
    Assert(Norm(m3-m4) <= eps*Norm(m4),label+" m=b^-1*m");
  }
  if (CanRDivEq(m1,b)) {
    tmv::Matrix<T> m3 = m1;
    tmv::Matrix<T> m4 = m1;
    m3 %= b;
    m4 = tmv::Matrix<T>(m4%m2);
    if (XXDEBUG2) {
      std::cout<<"m3 = "<<m3<<std::endl;
      std::cout<<"m4 = "<<m4<<std::endl;
      std::cout<<"m3*b = "<<m3*b<<std::endl;
      std::cout<<"m4*m2 = "<<m4*m2<<std::endl;
      std::cout<<"Norm(m3-m4) = "<<Norm(m3-m4)<<std::endl;
      std::cout<<"eps*Norm(m4) = "<<eps*Norm(m4)<<std::endl;
    }
    Assert(Norm(m3-m4) <= eps*Norm(m4),label+" m%=b");
    m3 = m1;
    m3 *= b.Inverse();
    Assert(Norm(m3-m4) <= eps*Norm(m4),label+" m*=b^-1");
    m3 = m1;
    m3 = m3 * b.Inverse();
    Assert(Norm(m3-m4) <= eps*Norm(m4),label+" m=m*b^-1");
  }

  if (showstartdone) 
    std::cout<<"Done MM2b"<<std::endl;
}

template <class SM1, class SM2, class T, class T2> 
inline void DoTestMatrixDivArithMM2a(
    tmv::DivType dt, const SM1& a, const SM2& b, const tmv::Matrix<T>& m1,
    const tmv::Matrix<T2>& m2, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MM2a: "<<label<<std::endl;
  }
  b.DivideUsing(dt);
  m2.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  DoTestMatrixDivArithMM2b(a,b,m1,m2,label+" "+tmv::Text(dt));
  if (showstartdone) 
    std::cout<<"Done MM2a"<<std::endl;
}

template <class SM1, class SM2, class T, class T2> 
inline void DoTestMatrixDivArithMM2(
    tmv::DivType dt, const SM1& a, const SM2& b,
    const tmv::Matrix<T>& m1,
    const tmv::Matrix<T2>& m2, std::string label)
{
  Assert(Norm(a-m1) <= EPS*Norm(m1),label+" a!=m1");
  Assert(Norm(b-m2) <= EPS*Norm(m2),label+" b!=m2");


  tmv::Matrix<T> m1t = m1.Transpose();
  tmv::Matrix<T2> m2t = m2.Transpose();
  DoTestMatrixDivArithMM2a(dt,a.View(),b.View(),m1,m2,label);
  DoTestMatrixDivArithMM2a(dt,Transpose(a),b.View(),m1t,m2,label+" TransA");
  DoTestMatrixDivArithMM2a(dt,a.View(),Transpose(b),m1,m2t,label+" TransB");
  DoTestMatrixDivArithMM2a(dt,Transpose(a),Transpose(b),m1t,m2t,
      label+" TransA TransB");

  if (tmv::IsComplex(T())) {
    tmv::Matrix<T> m1c = m1.Conjugate();
    tmv::Matrix<T> m1a = m1.Adjoint();
    DoTestMatrixDivArithMM2a(dt,Conjugate(a),b.View(),m1c,m2,label+" ConjA");
    DoTestMatrixDivArithMM2a(dt,Adjoint(a),b.View(),m1a,m2,label+" AdjA");
    DoTestMatrixDivArithMM2a(dt,Conjugate(a),Transpose(b),m1c,m2t,
	label+" ConjA TransB");
    DoTestMatrixDivArithMM2a(dt,Adjoint(a),Transpose(b),m1a,m2t,
	label+" AdjA TransB");
  }

  if (tmv::IsComplex(T2())) {
    tmv::Matrix<T2> m2c = m2.Conjugate();
    tmv::Matrix<T2> m2a = m2.Adjoint();
    DoTestMatrixDivArithMM2a(dt,a.View(),Conjugate(b),m1,m2c,label+" ConjB");
    DoTestMatrixDivArithMM2a(dt,Transpose(a),Conjugate(b),m1t,m2c,
	label+" TransA ConjB");
    DoTestMatrixDivArithMM2a(dt,a.View(),Adjoint(b),m1,m2a,label+" AdjB");
    DoTestMatrixDivArithMM2a(dt,Transpose(a),Adjoint(b),m1t,m2a,
	label+" TransA AdjB");
  }

  if (tmv::IsComplex(T()) && tmv::IsComplex(T2())) {
    tmv::Matrix<T> m1c = m1.Conjugate();
    tmv::Matrix<T> m1a = m1.Adjoint();
    tmv::Matrix<T2> m2c = m2.Conjugate();
    tmv::Matrix<T2> m2a = m2.Adjoint();
    DoTestMatrixDivArithMM2a(dt,Conjugate(a),Conjugate(b),m1c,m2c,
	label+" ConjA ConjB");
    DoTestMatrixDivArithMM2a(dt,Adjoint(a),Conjugate(b),m1a,m2c,
	label+" AdjA ConjB");
    DoTestMatrixDivArithMM2a(dt,Conjugate(a),Adjoint(b),m1c,m2a,
	label+" ConjA AdjB");
    DoTestMatrixDivArithMM2a(dt,Adjoint(a),Adjoint(b),m1a,m2a,
	label+" AdjA AdjB");
  }
}

template <class SM, class T, class T2> inline void DoTestMatrixDivArithMV1b(
    const tmv::Vector<T>& v, const SM& b,
    const tmv::Matrix<T2>& m, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV1b: "<<label<<std::endl;
  }

  m.SaveDiv();
  b.SaveDiv();

  Assert(Norm(b-m) <= EPS*Norm(m),label+" b!=m");

  double eps = EPS*Norm(m)*Norm(m.Inverse())*v.size();

  if (CanLDiv(v,b)) {
    tmv::Vector<ComplexType(T)> frac = v/m;
    double normfrac = Norm(frac);
    if (XXDEBUG3) {
      std::cout<<"v = "<<tmv::Type(v)<<"  step "<<v.step()<<"  "<<v<<std::endl;
      std::cout<<"b = "<<tmv::Type(b)<<"  "<<b<<std::endl;
      std::cout<<"m = "<<tmv::Type(m)<<"  "<<m<<std::endl;
      std::cout<<"v/b = "<<v/b<<std::endl;
      std::cout<<"v/m = "<<frac<<std::endl;
      std::cout<<"Norm(diff) = "<<Norm((v/b)-frac)<<std::endl;
      std::cout<<"eps = "<<eps*normfrac<<std::endl;
    }
    Assert(Norm((v/b)-frac) <= eps*normfrac,label+" v/m");
    Assert(Norm((b.Inverse()*v)-frac) <= eps*normfrac,label+" m^-1*v");
  }
  if (CanRDiv(v,b)) {
    tmv::Vector<ComplexType(T)> frac = v%m;
    double normfrac = Norm(frac);
    if (XXDEBUG4) {
      std::cout<<"v = "<<tmv::Type(v)<<"  step "<<v.step()<<"  "<<v<<std::endl;
      std::cout<<"b = "<<tmv::Type(b)<<"  "<<b<<std::endl;
      std::cout<<"m = "<<tmv::Type(m)<<"  "<<m<<std::endl;
      std::cout<<"v%b = "<<v%b<<std::endl;
      std::cout<<"v%m = "<<frac<<std::endl;
      std::cout<<"(v%b)*m = "<<(v%b)*m<<std::endl;
      std::cout<<"(v%b)*b = "<<(v%b)*b<<std::endl;
      std::cout<<"(v%m)*m = "<<frac*m<<std::endl;
      std::cout<<"(v%m)*b = "<<frac*b<<std::endl;
      std::cout<<"diff = "<<(v%b)-frac<<std::endl;
      std::cout<<"Norm(diff) = "<<Norm((v%b)-frac)<<std::endl;
      std::cout<<"eps = "<<eps*normfrac<<std::endl;
    }
    Assert(Norm((v%b)-frac) <= eps*normfrac,label+" v%m");
    Assert(Norm((v*b.Inverse())-frac) <= eps*normfrac,label+" v*m^-1");
  }

  if (showstartdone) 
    std::cout<<"Done MV1b"<<std::endl;
}

template <class SM, class T, class T2> inline void DoTestMatrixDivArithMV1a(
    tmv::DivType dt, const tmv::Vector<T>& v, const SM& b, 
    const tmv::Matrix<T2>& m, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV1a: "<<label<<std::endl;
  }
  b.DivideUsing(dt);
  m.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  DoTestMatrixDivArithMV1b(v,b,m,label+" "+tmv::Text(dt));
  if (showstartdone) 
    std::cout<<"Done MV1a"<<std::endl;
}

template <class SM, class T, class T2> inline void DoTestMatrixDivArithMV1(
    tmv::DivType dt, const tmv::Vector<T>& v, const SM& b,
    const tmv::Matrix<T2>& m, std::string label)
{
  Assert(Norm(b-m) <= EPS*Norm(m),label+" b!=m");

  tmv::Matrix<T2> mt = m.Transpose();
  DoTestMatrixDivArithMV1a(dt,v,b.View(),m,label);
  DoTestMatrixDivArithMV1a(dt,v,Transpose(b),mt,label+" Trans");

  if (tmv::IsComplex(T2())) {
  tmv::Matrix<T2> mc = m.Conjugate();
  tmv::Matrix<T2> ma = m.Adjoint();
    DoTestMatrixDivArithMV1a(dt,v,Conjugate(b),mc,label+" Conj");
    DoTestMatrixDivArithMV1a(dt,v,Adjoint(b),ma,label+" Adj");
  }
}

template <class SM, class T, class T2> inline void DoTestMatrixDivArithMV2b(
    const tmv::Vector<T>& v, const SM& b,
    const tmv::Matrix<T2>& m, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV2b: "<<label<<std::endl;
  }
  Assert(Norm(b-m) <= EPS*Norm(m),label+" b!=m");

  m.SaveDiv();
  b.SaveDiv();

  tmv::Matrix<T2,tmv::ColMajor> minv = Inverse(m);
  double eps = EPS*Norm(m)*Norm(minv)*v.size();

  if (CanLDivEq(v,m)) {
    tmv::Vector<T> v1 = v;
    tmv::Vector<T> v2 = v/m;
    v1 /= b;
    Assert(Norm(v1-v2) <= eps*Norm(v1),label+" v/=m");
  }
  if (CanRDivEq(v,m)) {
    tmv::Vector<T> v1 = v;
    tmv::Vector<T> v2 = v%m;
    v1 %= b;
    Assert(Norm(v1-v2) <= eps*Norm(v1),label+" v%=m");
    v1 = v;
    v1 *= b.Inverse();
    Assert(Norm(v1-v2) <= eps*Norm(v1),label+" v*=m^-1");
  }
  if (showstartdone) 
    std::cout<<"Done MV2b"<<std::endl;
}

template <class SM, class T, class T2> inline void DoTestMatrixDivArithMV2a(
    tmv::DivType dt, const tmv::Vector<T>& v, const SM& b,
    const tmv::Matrix<T2>& m, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MV2a: "<<label<<std::endl;
  }
  b.DivideUsing(dt);
  m.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  DoTestMatrixDivArithMV2b(v,b,m,label+" "+tmv::Text(dt));
  if (showstartdone) 
    std::cout<<"Done MV2a"<<std::endl;
}

template <class SM, class T, class T2> inline void DoTestMatrixDivArithMV2(
    tmv::DivType dt, const tmv::Vector<T>& v, const SM& b,
    const tmv::Matrix<T2>& m, std::string label)
{
  Assert(Norm(b-m) <= EPS*Norm(m),label+" b!=m");

  tmv::Matrix<T2> mt = m.Transpose();
  DoTestMatrixDivArithMV2a(dt,v,b.View(),m,label);
  DoTestMatrixDivArithMV2a(dt,v,Transpose(b),mt,label+" TransA");

  if (tmv::IsComplex(T2())) {
    tmv::Matrix<T2> mc = m.Conjugate();
    tmv::Matrix<T2> ma = m.Adjoint();
    DoTestMatrixDivArithMV2a(dt,v,Conjugate(b),mc,label+" ConjA");
    DoTestMatrixDivArithMV2a(dt,v,Adjoint(b),ma,label+" AdjA");
  }
}

template <class SM1, class T, class T2> inline void DoTestMatrixDivArithMXb(
    const SM1& a, const tmv::Matrix<T>& m, T2 x, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MXb: "<<label<<std::endl;
  }

  m.SaveDiv();
  a.SaveDiv();

  double eps = EPS*Norm(m)*Norm(m.Inverse())*std::max(m.colsize(),m.rowsize());
  tmv::Matrix<ComplexType(T)> frac = x/m;
  double normfrac = Norm(frac);
  double normm = Norm(m);

  if (XXDEBUG5) {
    std::cout<<"eps = "<<eps<<std::endl;
    std::cout<<"x = "<<x<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<"  "<<a<<std::endl;
    std::cout<<"m = "<<tmv::Type(m)<<"  "<<m<<std::endl;
    std::cout<<"x/a = "<<x/a<<std::endl;
    std::cout<<"x/m = "<<frac<<std::endl;
    std::cout<<"a*(x/a) = "<<a*(x/a)<<std::endl;
    std::cout<<"m*(x/m) = "<<m*frac<<std::endl;
    std::cout<<"Norm(diff) = "<<Norm((x/a)-frac)<<std::endl;
    std::cout<<"eps*Norm(diff) = "<<eps*normfrac<<std::endl;
  }
  Assert(Norm((x/a)-frac) <= eps*normfrac,label+" x/a");
  Assert(Norm((a.Inverse()*x)-frac) <= eps*normfrac,label+" a^-1*x");

  if (XXDEBUG5) {
    std::cout<<"x%a = "<<x%a<<std::endl;
    std::cout<<"x%m = "<<frac<<std::endl;
    std::cout<<"x%a*a = "<<(x%a)*a<<std::endl;
    std::cout<<"x%m*m = "<<frac*m<<std::endl;
    std::cout<<"Norm(diff) = "<<Norm((x%a)-frac)<<std::endl;
  }
  Assert(Norm((x%a)-frac) <= eps*normfrac,label+" x%a");
  Assert(Norm((x*a.Inverse())-frac) <= eps*normfrac,label+" x*a^-1");
  if (XXDEBUG5) {
    std::cout<<"a = "<<tmv::Type(a)<<"  "<<a<<std::endl;
    std::cout<<"1/a = "<<T2(1)/a<<std::endl;
    std::cout<<"a*(1/a) = "<<a*(T2(1)/a)<<std::endl;
    std::cout<<"(1/a)*a = "<<(T2(1)/a)*a<<std::endl;
    std::cout<<"(1/a)*a-((1/a)*a)T = "<<(T2(1)/a)*a-Adjoint((T2(1)/a)*a)<<std::endl;
    std::cout<<"(a*(1/a))*a = "<<(a*(T2(1)/a))*a<<std::endl;
    std::cout<<"(a*(1/a))*a-a = "<<(a*(T2(1)/a))*a-a<<std::endl;
    std::cout<<Norm(a*((T2(1)/a)*a)-a)<<"  "<<eps*normm<<std::endl;
    std::cout<<Norm((a*(T2(1)/a))*a-a)<<"  "<<eps*normm<<std::endl;
    std::cout<<Norm(a*((T2(1)%a)*a)-a)<<"  "<<eps*normm<<std::endl;
    std::cout<<Norm((a*(T2(1)%a))*a-a)<<"  "<<eps*normm<<std::endl;
    std::cout<<Norm((T2(1)/a)*a-Adjoint((T2(1)/a)*a))<<"  "<<eps<<std::endl;
    std::cout<<Norm((T2(1)%a)*a-Adjoint((T2(1)%a)*a))<<"  "<<eps<<std::endl;
    std::cout<<Norm(a*(T2(1)/a)-Adjoint(a*(T2(1)/a)))<<"  "<<eps<<std::endl;
    std::cout<<Norm(a*(T2(1)%a)-Adjoint(a*(T2(1)%a)))<<"  "<<eps<<std::endl;
  }
  Assert(Norm(a*((T2(1)/a)*a)-a) <= eps*normm,label+" a*(1/a)*a");
  Assert(Norm((a*(T2(1)/a))*a-a) <= eps*normm,label+" a*(1/a)*a");
  Assert(Norm(a*((T2(1)%a)*a)-a) <= eps*normm,label+" a*(1%a)*a");
  Assert(Norm((a*(T2(1)%a))*a-a) <= eps*normm,label+" a*(1%a)*a");
  Assert(Norm((T2(1)/a)*a-Adjoint((T2(1)/a)*a)) <= eps,
      label+" (1/a)*a-((1/a)*a)T");
  Assert(Norm((T2(1)%a)*a-Adjoint((T2(1)%a)*a)) <= eps,
      label+" (1%a)*a-((1%a)*a)T");
  Assert(Norm(a*(T2(1)/a)-Adjoint(a*(T2(1)/a))) <= eps,
      label+" a*(1/a)-(a*(1/a))T");
  Assert(Norm(a*(T2(1)%a)-Adjoint(a*(T2(1)%a))) <= eps,
      label+" a*(1%a)-(a*(1%a))T");

  if (showstartdone) 
    std::cout<<"Done MXb"<<std::endl;
}

template <class SM1, class T, class T2> inline void DoTestMatrixDivArithMXa(
    tmv::DivType dt, const SM1& a,
    const tmv::Matrix<T>& m, T2 x, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MXa: "<<label<<std::endl;
  }
  a.DivideUsing(dt);
  m.DivideUsing(dt==tmv::CH?tmv::LU:dt);
  DoTestMatrixDivArithMXb(a,m,x,label+" "+tmv::Text(dt));
  if (showstartdone) 
    std::cout<<"Done MXa"<<std::endl;
}

template <class SM1, class T, class T2> inline void DoTestMatrixDivArithMX(
    tmv::DivType dt, const SM1& a,
    const tmv::Matrix<T>& m, T2 x, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start MX: "<<label<<std::endl;
  }
  tmv::Matrix<T> mt = m.Transpose();
  DoTestMatrixDivArithMXa(dt,a.View(),m,x,label);
  DoTestMatrixDivArithMXa(dt,Transpose(a),mt,x,label+" Trans");

  if (tmv::IsComplex(T())) {
    tmv::Matrix<T> mc = m.Conjugate();
    tmv::Matrix<T> ma = m.Adjoint();
    DoTestMatrixDivArithMXa(dt,Conjugate(a),mc,x,label+" Conj");
    DoTestMatrixDivArithMXa(dt,Adjoint(a),ma,x,label+" Adj");
  }
  if (showstartdone) 
    std::cout<<"Done MX"<<std::endl;
}

template <class T, class SM1, class SM2, class CSM1, class CSM2> 
inline void TestMatrixDivArith1(
    tmv::DivType dt, const SM1& a, const SM2& b,
    const CSM1& ca, const CSM2& cb, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start Test Div 1: "<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<"  "<<b<<std::endl;
    std::cout<<"ca = "<<tmv::Type(ca)<<"  "<<ca<<std::endl;
    std::cout<<"cb = "<<tmv::Type(cb)<<"  "<<cb<<std::endl;
  }
  
  a.SaveDiv();
  b.SaveDiv();
  ca.SaveDiv();
  cb.SaveDiv();

  tmv::Matrix<T> m1 = a;
  tmv::Matrix<T> m2 = b;
  tmv::Matrix<CT> cm1 = ca;
  tmv::Matrix<CT> cm2 = cb;

  m1.SaveDiv();
  m2.SaveDiv();
  cm1.SaveDiv();
  cm2.SaveDiv();

#ifdef XTEST
  DoTestMatrixDivArithMM1(dt,b.View(),a.View(),m2,m1,label+" R,R");
  DoTestMatrixDivArithMM2(dt,b.View(),a.View(),m2,m1,label+" R,R");
  DoTestMatrixDivArithMM1(dt,b.View(),ca.View(),m2,cm1,label+" R,C");
  DoTestMatrixDivArithMM1(dt,cb.View(),a.View(),cm2,m1,label+" C,R");
  DoTestMatrixDivArithMM2(dt,cb.View(),a.View(),cm2,m1,label+" C,R");
#endif
  DoTestMatrixDivArithMM1(dt,cb.View(),ca.View(),cm2,cm1,label+" C,C");
  DoTestMatrixDivArithMM2(dt,cb.View(),ca.View(),cm2,cm1,label+" C,C");
}

template <class T, class SM1, class SM2, class CSM1, class CSM2> 
inline void TestMatrixDivArith2(
    tmv::DivType dt, const SM1& a, const SM2& b,
    const CSM1& ca, const CSM2& cb, std::string label)
{
  if (showstartdone) {
    std::cout<<"Start Test Div 2: "<<label<<std::endl;
    std::cout<<"a = "<<tmv::Type(a)<<"  "<<a<<std::endl;
    std::cout<<"b = "<<tmv::Type(b)<<"  "<<b<<std::endl;
    std::cout<<"ca = "<<tmv::Type(ca)<<"  "<<ca<<std::endl;
    std::cout<<"cb = "<<tmv::Type(cb)<<"  "<<cb<<std::endl;
  }
  
  a.SaveDiv();
  b.SaveDiv();
  ca.SaveDiv();
  cb.SaveDiv();

  tmv::Matrix<T> m1 = a;
  tmv::Matrix<T> m2 = b;
  tmv::Matrix<CT> cm1 = ca;
  tmv::Matrix<CT> cm2 = cb;

  m1.SaveDiv();
  m2.SaveDiv();
  cm1.SaveDiv();
  cm2.SaveDiv();

  CT z1(9,-2);
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

