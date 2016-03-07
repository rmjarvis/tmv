

// Preprocessor options:
// NONSQUARE    Don't do anything that requires a square matrix to compile
//              such as += x, etc.
// NODIV        Don't do Det, Norm2, etc.
// NOSV         Only skip division stuff that requires SVD, like Norm2.
// ALIASOK      Include alised arithmetic, like m = m + x.
// SYMOPROD     Require outer product statements to be symmetric (v1 ^ v1).
// NOADDEQ      Skip m += m2 and similar.
// NOADDEQX     Skip m += x and similar.
// NOMULTEQ     Skip m *= m2 and similar.
// NOMULTEQX    Skip m *= x and similar.
// NOELEMMULT   Skip ElemMult(m1,m2)
// BASIC_MULTMM_ONLY Skip all matrix multiplication other than simple c=a*b.
// NO_REAL_ARITH     Skip all real tests.
// NO_COMPLEX_ARITH  Skip all complex tests.

template <class M1, class M2> 
inline bool CanAdd(const M1& a, const M2& b)
{ return a.colsize() == b.colsize() && a.rowsize() == b.rowsize(); }

template <class M1, class M2> 
inline bool CanAddEq(const M1& a, const M2& b)
{ return CanAdd(a,b); }

template <class M, class T2> 
inline bool CanAddX(const M& a, const T2)
{ return a.isSquare(); }

template <class M, class T2> 
inline bool CanAddEqX(const M& a, const T2 x)
{ return CanAddX(a,x); }

template <class M, class T2> 
inline bool CanMultX(const M&, const T2)
{ return true; }

template <class M, class T2> 
inline bool CanMultEqX(const M& a, const T2 x)
{ return CanMultX(a,x); }

template <class M1, class M2> 
inline bool CanElemMultMM(const M1& a, const M2& b)
{ return a.colsize() == b.colsize() && a.rowsize() == b.rowsize(); }

template <class M1, class M2, class M3> 
inline bool CanAddElemMultMM(const M1& a, const M2& b, const M3& c)
{
    return CanElemMultMM(a,b) && 
        c.colsize() == b.colsize() && c.rowsize() == b.rowsize(); 
}

template <class M1, class M2> 
inline bool CanMultMM(const M1& a, const M2& b)
{ return a.rowsize() == b.colsize(); }

template <class M, class V> 
inline bool CanMultMV(const M& m, const V& v)
{ return m.rowsize() == v.size(); }

template <class M, class V> 
inline bool CanMultVM(const V& v, const M& m)
{ return m.colsize() == v.size(); }

template <class M1, class M2, class M3> 
inline bool CanMultMM(const M1& a, const M2& b, const M3& c)
{
    return CanMultMM(a,b) && c.colsize() == a.colsize() &&
        c.rowsize() == b.rowsize();
}

template <class M1, class V2, class V3> 
inline bool CanMultMV(const M1& a, const V2& b, const V3& c)
{ return CanMultMV(a,b) && c.size() == a.colsize(); }

template <class V1, class M2, class V3> 
inline bool CanMultVM(const V1& a, const M2& b, const V3& c)
{ return CanMultVM(a,b) && c.size() == b.rowsize(); }

template <class M1, class M2, class M3> 
inline bool CanMultXMM(const M1& a, const M2& b, const M3& c)
{ return CanMultMM(a,b,c); }

template <class V0, class V1>
inline void CopyBackV(const V0& v0, V1& v1)
{ v1 = v0; }

template <class M0, class M1>
inline void CopyBackM(const M0& m0, M1& m1)
{ m1 = m0; }

#define ProductType(T1,T2) typename tmv::Traits2<T1,T2>::type

// Every one of these inner routines has a _Basic and a _Full version.
// The _Basic is used to test all the various views of the matrices
// and vectors.  
// The _Full is used to test alternative names for funcions and
// some of the composite manipulations of combining the different
// coefficients and such.
// The _Full routines are a bit taxing for the compiler, so we only 
// call them on the regular view, with everything complex.  The other
// views and real types just get the _Basic calls.

template <class MM> 
static void DoTestMa_Basic(const MM& a, std::string label)
{
    typedef typename MM::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start Ma "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
    }

    tmv::Matrix<T> m = a;
    RT eps = EPS * RT(m.colsize()+m.rowsize());

#ifdef XXD
    if (XXDEBUG1) {
        std::cout<<"a = "<<tmv::TMV_Text(a)<<" = "<<a<<std::endl;
        std::cout<<"m = "<<tmv::TMV_Text(m)<<" = "<<m<<std::endl;
        std::cout<<"a-m = "<<a-m<<std::endl;
        std::cout<<"Norm(a-m) = "<<Norm(a-m)<<std::endl;
#ifndef NONSQUARE
        if (m.isSquare()) {
            std::cout<<"Trace(a) = "<<Trace(a)<<"  "<<Trace(m)<<std::endl;
#ifndef NODIV
            std::cout<<"Det(a) = "<<Det(a)<<"  "<<Det(m)<<std::endl;
            std::cout<<"diff = "<<tmv::TMV_ABS2(Det(a)-Det(m))<<"  "<<
                eps*Norm(m)*Norm(m.inverse())*tmv::TMV_ABS2(Det(m))<<std::endl;
            std::cout<<"LogDet(a) = "<<LogDet(a)<<"  "<<LogDet(m)<<std::endl;
            std::cout<<"diff = "<<tmv::TMV_ABS2(LogDet(a)-LogDet(m))<<"  "<<
                eps*Norm(m)*Norm(m.inverse())<<std::endl;
#endif
        }
#endif
        std::cout<<"NormF(a) = "<<NormF(a)<<"  "<<NormF(m)<<std::endl;
        std::cout<<"Norm(a) = "<<Norm(a)<<"  "<<Norm(m)<<std::endl;
        std::cout<<"NormSq(a) = "<<NormSq(a)<<"  "<<NormSq(m)<<std::endl;
        std::cout<<"abs(diff) = "<<tmv::TMV_ABS2(NormSq(a)-NormSq(m))<<std::endl;
        std::cout<<"eps*normsq = "<<eps*NormSq(m)<<std::endl;
        std::cout<<"Norm1(a) = "<<Norm1(a)<<"  "<<Norm1(m)<<std::endl;
        std::cout<<"NormInf(a) = "<<NormInf(a)<<"  "<<NormInf(m)<<std::endl;
        std::cout<<"abs(diff) = "<<tmv::TMV_ABS2(NormInf(a)-NormInf(m))<<std::endl;
        std::cout<<"eps*norminf = "<<eps*NormInf(m)<<std::endl;
        std::cout<<"MaxAbsElement(a) = "<<MaxAbsElement(a)<<"  "<<
            MaxAbsElement(m)<<std::endl;
    }
#endif

    Assert(Equal(a,m,eps),label+" a != m");
#ifndef NONSQUARE
    if (m.isSquare()) {
        Assert(Equal2(Trace(a),Trace(m),eps*tmv::TMV_ABS2(Trace(m))),
               label+" Trace");
#ifndef NODIV
        T d = Det(m);
        if (tmv::TMV_ABS2(d) > 0.5) {
            RT eps1 = eps;
            if (!std::numeric_limits<RT>::is_integer) 
                eps1 *= Norm(m) * Norm(m.inverse());
            Assert(Equal2(Det(a),d,eps1*tmv::TMV_ABS2(d)),label+" Det");
            if (!std::numeric_limits<RT>::is_integer) {
                Assert(Equal2(LogDet(a),LogDet(m),eps1),label+" LogDet");
            }
        } else if (tmv::TMV_ABS2(d) != 0.0) {
            RT eps1 = eps;
            if (!std::numeric_limits<RT>::is_integer) 
                eps1 *= Norm(m) * Norm(m.inverse());
            Assert(Equal2(Det(a),d,eps1*(1+tmv::TMV_ABS2(d))),label+" Det");
        } else {
            Assert(Equal2(Det(a),d,eps),label+" Det");
        }
#endif
    }
#endif

    if (!std::numeric_limits<RT>::is_integer) {
        Assert(Equal2(NormF(a),NormF(m),eps*NormF(m)),label+" NormF");
        Assert(Equal2(Norm(a),Norm(m),eps*Norm(m)),label+" Norm");
    }
    Assert(Equal2(NormSq(a),NormSq(m),eps*NormSq(m)),label+" NormSq");
    if (!std::numeric_limits<RT>::is_integer && 
        tmv::Traits<T>::iscomplex) {
        Assert(Equal2(Norm1(a),Norm1(m),eps*Norm1(m)),label+" Norm1");
        Assert(Equal2(NormInf(a),NormInf(m),eps*NormInf(m)),label+" NormInf");
        Assert(Equal2(MaxAbsElement(a),MaxAbsElement(m),eps*MaxAbsElement(m)),
               label+" MaxAbsElement");
    }
    Assert(Equal2(MaxAbs2Element(a),MaxAbs2Element(m),eps*MaxAbs2Element(m)),
           label+" MaxAbs2Element");
#ifndef NODIV
#ifndef NOSV
    if (donorm2 && !std::numeric_limits<RT>::is_integer) {
#ifdef XXD
        if (XXDEBUG1) {
            std::cout<<"Norm2(a) = "<<a.doNorm2()<<"  "<<m.doNorm2()<<std::endl;
            std::cout<<"abs(diff) = "<<tmv::TMV_ABS2(a.doNorm2()-m.doNorm2())<<std::endl;
            std::cout<<"eps*kappa = "<<eps*m.doCondition()<<std::endl;
        }
#endif
        Assert(Equal2(a.doNorm2(),m.doNorm2(),eps*m.doCondition()*m.doNorm2()),
               label+" DoNorm2");
        Assert(Equal2(Norm2(a),Norm2(m),eps*m.condition()*m.doNorm2()),
               label+" Norm2");
    }
#endif
#endif
#ifdef XXD
    if (XXDEBUG1) {
        std::cout<<"Norm(aT-mT) = "<<Norm(Transpose(a)-Transpose(m))<<std::endl;
        std::cout<<"Conjugate(a) = "<<Conjugate(a)<<std::endl;
        std::cout<<"Conjugate(m) = "<<Conjugate(m)<<std::endl;
        std::cout<<"a*-m* = "<<Conjugate(a)-Conjugate(m)<<std::endl;
        std::cout<<"Norm(a*-m*) = "<<Norm(Conjugate(a)-Conjugate(m))<<std::endl;
        std::cout<<"Adjoint(a) = "<<Adjoint(a)<<std::endl;
        std::cout<<"Adjoint(m) = "<<Adjoint(m)<<std::endl;
        std::cout<<"Norm(at-mt) = "<<Norm(Adjoint(a)-Adjoint(m))<<std::endl;
    }
#endif
    if (!std::numeric_limits<RT>::is_integer) eps *= Norm(m);
    Assert(Equal(Transpose(a),Transpose(m),eps),label+" Transpose");
    Assert(Equal(Conjugate(a),Conjugate(m),eps),label+" Conjugate");
    Assert(Equal(Adjoint(a),Adjoint(m),eps),label+" Adjoint");

    if (showstartdone) std::cout<<"Done Ma"<<std::endl;
}

template <class MM> 
static void DoTestMa_Full(const MM& a, std::string label)
{
    DoTestMa_Basic(a,label);

#if (XTEST & 16)
    typedef typename MM::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Continuing on to MMa Full"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
    }
    tmv::Matrix<T> m = a;
    RT eps = EPS * (a.colsize() + a.rowsize());

#ifdef XXD
    if (XXDEBUG1) {
        std::cout<<"a = "<<tmv::TMV_Text(a)<<" = "<<a<<std::endl;
        std::cout<<"m = "<<tmv::TMV_Text(m)<<" = "<<m<<std::endl;
        std::cout<<"a.normF() = "<<a.normF()<<"  "<<NormF(m)<<std::endl;
        std::cout<<"a.norm() = "<<a.norm()<<"  "<<Norm(m)<<std::endl;
        std::cout<<"a.normSq() = "<<a.normSq()<<"  "<<NormSq(m)<<std::endl;
        std::cout<<"a.normSq(1.e-3) = "<<a.normSq(RT(1.e-3))<<"  "<<NormSq(m*RT(1.e-3))<<std::endl;
        std::cout<<"abs(diff) = "<<tmv::TMV_ABS2(a.normSq(RT(1.e-3))-NormSq(m*RT(1.e-3)))<<std::endl;
        std::cout<<"eps*1.e-6*normsq = "<<eps*1.e-6*NormSq(m)<<std::endl;
        std::cout<<"a.norm1() = "<<a.norm1()<<"  "<<Norm1(m)<<std::endl;
        std::cout<<"a.normInf() = "<<a.normInf()<<"  "<<NormInf(m)<<std::endl;
        std::cout<<"a.sumAbsElements() = "<<a.sumAbsElements()<<"  "<<SumAbsElements(m)<<std::endl;
        std::cout<<"a.sumAbs2Elements() = "<<a.sumAbs2Elements()<<"  "<<SumAbs2Elements(m)<<std::endl;
    }
#endif
    if (!(std::numeric_limits<RT>::is_integer)) {
        Assert(Equal2(a.normF(),NormF(m),eps*NormF(m)),label+" NormF");
        Assert(Equal2(a.norm(),Norm(m),eps*Norm(m)),label+" Norm");
    }
    Assert(Equal2(a.sumElements(),SumElements(m),
                  eps*tmv::TMV_ABS2(SumElements(m))),
           label+" SumElements");
    Assert(Equal2(a.normSq(),NormSq(m),eps*NormSq(m)),label+" NormSq");
    Assert(Equal2(a.normSq(1000),NormSq(m*RT(1000)),eps*1000000*NormSq(m)),
           label+" NormSq,scale");
    if (!std::numeric_limits<RT>::is_integer && tmv::Traits<T>::iscomplex) {
        Assert(Equal2(a.sumAbsElements(),SumAbsElements(m),
                      eps*SumAbsElements(m)),
               label+" SumAbsElements");
        Assert(Equal2(a.normSq(RT(1.e-3)),NormSq(m*RT(1.e-3)),
                      eps*RT(1.e-6)*NormSq(m)),
               label+" NormSq,scale");
        Assert(Equal2(a.norm1(),Norm1(m),eps*Norm1(m)),label+" Norm1");
        Assert(Equal2(a.normInf(),NormInf(m),eps*NormInf(m)),label+" NormInf");
        Assert(Equal2(a.maxAbsElement(),MaxAbsElement(m),eps*MaxAbsElement(m)),
               label+" MaxAbsElement");
    }
    Assert(Equal2(a.sumAbs2Elements(),SumAbs2Elements(m),eps*SumAbs2Elements(m)),
           label+" SumAbs2Elements");
    Assert(Equal2(a.maxAbs2Element(),MaxAbs2Element(m),eps*MaxAbs2Element(m)),
           label+" MaxAbs2Element");
#ifndef NODIV
#ifndef NOSV
    if (donorm2 && !std::numeric_limits<RT>::is_integer) {
        Assert(Equal2(a.norm2(),m.doNorm2(),eps*m.condition()*m.doNorm2()),
               label+" Norm2");
    }
#endif
#endif

    if (showstartdone) std::cout<<"Done Ma_Full"<<std::endl;
#endif
}

template <class MM> 
static void DoTestMR(const MM& a, std::string label)
{
    DoTestMa_Basic(a,label);
    DoTestMa_Basic(Transpose(a),label+" Trans");
}

template <class MM> 
static void DoTestMC(const MM& a, std::string label)
{
    DoTestMa_Full(a,label);
    DoTestMa_Basic(Transpose(a),label+" Trans");
    DoTestMa_Basic(Conjugate(a),label+" Conj");
    DoTestMa_Basic(Adjoint(a),label+" Adj");
}

template <class MM, class T2> 
static void DoTestMX1a_Basic(const MM& a, T2 x, std::string label)
{
    typedef typename MM::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start MX1a"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"x = "<<tmv::TMV_Text(x)<<" "<<x<<std::endl;
    }

    tmv::Matrix<T> m = a;

    RT eps = EPS * (a.colsize() + a.rowsize());
    if (!std::numeric_limits<RT>::is_integer) eps *= Norm(m);

#ifndef NONSQUARE
    if (CanAddX(a,x)) {
#ifdef XXD
        if (XXDEBUG2) {
            std::cout<<"CanAddX("<<tmv::TMV_Text(a)<<","<<tmv::TMV_Text(x)<<")\n";
            std::cout<<"x = "<<tmv::TMV_Text(x)<<"  "<<x<<std::endl;
            std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
            std::cout<<"x-a = "<<(x-a)<<std::endl;
            std::cout<<"x-m = "<<x-m<<std::endl;
            std::cout<<"(x-a)-(x-m) = "<<(x-a)-(x-m)<<std::endl;
            std::cout<<"Norm((x-a)-(x-m)) = "<<Norm((x-a)-(x-m))<<std::endl;
        }
#endif
        Assert(Equal(a,m,eps),label+" a != m1");
        Assert(Equal(x-a,x-m,eps),label+" x-a");
    }
#endif
    if (CanMultX(a,x)) {
#ifdef XXD
        if (XXDEBUG2) {
            std::cout<<"CanMultX("<<tmv::TMV_Text(a)<<","<<tmv::TMV_Text(x)<<")\n";
        }
#endif
        Assert(Equal(x*a,x*m,eps*tmv::TMV_ABS2(x)),label+" x*a");
    }
    if (showstartdone) std::cout<<"Done MX1a"<<std::endl;
}

template <class MM, class T2> 
static void DoTestMX1a_Full(const MM& a, T2 x, std::string label)
{
    DoTestMX1a_Basic(a,x,label);

#if (XTEST & 16)
    typedef typename MM::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Continuing on to MX1a Full"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"x = "<<tmv::TMV_Text(x)<<" "<<x<<std::endl;
    }
    tmv::Matrix<T> m = a;

    RT eps = EPS * (a.colsize() + a.rowsize());
    if (!std::numeric_limits<RT>::is_integer) eps *= Norm(m);

#ifndef NONSQUARE
    if (CanAddX(a,x)) {
        Assert(Equal(a-x,m-x,eps),label+" a-x");
        Assert(Equal(x+a,x+m,eps),label+" x+a");
        Assert(Equal(a+x,m+x,eps),label+" a+x");
    }
#endif
    if (CanMultX(a,x)) {
        Assert(Equal(a*x,x*m,eps*tmv::TMV_ABS2(x)),label+" a*x");
        if (!std::numeric_limits<RT>::is_integer) {
            Assert(Equal(a/x,m/x,eps/tmv::TMV_ABS2(x)),label+" a/x");
        }
    }
    if (showstartdone) std::cout<<"Done MX1a_Full"<<std::endl;
#endif
}

template <class MM, class T2> 
static void DoTestMX1R(const MM& a, T2 x, std::string label)
{
    DoTestMX1a_Basic(a,x,label);
}

template <class MM, class T2> 
static void DoTestMX1C(const MM& a, T2 x, std::string label)
{
    DoTestMX1a_Full(a,x,label);
    DoTestMX1a_Basic(Conjugate(a),x,label+" Conj");
}

template <class MM, class T2> 
static void DoTestMX2a_Basic(MM& a, T2 x, std::string label)
{
    typedef typename MM::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start MX2a"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"x = "<<tmv::TMV_Text(x)<<" "<<x<<std::endl;
    }
#ifdef XXD
    if (XXDEBUG3) {
        std::cout<<"a = "<<a<<std::endl;
    }
#endif

    tmv::Matrix<T> m1 = a;
    tmv::Matrix<T> temp = a;

    RT eps = EPS * (a.colsize() + a.rowsize());
    if (!std::numeric_limits<RT>::is_integer) eps *= Norm(m1);

    Assert(Equal(a,m1,eps),label+" a = m1");

#ifndef NOADDEQX
#ifndef NONSQUARE
    if (CanAddEqX(a,x)) {
        typename MM::copy_type a0(a);
#ifdef XXD
        if (XXDEBUG3) {
            std::cout<<"CanAddEqX("<<tmv::TMV_Text(a)<<","<<tmv::TMV_Text(x)<<")\n";
        }
#endif
        tmv::Matrix<T> m2 = m1;
        a += x;
        m2 = m1+x;
        Assert(Equal(a,m2,eps),label+" a += x");
        Assert(Equal(a+=x,m2+=x,eps),label+" a += x (2)");
        CopyBackM(a0,a);
        a = a+x; 
        m2 = m1+x;
        Assert(Equal(a,m2,eps),label+" a = a+x");
        CopyBackM(a0,a);
    }
#endif
#endif
#ifndef NOMULTEQX
    if (CanMultEqX(a,x)) {
        typename MM::copy_type a0(a);
#ifdef XXD
        if (XXDEBUG3) {
            std::cout<<"CanMultEqX("<<tmv::TMV_Text(a)<<","<<tmv::TMV_Text(x)<<")\n";
        }
#endif
        tmv::Matrix<T> m2 = m1;
        a *= x;
        m2 = m1*x;
#ifdef XXD
        if (XXDEBUG3) {
            std::cout<<"a *= x = "<<a<<std::endl;
            std::cout<<"m2 = "<<m2<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(a-m2)<<"  <=? "<<eps*tmv::TMV_ABS2(x)<<std::endl;
        }
#endif
        Assert(Equal(a,m2,eps*tmv::TMV_ABS2(x)),label+" a *= x");
        Assert(Equal(a*=x,m2*=x,eps*tmv::TMV_ABS2(x*x)),label+" a *= x");
        CopyBackM(a0,a);
        a = a*x;
        m2 = m1*x;
        Assert(Equal(a,m2,eps*tmv::TMV_ABS2(x)),label+" a = a*x");
        CopyBackM(a0,a);
    }
#endif

    if (showstartdone) std::cout<<"Done MX2a"<<std::endl;
}

template <class MM, class T2> 
static void DoTestMX2a_Full(MM& a, T2 x, std::string label)
{
    DoTestMX2a_Basic(a,x,label);

#if (XTEST & 16)
    typedef typename MM::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Continuing on to MX2a Full"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"x = "<<tmv::TMV_Text(x)<<" "<<x<<std::endl;
    }
    tmv::Matrix<T> m1 = a;
    tmv::Matrix<T> temp = a;

    RT eps = EPS * (a.colsize() + a.rowsize());
    if (!std::numeric_limits<RT>::is_integer) eps *= Norm(m1);

    Assert(Equal(a,m1,eps),label+" a = m1");

#ifndef NOADDEQX
#ifndef NONSQUARE
    if (CanAddEqX(a,x)) {
        typename MM::copy_type a0(a);
        tmv::Matrix<T> m2 = m1;
        a += -x;
        m2 = m1-x;
        Assert(Equal(a,m2,eps),label+" a += x");
        Assert(Equal(a+=-x,m2+=-x,eps),label+" a += -x (2)");
        CopyBackM(a0,a);
        a -= x;
        m2 = m1-x;
        Assert(Equal(a,m2,eps),label+" a -= x");
        Assert(Equal(a-=x,m2-=x,eps),label+" a -= x (2)");
        CopyBackM(a0,a);
        a -= -x;
        m2 = m1+x;
        Assert(Equal(a,m2,eps),label+" a -= x");
        Assert(Equal(a-=-x,m2-=-x,eps),label+" a -= -x (2)");
        CopyBackM(a0,a);
        a = a-x;
        m2 = m1-x;
        Assert(Equal(a,m2,eps),label+" a = a-x");
        CopyBackM(a0,a);
        a = x+a;
        m2 = x+m1;
        Assert(Equal(a,m2,eps),label+" a = x+a");
        CopyBackM(a0,a);
        a = x-a; 
        m2 = x-m1;
        Assert(Equal(a,m2,eps),label+" a = x-a");
        CopyBackM(a0,a);
    }
#endif
#endif
#ifndef NOMULTEQX
    if (CanMultEqX(a,x)) {
        typename MM::copy_type a0(a);
        tmv::Matrix<T> m2 = m1;
        a *= -x;
        m2 = -m1*x;
        Assert(Equal(a,m2,eps*tmv::TMV_ABS2(x)),label+" a *= -x");
        Assert(Equal(a*=-x,m2*=-x,eps*tmv::TMV_ABS2(x*x)),label+" a *= -x");
        CopyBackM(a0,a);
        if (!std::numeric_limits<RT>::is_integer) {
            a /= x;
            m2 = m1/x;
            Assert(Equal(a,m2,eps*tmv::TMV_ABS2(x)),label+" a /= x");
            Assert(Equal(a/=x,m2/=x,eps),label+" a /= x");
            CopyBackM(a0,a);
            a /= -x;
            m2 = -m1/x;
            Assert(Equal(a,m2,eps*tmv::TMV_ABS2(x)),label+" a /= -x");
            Assert(Equal(a/=-x,m2/=-x,eps),label+" a /= -x");
            CopyBackM(a0,a);
            a = a/x;
            m2 = m1/x;
            Assert(Equal(a,m2,eps),label+" a = a/x");
            CopyBackM(a0,a);
        }
        a = x*a; 
        m2 = x*m1;
        Assert(Equal(a,m2,eps*tmv::TMV_ABS2(x)),label+" a = x*a");
        CopyBackM(a0,a);
    }
#endif

    if (showstartdone) std::cout<<"Done MX2a_Full"<<std::endl;
#endif
}

template <class MM, class T2> 
static void DoTestMX2R(MM& a, T2 x, std::string label)
{
    DoTestMX2a_Basic(a,x,label);
}

template <class MM, class T2> 
static void DoTestMX2C(MM& a, T2 x, std::string label)
{
    typename MM::conjugate_type ac = a.conjugate();
    DoTestMX2a_Full(a,x,label);
    DoTestMX2a_Basic(ac,x,label+" Conj");
}

template <class MM, class V> 
static void DoTestMV1a_Basic(const MM& a, const V& b, std::string label)
{
    typedef typename MM::value_type Ta;
    typedef typename V::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    typedef typename tmv::Traits<T>::complex_type CT;
    if (showstartdone) {
        std::cout<<"Start MV1a"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }
#ifdef XXD
    if (XXDEBUG4) {
        std::cout<<"a = "<<a<<std::endl;
        std::cout<<"b = "<<b<<std::endl;
    }
#endif

    tmv::Matrix<Ta> m = a;
    tmv::Vector<T> v = b;

    RT eps = EPS * 
        (a.colsize() + a.rowsize()) * (a.colsize() + a.rowsize());
    if (!std::numeric_limits<RT>::is_integer) 
        eps *= Norm(m) * Norm(v);

    if (CanMultMV(a,b)) {
#ifdef XXD
        if (XXDEBUG4) {
            std::cout<<"CanMult("<<tmv::TMV_Text(m)<<","<<tmv::TMV_Text(b)<<")\n";
            std::cout<<"a*b = "<<(a*b)<<std::endl;
            std::cout<<"m*v = "<<m*v<<std::endl;
            std::cout<<"a*b-m*v = "<<(a*b-(m*v))<<std::endl;
            std::cout<<"Norm(a*b-m*v) = "<<Norm(a*b-(m*v))<<std::endl;
            std::cout<<"eps = "<<eps<<std::endl;
        }
#endif
        Assert(Equal(a*b,m*v,eps),label+" a*b");
        RT x(5);
        CT z(3,4);
#ifdef XXD
        if (XXDEBUG4) {
            std::cout<<"x*a*b = "<<(x*a*b)<<std::endl;
            std::cout<<"x*m*v = "<<x*m*v<<std::endl;
            std::cout<<"diff = "<<(x*a*b-x*m*v)<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(x*a*b-x*m*v)<<std::endl;
            std::cout<<"x*eps = "<<x*eps<<std::endl;
        }
#endif
        Assert(Equal(x*a*b,x*m*v,x*eps),label+" x*a*b");
#ifdef XXD
        if (XXDEBUG4) {
            std::cout<<"z*a*b = "<<(z*a*b)<<std::endl;
            std::cout<<"z*m*v = "<<z*m*v<<std::endl;
            std::cout<<"diff = "<<(z*a*b-z*m*v)<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(z*a*b-z*m*v)<<std::endl;
            std::cout<<"x*eps = "<<x*eps<<std::endl;
        }
#endif
        Assert(Equal(z*a*b,z*m*v,x*eps),label+" z*a*b");
    }
    if (showstartdone) std::cout<<"Done MV1a"<<std::endl;
}

template <class MM, class V> 
static void DoTestMV1a_Full(const MM& a, const V& b, std::string label)
{
    DoTestMV1a_Basic(a,b,label);

#if (XTEST & 16)
    typedef typename MM::value_type Ta;
    typedef typename V::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    typedef typename tmv::Traits<T>::complex_type CT;
    if (showstartdone) {
        std::cout<<"Continuing on to MV1a Full"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }
    tmv::Matrix<Ta> m = a;
    tmv::Vector<T> v = b;

    RT eps = EPS * (a.colsize() + a.rowsize());
    if (!std::numeric_limits<RT>::is_integer) 
        eps *= Norm(m) * Norm(v);

    if (CanMultMV(a,b)) {
        RT x(5);
        std::complex<RT> z(3,4);
        Assert(Equal((x*a)*b,(x*m*v),x*eps),label+" (x*a)*b");
        Assert(Equal(x*(a*b),(x*m*v),x*eps),label+" x*(a*b)");
        Assert(Equal(a*(x*b),(x*m*v),x*eps),label+" a*(x*b)");
        Assert(Equal((z*a)*b,(z*m*v),x*eps),label+" (z*a)*b");
        Assert(Equal(z*(a*b),(z*m*v),x*eps),label+" z*(a*b)");
        Assert(Equal(a*(z*b),(z*m*v),x*eps),label+" a*(z*b)");
    }
    if (showstartdone) std::cout<<"Done MV1a_Full"<<std::endl;
#endif
}

template <class MM, class V> 
static void DoTestVM1a_Basic(const MM& a, const V& b, std::string label)
{
    typedef typename MM::value_type Ta;
    typedef typename V::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    typedef typename tmv::Traits<T>::complex_type CT;
    if (showstartdone) {
        std::cout<<"Start VM1a"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }
#ifdef XXD
    if (XXDEBUG4) {
        std::cout<<"a = "<<a<<std::endl;
        std::cout<<"b = "<<b<<std::endl;
    }
#endif

    tmv::Matrix<Ta> m = a;
    tmv::Vector<T> v = b;

    RT eps = EPS * (a.colsize() + a.rowsize());
    if (!std::numeric_limits<RT>::is_integer) 
        eps *= Norm(m) * Norm(v);

    if (CanMultVM(b,a)) {
        RT x(5);
        CT z(3,4);
#ifdef XXD
        if (XXDEBUG4) {
            std::cout<<"CanMult("<<tmv::TMV_Text(v)<<","<<tmv::TMV_Text(m)<<")\n";
            std::cout<<"v*m = "<<v*m<<std::endl;
            std::cout<<"b*a = "<<(b*a)<<std::endl;
            std::cout<<"x*b*a = "<<(x*b*a)<<std::endl;
            std::cout<<"z*b*a = "<<(z*b*a)<<std::endl;
            std::cout<<"b*a Norm(diff) = "<<Norm(b*a-v*m)<<std::endl;
            std::cout<<"x*b*a Norm(diff) = "<<Norm(x*b*a-x*v*m)<<std::endl;
            std::cout<<"z*b*a Norm(diff) = "<<Norm(z*b*a-z*v*m)<<std::endl;
            std::cout<<"eps = "<<eps<<std::endl;
            std::cout<<"x*eps = "<<x*eps<<std::endl;
        }
#endif
        Assert(Equal(b*a,v*m,eps),label+" b*a");
        Assert(Equal(x*b*a,x*v*m,x*eps),label+" x*b*a");
        Assert(Equal(z*b*a,z*v*m,x*eps),label+" z*b*a");
    }
    if (showstartdone) std::cout<<"Done VM1a"<<std::endl;
}

template <class MM, class V> 
static void DoTestVM1a_Full(const MM& a, const V& b, std::string label)
{
    DoTestVM1a_Basic(a,b,label);

#if (XTEST & 16)
    typedef typename MM::value_type Ta;
    typedef typename V::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    typedef typename tmv::Traits<T>::complex_type CT;
    if (showstartdone) {
        std::cout<<"Continuing on to VM1a Full"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }
    tmv::Matrix<Ta> m = a;
    tmv::Vector<T> v = b;

    RT eps = EPS * (a.colsize() + a.rowsize());
    if (!std::numeric_limits<RT>::is_integer) 
        eps *= Norm(m) * Norm(v);

    if (CanMultVM(b,a)) {
        RT x(5);
        CT z(3,4);
        Assert(Equal((x*b)*a,(x*v*m),x*eps),label+" (x*b)*a");
        Assert(Equal(x*(b*a),(x*v*m),x*eps),label+" x*(b*a)");
        Assert(Equal(b*(x*a),(x*v*m),x*eps),label+" b*(x*a)");
        Assert(Equal((z*b)*a,(z*v*m),x*eps),label+" (z*b)*a");
        Assert(Equal(z*(b*a),(z*v*m),x*eps),label+" z*(b*a)");
        Assert(Equal(b*(z*a),(z*v*m),x*eps),label+" b*(z*a)");
    }
    if (showstartdone) std::cout<<"Done VM1a_Full"<<std::endl;
#endif
}

template <class MM, class V> 
static void DoTestMV1R(const MM& a, V& b, std::string label)
{
    DoTestMV1a_Basic(a,b,label);
    DoTestMV1a_Basic(a,b.reverse(),label);

#if (XTEST & 8)
    typename V::copy_type b0 = b;

    b.setZero();
    DoTestMV1a_Basic(a,b,label);

    CopyBackV(b0,b);
    b.cSubVector(0,b.size()/2).setZero();
    DoTestMV1a_Basic(a,b,label);

    CopyBackV(b0,b);
    b.cSubVector(b.size()/2,b.size()).setZero();
    DoTestMV1a_Basic(a,b,label);

    CopyBackV(b0,b);
    b.cSubVector(0,b.size()/4).setZero();
    b.cSubVector(3*b.size()/4,b.size()).setZero();
    DoTestMV1a_Basic(a,b,label);

    if (b.size() > 1) {
        CopyBackV(b0,b);
        b.cSubVector(0,1).setZero();
        DoTestMV1a_Basic(a,b,label);

        CopyBackV(b0,b);
        b.cSubVector(b.size()-1,b.size()).setZero();
        DoTestMV1a_Basic(a,b,label);

        CopyBackV(b0,b);
        b.cSubVector(0,1).setZero();
        b.cSubVector(b.size()-1,b.size()).setZero();
        DoTestMV1a_Basic(a,b,label);
    }
    CopyBackV(b0,b);
#endif
}

template <class MM, class V> 
static void DoTestVM1R(const MM& a, V& b, std::string label)
{
    DoTestVM1a_Basic(a,b,label);
    DoTestVM1a_Basic(a,b.reverse(),label);

#if (XTEST & 8)
    typename V::copy_type b0 = b;

    b.setZero();
    DoTestVM1a_Basic(a,b,label);

    CopyBackV(b0,b);
    b.cSubVector(0,b.size()/2).setZero();
    DoTestVM1a_Basic(a,b,label);

    CopyBackV(b0,b);
    b.cSubVector(b.size()/2,b.size()).setZero();
    DoTestVM1a_Basic(a,b,label);

    CopyBackV(b0,b);
    b.cSubVector(0,b.size()/4).setZero();
    b.cSubVector(3*b.size()/4,b.size()).setZero();
    DoTestVM1a_Basic(a,b,label);

    if (b.size() > 1) {
        CopyBackV(b0,b);
        b.cSubVector(0,1).setZero();
        DoTestVM1a_Basic(a,b,label);

        CopyBackV(b0,b);
        b.cSubVector(b.size()-1,b.size()).setZero();
        DoTestVM1a_Basic(a,b,label);

        CopyBackV(b0,b);
        b.cSubVector(0,1).setZero();
        b.cSubVector(b.size()-1,b.size()).setZero();
        DoTestVM1a_Basic(a,b,label);
    }
    CopyBackV(b0,b);
#endif
}

template <class MM, class V> 
static void DoTestMV1C(const MM& a, V& b, std::string label)
{
    DoTestMV1a_Full(a,b,label);
    DoTestMV1a_Basic(a,b.reverse(),label);

    DoTestMV1a_Basic(Conjugate(a),b,label+" Conj");
#if (XTEST & 1)
    DoTestMV1a_Basic(Conjugate(a),b.reverse(), label+" Conj");
#endif

#if (XTEST & 8)
    typename V::copy_type b0 = b;

    b.setZero();
    DoTestMV1a_Basic(a,b,label);
    DoTestMV1a_Basic(Conjugate(a),b,label+" Conj");

    CopyBackV(b0,b);
    b.cSubVector(0,b.size()/2).setZero();
    DoTestMV1a_Basic(a,b,label);
    DoTestMV1a_Basic(Conjugate(a),b,label+" Conj");

    CopyBackV(b0,b);
    b.cSubVector(b.size()/2,b.size()).setZero();
    DoTestMV1a_Basic(a,b,label);
    DoTestMV1a_Basic(Conjugate(a),b,label+" Conj");

    CopyBackV(b0,b);
    b.cSubVector(0,b.size()/4).setZero();
    b.cSubVector(3*b.size()/4,b.size()).setZero();
    DoTestMV1a_Basic(a,b,label);
    DoTestMV1a_Basic(Conjugate(a),b,label+" Conj");

    if (b.size() > 1) {
        CopyBackV(b0,b);
        b.cSubVector(0,1).setZero();
        DoTestMV1a_Basic(a,b,label);
        DoTestMV1a_Basic(Conjugate(a),b, label+" Conj");

        CopyBackV(b0,b);
        b.cSubVector(b.size()-1,b.size()).setZero();
        DoTestMV1a_Basic(a,b,label);
        DoTestMV1a_Basic(Conjugate(a),b, label+" Conj");

        CopyBackV(b0,b);
        b.cSubVector(0,1).setZero();
        b.cSubVector(b.size()-1,b.size()).setZero();
        DoTestMV1a_Basic(a,b,label);
        DoTestMV1a_Basic(Conjugate(a),b, label+" Conj");
    }
    CopyBackV(b0,b);
#endif
}

template <class MM, class V> 
static void DoTestVM1C(const MM& a, V& b, std::string label)
{
    DoTestVM1a_Full(a,b,label);
    DoTestVM1a_Basic(a,b.reverse(),label);

    DoTestVM1a_Basic(Conjugate(a),b, label+" Conj");
#if (XTEST & 1)
    DoTestVM1a_Basic(Conjugate(a),b.reverse(), label+" Conj");
#endif

#if (XTEST & 8)
    typename V::copy_type b0 = b;

    b.setZero();
    DoTestVM1a_Basic(a,b,label);
    DoTestVM1a_Basic(Conjugate(a),b,label+" Conj");

    CopyBackV(b0,b);
    b.cSubVector(0,b.size()/2).setZero();
    DoTestVM1a_Basic(a,b,label);
    DoTestVM1a_Basic(Conjugate(a),b,label+" Conj");

    CopyBackV(b0,b);
    b.cSubVector(b.size()/2,b.size()).setZero();
    DoTestVM1a_Basic(a,b,label);
    DoTestVM1a_Basic(Conjugate(a),b,label+" Conj");

    CopyBackV(b0,b);
    b.cSubVector(0,b.size()/4).setZero();
    b.cSubVector(3*b.size()/4,b.size()).setZero();
    DoTestVM1a_Basic(a,b,label);
    DoTestVM1a_Basic(Conjugate(a),b,label+" Conj");

    if (b.size() > 1) {
        CopyBackV(b0,b);
        b.cSubVector(0,1).setZero();
        DoTestVM1a_Basic(a,b,label);
        DoTestVM1a_Basic(Conjugate(a),b, label+" Conj");

        CopyBackV(b0,b);
        b.cSubVector(b.size()-1,b.size()).setZero();
        DoTestVM1a_Basic(a,b,label);
        DoTestVM1a_Basic(Conjugate(a),b, label+" Conj");

        CopyBackV(b0,b);
        b.cSubVector(0,1).setZero();
        b.cSubVector(b.size()-1,b.size()).setZero();
        DoTestVM1a_Basic(a,b,label);
        DoTestVM1a_Basic(Conjugate(a),b, label+" Conj");
    }
    CopyBackV(b0,b);
#endif
}

template <class T> 
static inline void SetZ(T& z) { z = T(5); }
template <class T>
static inline void SetZ(std::complex<T>& z) { z = std::complex<T>(3,4); }

template <class MM, class V> 
static void DoTestMV2a_Basic(const MM& a, V& b, std::string label)
{
    typedef typename MM::value_type Ta;
    typedef typename V::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start MV2a"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }
#ifdef XXD
    if (XXDEBUG5) {
        std::cout<<"a = "<<a<<std::endl;
        std::cout<<"b = "<<b<<std::endl;
    }
#endif

    tmv::Matrix<Ta> m = a;
    tmv::Vector<T> v = b;

    RT eps = EPS * (a.colsize() + a.rowsize());
    RT eps2 = eps;
    if (!std::numeric_limits<RT>::is_integer) {
        eps *= Norm(m) * Norm(v);
        eps2 *= Norm(m*v) + Norm(m)*Norm(v);
    }

    if (CanMultMV(a,b)) {
        tmv::Vector<T> prod = m*v;
        tmv::Vector<T> c = prod;
        tmv::Vector<T> c0 = c;
        c0 /= T(2);
#ifdef XXD
        if (XXDEBUG5) {
            std::cout<<"CanMult("<<tmv::TMV_Text(m)<<","<<tmv::TMV_Text(b)<<","<<tmv::TMV_Text(v)<<")\n";
            std::cout<<"v = "<<v<<std::endl;
            std::cout<<"m = "<<m<<std::endl;
            std::cout<<"prod = "<<prod<<std::endl;
        }
#endif
        c = a*b;
#ifdef XXD
        if (XXDEBUG5) {
            std::cout<<"Norm(c-prod) = "<<Norm(c-prod)<<std::endl;
        }
#endif
        Assert(Equal(c,prod,eps),label+" c=a*b");
        c = c0;
        c += a*b;
#ifdef XXD
        if (XXDEBUG5) {
            std::cout<<"c-(c0+prod) = "<<(c-(c0+prod))<<std::endl;
            std::cout<<"Norm(c-(c0+prod)) = "<<Norm(c-(c0+prod))<<std::endl;
        }
#endif
        Assert(Equal(c,(c0+prod),eps2),label+" c+=a*b");
        c = c0;
        RT x(5);
        T z; SetZ(z);
        c = x*a*b;
#ifdef XXD
        if (XXDEBUG5) {
            std::cout<<"c = x*a*b = "<<c<<std::endl;
            std::cout<<"c-x*prod = "<<(c-x*prod)<<std::endl;
            std::cout<<"Norm(c-x*prod) = "<<Norm(c-x*prod)<<std::endl;
            std::cout<<"x*eps = "<<x*eps<<std::endl;
        }
#endif
        Assert(Equal(c,x*prod,x*eps),label+" c=x*a*b");
        c = z*a*b;
#ifdef XXD
        if (XXDEBUG5) {
            std::cout<<"c = z*a*b = "<<c<<std::endl;
            std::cout<<"c-z*prod = "<<(c-z*prod)<<std::endl;
            std::cout<<"Norm(c-z*prod) = "<<Norm(c-z*prod)<<std::endl;
        }
#endif
        Assert(Equal(c,z*prod,x*eps),label+" c=z*a*b");
        c = c0;
    }

#ifdef ALIASOK
    if (CanMultMV(a,b,b)) {
        tmv::Vector<T> prod = m*v;
        b = a*b;
#ifdef XXD
        if (XXDEBUG5) {
            std::cout<<"CanMult("<<tmv::TMV_Text(a)<<","<<tmv::TMV_Text(b)<<","<<tmv::TMV_Text(b)<<")\n";
            std::cout<<"b = a*b = "<<b<<std::endl;
            std::cout<<"b-prod = "<<(b-prod)<<std::endl;
            std::cout<<"Norm(b-prod) = "<<Norm(b-prod)<<std::endl;
            std::cout<<"eps = "<<eps<<std::endl;
        }
#endif
        Assert(Equal(b,prod,eps),label+" b=a*b");
        b = v;
        b += a*b;
        Assert(Equal(b,(v+prod),eps2),label+" b+=a*b");
        b = v;
        b -= a*b;
        Assert(Equal(b,(v-prod),eps2),label+" b-=a*b");
        b = v;
        RT x(5);
        T z; SetZ(z);
        b = x*a*b;
        Assert(Equal(b,x*prod,x*eps),label+" b=x*a*b");
        b = v;
        b = z*a*b;
        Assert(Equal(b,z*prod,x*eps),label+" b=z*a*b");
        b = v;
    }
#endif // ALIAS
    if (showstartdone) std::cout<<"Done MV2a"<<std::endl;
}

template <class MM, class V> 
static void DoTestMV2a_Full(const MM& a, V& b, std::string label)
{
    DoTestMV2a_Basic(a,b,label);

#if (XTEST & 16)
    typedef typename MM::value_type Ta;
    typedef typename V::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Continuing on to MV2a Full"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }
    tmv::Matrix<Ta> m = a;
    tmv::Vector<T> v = b;

    RT eps = EPS * (a.colsize() + a.rowsize());
    RT eps2 = eps;
    if (!std::numeric_limits<RT>::is_integer) {
        eps *= Norm(m) * Norm(v);
        eps2 *= Norm(m*v) + Norm(m)*Norm(v);
    }

    if (CanMultMV(a,b)) {
        tmv::Vector<T> prod = m*v;
        tmv::Vector<T> c = prod;
        tmv::Vector<T> c0 = c;
        c0 /= T(2);
        RT x(5);
        T z; SetZ(z);
        c = c0;
        c += x*a*b;
        Assert(Equal(c,(c0+x*prod),x*eps2),label+" c+=x*a*b");
        c = c0;
        c += z*a*b;
        Assert(Equal(c,(c0+z*prod),x*eps2),label+" c+=z*a*b");
        c = c0;
        c = -a*b;
        Assert(Equal(c,(-prod),eps),label+" c=-a*b");
        c = -x*a*b;
        Assert(Equal(c,(-x*prod),x*eps),label+" c=-x*a*b");
        c = -z*a*b;
        Assert(Equal(c,(-z*prod),x*eps),label+" c=-z*a*b");
        c = c0;
        c += -a*b;
        Assert(Equal(c,(c0-prod),eps2),label+" c+=-a*b");
        c = c0;
        c -= a*b;
        Assert(Equal(c,(c0-prod),eps2),label+" c-=a*b");
        c = c0;
        c -= -a*b;
        Assert(Equal(c,(c0+prod),eps2),label+" c-=-a*b");
        c = c0;
        c += -x*a*b;
        Assert(Equal(c,(c0-x*prod),x*eps2),label+" c+=-x*a*b");
        c = c0;
        c -= x*a*b;
        Assert(Equal(c,(c0-x*prod),x*eps2),label+" c-=x*a*b");
        c = c0;
        c -= -x*a*b;
        Assert(Equal(c,(c0+x*prod),x*eps2),label+" c-=-x*a*b");
        c = c0;
        c += -z*a*b;
        Assert(Equal(c,(c0-z*prod),x*eps2),label+" c+=-z*a*b");
        c = c0;
        c -= z*a*b;
        Assert(Equal(c,(c0-z*prod),x*eps2),label+" c-=z*a*b");
        c = c0;
        c -= -z*a*b;
        Assert(Equal(c,(c0+z*prod),x*eps2),label+" c-=-z*a*b");
        c = c0;
    }

#ifdef ALIASOK
    if (CanMultMV(a,b,b)) {
        tmv::Vector<T> prod = m*v;
        RT x(5);
        T z; SetZ(z);
        b = v;
        b += x*a*b;
        Assert(Equal(b,(v+x*prod),x*eps2),label+" b+=x*a*b");
        b = v;
        b += z*a*b;
        Assert(Equal(b,(v+z*prod),x*eps2),label+" b+=z*a*b");
        b = v;
        b = -a*b;
        Assert(Equal(b,(-prod),eps),label+" b=-a*b");
        b = v;
        b = -x*a*b;
        Assert(Equal(b,(-x*prod),x*eps),label+" b=-x*a*b");
        b = v;
        b = -z*a*b;
        Assert(Equal(b,(-z*prod),x*eps),label+" b=-z*a*b");
        b = v;
        b += -a*b;
        Assert(Equal(b,(v-prod),eps2),label+" b+=-a*b");
        b = v;
        b -= -a*b;
        Assert(Equal(b,(v+prod),eps2),label+" b-=-a*b");
        b = v;
        b += -x*a*b;
        Assert(Equal(b,(v-x*prod),x*eps2),label+" b+=-x*a*b");
        b = v;
        b -= x*a*b;
        Assert(Equal(b,(v-x*prod),x*eps2),label+" b-=x*a*b");
        b = v;
        b -= -x*a*b;
        Assert(Equal(b,(v+x*prod),x*eps2),label+" b-=-x*a*b");
        b = v;
        b += -z*a*b;
        Assert(Equal(b,(v-z*prod),x*eps2),label+" b+=-z*a*b");
        b = v;
        b -= z*a*b;
        Assert(Equal(b,(v-z*prod),x*eps2),label+" b-=z*a*b");
        b = v;
        b -= -z*a*b;
        Assert(Equal(b,(v+z*prod),x*eps2),label+" b-=-z*a*b");
        b = v;
    }
#endif // ALIAS
    if (showstartdone) std::cout<<"Done MV2a_Full"<<std::endl;
#endif
}

template <class MM, class V> 
static void DoTestVM2a_Basic(const MM& a, V& b, std::string label)
{
    typedef typename MM::value_type Ta;
    typedef typename V::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start VM2a"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }
#ifdef XXD
    if (XXDEBUG5) {
        std::cout<<"a = "<<a<<std::endl;
        std::cout<<"b = "<<b<<std::endl;
    }
#endif

    tmv::Matrix<Ta> m = a;
    tmv::Vector<T> v = b;

    RT eps = EPS * (a.colsize() + a.rowsize());
    RT eps2 = eps;
    if (!std::numeric_limits<RT>::is_integer) {
        eps *= Norm(m) * Norm(v);
        eps2 *= Norm(v*m) + Norm(m)*Norm(v);
    }

    if (CanMultVM(b,a)) {
        tmv::Vector<T> prod = v*m;
        tmv::Vector<T> c = prod;
        tmv::Vector<T> c0 = c;
        c0 /= T(2);

        c = b*a;
#ifdef XXD
        if (XXDEBUG5) {
            std::cout<<"CanMult("<<tmv::TMV_Text(b)<<","<<tmv::TMV_Text(a)<<")\n";
            std::cout<<"v = "<<v<<std::endl;
            std::cout<<"m = "<<m<<std::endl;
            std::cout<<"prod = "<<prod<<std::endl;
            std::cout<<"c = "<<c<<std::endl;
        }
#endif
        Assert(Equal(c,prod,eps),label+" c=b*a");
        c = c0;
        c += b*a;
#ifdef XXD
        if (XXDEBUG5) {
            std::cout<<"c+=b*a = "<<c<<std::endl;
            std::cout<<"c0+prod = "<<c0+prod<<std::endl;

        }
#endif
        Assert(Equal(c,(c0+prod),eps2),label+" c+=b*a");
        c = c0;
        RT x(5);
        T z; SetZ(z);
        c = x*b*a;
        Assert(Equal(c,x*prod,x*eps),label+" c=x*b*a");
        c = z*b*a;
        Assert(Equal(c,z*prod,x*eps),label+" c=z*b*a");
        c = c0;
    }

    if (CanMultVM(b,a,b)) {
        tmv::Vector<T> prod = v*m;
#ifdef XXD
        if (XXDEBUG5) {
            std::cout<<"CanMult("<<tmv::TMV_Text(b)<<","<<tmv::TMV_Text(a)<<","<<tmv::TMV_Text(b)<<")\n";
            std::cout<<"v = "<<v<<std::endl;
            std::cout<<"m = "<<m<<std::endl;
            std::cout<<"prod = "<<prod<<std::endl;
        }
#endif
        b *= a;
        Assert(Equal(b,prod,eps),label+" b*=a");
        b = v;
#ifdef ALIASOK
        b = b*a;
#ifdef XXD
        if (XXDEBUG5) {
            std::cout<<"b = b*a = "<<b<<std::endl;
            std::cout<<"b-prod = "<<(b-prod)<<std::endl;
            std::cout<<"Norm(b-prod) = "<<Norm(b-prod)<<std::endl;
            std::cout<<"eps = "<<eps<<std::endl;
        }
#endif
        Assert(Equal(b,prod,eps),label+" b=b*a");
        b = v;
        b += b*a;
        Assert(Equal(b,(v+prod),eps2),label+" b+=b*a");
        b = v;
#endif
    }
    if (showstartdone) std::cout<<"Done VM2a"<<std::endl;
}

template <class MM, class V> 
static void DoTestVM2a_Full(const MM& a, V& b, std::string label)
{
    DoTestVM2a_Basic(a,b,label);

#if (XTEST & 16)
    typedef typename MM::value_type Ta;
    typedef typename V::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Continuing on to VM2a Full"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }
    tmv::Matrix<Ta> m = a;
    tmv::Vector<T> v = b;

    RT eps = EPS * (a.colsize() + a.rowsize());
    RT eps2 = eps;
    if (!std::numeric_limits<RT>::is_integer) {
        eps *= Norm(m) * Norm(v);
        eps2 *= Norm(v*m) + Norm(m)*Norm(v);
    }

    if (CanMultVM(b,a)) {
        tmv::Vector<T> prod = v*m;
        tmv::Vector<T> c = prod;
        tmv::Vector<T> c0 = c;
        c0 /= T(2);
        RT x(5);
        T z; SetZ(z);
        c = c0;
        c += -b*a;
        Assert(Equal(c,(c0-prod),eps2),label+" c+=-b*a");
        c = c0;
        c += x*b*a;
        Assert(Equal(c,(c0+x*prod),x*eps2),label+" c+=x*b*a");
        c = c0;
        c += z*b*a;
        Assert(Equal(c,(c0+z*prod),x*eps2),label+" c+=z*b*a");
        c = c0;
        c = -b*a;
        Assert(Equal(c,(-prod),eps),label+" c=-b*a");
        c = -x*b*a;
        Assert(Equal(c,(-x*prod),x*eps),label+" c=-x*b*a");
        c = -z*b*a;
        Assert(Equal(c,(-z*prod),x*eps),label+" c=-z*b*a");
        c = c0;
        c -= b*a;
        Assert(Equal(c,(c0-prod),eps2),label+" c-=b*a");
        c = c0;
        c -= -b*a;
        Assert(Equal(c,(c0+prod),eps2),label+" c-=-b*a");
        c = c0;
        c += -x*b*a;
        Assert(Equal(c,(c0-x*prod),x*eps2),label+" c+=-x*b*a");
        c = c0;
        c -= x*b*a;
        Assert(Equal(c,(c0-x*prod),x*eps2),label+" c-=x*b*a");
        c = c0;
        c -= -x*b*a;
        Assert(Equal(c,(c0+x*prod),x*eps2),label+" c-=-x*b*a");
        c = c0;
        c += -z*b*a;
        Assert(Equal(c,(c0-z*prod),x*eps2),label+" c+=-z*b*a");
        c = c0;
        c -= z*b*a;
        Assert(Equal(c,(c0-z*prod),x*eps2),label+" c-=z*b*a");
        c = c0;
        c -= -z*b*a;
        Assert(Equal(c,(c0+z*prod),x*eps2),label+" c-=-z*b*a");
        c = c0;
    }

    if (CanMultVM(b,a,b)) {
        tmv::Vector<T> prod = v*m;
        b = v;
        RT x(5);
        T z; SetZ(z);
        b *= x*a;
        Assert(Equal(b,x*prod,x*eps),label+" b*=(x*a)");
        b = v;
        b *= z*a;
        Assert(Equal(b,z*prod,x*eps),label+" b*=(z*a)");
        b = v;
        b *= -a;
        Assert(Equal(b,(-prod),eps),label+" b*=-a");
        b = v;
#ifdef ALIASOK
        b = x*b*a;
        Assert(Equal(b,x*prod,x*eps),label+" b=x*b*a");
        b = v;
        b = z*b*a;
        Assert(Equal(b,z*prod,x*eps),label+" b=z*b*a");
        b = v;
        b += x*b*a;
        Assert(Equal(b,(v+x*prod),x*eps2),label+" b+=x*b*a");
        b = v;
        b += z*b*a;
        Assert(Equal(b,(v+z*prod),x*eps2),label+" b+=z*b*a");
        b = v;
        b = -b*a;
        Assert(Equal(b,(-prod),eps),label+" b=-b*a");
        b = v;
        b = -x*b*a;
        Assert(Equal(b,(-x*prod),x*eps),label+" b=-x*b*a");
        b = v;
        b = -z*b*a;
        Assert(Equal(b,(-z*prod),x*eps),label+" b=-z*b*a");
        b = v;
        b += -b*a;
        Assert(Equal(b,(v-prod),eps2),label+" b+=-b*a");
        b = v;
        b -= b*a;
        Assert(Equal(b,(v-prod),eps2),label+" b-=b*a");
        b = v;
        b -= -b*a;
        Assert(Equal(b,(v+prod),eps2),label+" b-=-b*a");
        b = v;
        b += -x*b*a;
        Assert(Equal(b,(v-x*prod),x*eps2),label+" b+=-x*b*a");
        b = v;
        b -= x*b*a;
        Assert(Equal(b,(v-x*prod),x*eps2),label+" b-=x*b*a");
        b = v;
        b -= -x*b*a;
        Assert(Equal(b,(v+x*prod),x*eps2),label+" b-=-x*b*a");
        b = v;
        b += -z*b*a;
        Assert(Equal(b,(v-z*prod),x*eps2),label+" b+=-z*b*a");
        b = v;
        b -= z*b*a;
        Assert(Equal(b,(v-z*prod),x*eps2),label+" b-=z*b*a");
        b = v;
        b -= -z*b*a;
        Assert(Equal(b,(v+z*prod),x*eps2),label+" b-=-z*b*a");
        b = v;
#endif
    }
    if (showstartdone) std::cout<<"Done VM2a_Full"<<std::endl;
#endif
}

template <class MM, class V> 
static void DoTestMV2R(const MM& a, V& b, std::string label)
{
    typename V::reverse_type br = b.reverse();
    DoTestMV2a_Basic(a,b,label);
    DoTestMV2a_Basic(a,br,label+" Rev");

#if (XTEST & 8)
    typename V::copy_type b0 = b;
    b.setZero();
    DoTestMV2a_Basic(a,b,label+" 1");
    CopyBackV(b0,b);

    b.cSubVector(0,b.size()/2).setZero();
    DoTestMV2a_Basic(a,b,label+" 2");
    CopyBackV(b0,b);

    b.cSubVector(b.size()/2,b.size()).setZero();
    DoTestMV2a_Basic(a,b,label+" 3");
    CopyBackV(b0,b);

    b.cSubVector(0,b.size()/4).setZero();
    b.cSubVector(3*b.size()/4,b.size()).setZero();
    DoTestMV2a_Basic(a,b,label+" 4");
    CopyBackV(b0,b);
#endif
}

template <class MM, class V> 
static void DoTestVM2R(const MM& a, V& b, std::string label)
{
    typename V::reverse_type br = b.reverse();
    DoTestVM2a_Basic(a,b,label);
    DoTestVM2a_Basic(a,br,label+" Rev");

#if (XTEST & 8)
    typename V::copy_type b0 = b;
    b.setZero();
    DoTestVM2a_Basic(a,b,label+" 1");
    CopyBackV(b0,b);

    b.cSubVector(0,b.size()/2).setZero();
    DoTestVM2a_Basic(a,b,label+" 2");
    CopyBackV(b0,b);

    b.cSubVector(b.size()/2,b.size()).setZero();
    DoTestVM2a_Basic(a,b,label+" 3");
    CopyBackV(b0,b);

    b.cSubVector(0,b.size()/4).setZero();
    b.cSubVector(3*b.size()/4,b.size()).setZero();
    DoTestVM2a_Basic(a,b,label+" 4");
    CopyBackV(b0,b);
#endif
}

template <class MM, class V> 
static void DoTestMV2C(const MM& a, V& b, std::string label)
{
    typename V::reverse_type br = b.reverse();
    DoTestMV2a_Full(a,b,label);
    DoTestMV2a_Basic(a,br,label+" Rev");

    DoTestMV2a_Basic(Conjugate(a),b,label+" Conj");
#if (XTEST & 8)
    DoTestMV2a_Basic(Conjugate(a),br,label+" Rev,Conj");

    typename V::copy_type b0 = b;
    b.setZero();
    DoTestMV2a_Basic(a,b,label+" 1");
    DoTestMV2a_Basic(Conjugate(a),b,label+" Conj1");
    CopyBackV(b0,b);

    b.cSubVector(0,b.size()/2).setZero();
    DoTestMV2a_Basic(a,b,label+" 2");
    DoTestMV2a_Basic(Conjugate(a),b,label+" Conj2");
    CopyBackV(b0,b);

    b.cSubVector(b.size()/2,b.size()).setZero();
    DoTestMV2a_Basic(a,b,label+" 3");
    DoTestMV2a_Basic(Conjugate(a),b,label+" Conj3");
    CopyBackV(b0,b);

    b.cSubVector(0,b.size()/4).setZero();
    b.cSubVector(3*b.size()/4,b.size()).setZero();
    DoTestMV2a_Basic(a,b,label+" 4");
    DoTestMV2a_Basic(Conjugate(a),b,label+" Conj4");
    CopyBackV(b0,b);
#endif
}

template <class MM, class V> 
static void DoTestVM2C(const MM& a, V& b, std::string label)
{
    typename V::reverse_type br = b.reverse();
    DoTestVM2a_Full(a,b,label);
    DoTestVM2a_Basic(a,br,label+" Rev");

    DoTestVM2a_Basic(Conjugate(a),b,label+" Conj");
#if (XTEST & 8)
    DoTestVM2a_Basic(Conjugate(a),br,label+" Rev,Conj");

    typename V::copy_type b0 = b;
    b.setZero();
    DoTestVM2a_Basic(a,b,label+" 1");
    DoTestVM2a_Basic(Conjugate(a),b,label+" Conj1");
    CopyBackV(b0,b);

    b.cSubVector(0,b.size()/2).setZero();
    DoTestVM2a_Basic(a,b,label+" 2");
    DoTestVM2a_Basic(Conjugate(a),b,label+" Conj2");
    CopyBackV(b0,b);

    b.cSubVector(b.size()/2,b.size()).setZero();
    DoTestVM2a_Basic(a,b,label+" 3");
    DoTestVM2a_Basic(Conjugate(a),b,label+" Conj3");
    CopyBackV(b0,b);

    b.cSubVector(0,b.size()/4).setZero();
    b.cSubVector(3*b.size()/4,b.size()).setZero();
    DoTestVM2a_Basic(a,b,label+" 4");
    DoTestVM2a_Basic(Conjugate(a),b,label+" Conj4");
    CopyBackV(b0,b);
#endif
}

template <class MM, class V1, class V2> 
static void DoTestMV3a_Basic(
    const MM& a, const V1& b, V2& c, std::string label)
{
    typedef typename MM::value_type Ta;
    typedef typename V1::value_type Tb;
    typedef typename V2::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start MV3a"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
        std::cout<<"c = "<<tmv::TMV_Text(c)<<std::endl;
    }

    tmv::Matrix<Ta> m = a;
    tmv::Matrix<Ta> mt = Transpose(a);
    tmv::Vector<Tb> v1 = b;
    tmv::Vector<T> v2 = c;

    RT eps = EPS * (a.colsize() + a.rowsize());
    RT eps2 = eps;
    if (!std::numeric_limits<RT>::is_integer) {
        eps *= Norm(m) * Norm(v1);
        eps2 *= Norm(m*v1) + Norm(m)*Norm(v1);
    }

#ifdef XXD
    if (XXDEBUG6) {
        std::cout<<"a = "<<TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"b = "<<TMV_Text(b)<<"  "<<b.step()<<"  "<<b<<std::endl;
        std::cout<<"c = "<<TMV_Text(c)<<"  "<<c.step()<<"  "<<c<<std::endl;
        std::cout<<"a*b = "<<m*v1<<std::endl;
    }
#endif

    if (CanMultMV(a,b,c)) {
        typename V2::copy_type c0 = c;
        c = a*b;
        v2 = m*v1;
        Assert(Equal(c,v2,eps),label+" c=a*b");
        CopyBackV(c0,c);
        c += a*b;
        v2 = c0 + m*v1;
        Assert(Equal(c,v2,eps2),label+" c+=a*b");
        CopyBackV(c0,c);
        RT x(5);
        T z; SetZ(z);
        c = x*a*b;
        v2 = x*m*v1;
        Assert(Equal(c,v2,x*eps),label+" c=x*a*b");
        CopyBackV(c0,c);
        c = z*a*b;
        v2 = z*m*v1;
        Assert(Equal(c,v2,x*eps),label+" c=z*a*b");
        CopyBackV(c0,c);

        c = b*Transpose(a);
        v2 = v1*mt;
        Assert(Equal(c,v2,eps),label+" c=b*at");
        CopyBackV(c0,c);
        c += b*Transpose(a);
        v2 = c0 + v1*mt;
        Assert(Equal(c,v2,eps2),label+" c+=b*at");
        CopyBackV(c0,c);
        c = x*b*Transpose(a);
        v2 = x*v1*mt;
        Assert(Equal(c,v2,x*eps),label+" c=x*b*at");
        CopyBackV(c0,c);
        c = z*b*Transpose(a);
        v2 = z*v1*mt;
        Assert(Equal(c,v2,x*eps),label+" c=z*b*at");
        CopyBackV(c0,c);
    }

    if (showstartdone) std::cout<<"Done MV3a"<<std::endl;
}

template <class MM, class V1, class V2> 
static void DoTestMV3a_Full(const MM& a, const V1& b, V2& c, std::string label)
{
    DoTestMV3a_Basic(a,b,c,label);

#if (XTEST & 16)
    typedef typename MM::value_type Ta;
    typedef typename V1::value_type Tb;
    typedef typename V2::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Continuing on to MV3a Full"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
        std::cout<<"c = "<<tmv::TMV_Text(c)<<std::endl;
    }
    tmv::Matrix<Ta> m = a;
    tmv::Matrix<Ta> mt = Transpose(a);
    tmv::Vector<Tb> v1 = b;
    tmv::Vector<T> v2 = c;

    RT eps = EPS * (a.colsize() + a.rowsize());
    RT eps2 = eps;
    if (!std::numeric_limits<RT>::is_integer) {
        eps *= Norm(m) * Norm(v1);
        eps2 *= Norm(m*v1) + Norm(m)*Norm(v1);
    }

    if (CanMultMV(a,b,c)) {
        typename V2::copy_type c0 = c;
        RT x(5);
        T z; SetZ(z);
        c += x*a*b;
        v2 = c0 + x*m*v1;
        Assert(Equal(c,v2,x*eps2),label+" c+=x*a*b");
        CopyBackV(c0,c);
        c += z*a*b;
        v2 = c0 + z*m*v1;
        Assert(Equal(c,v2,x*eps2),label+" c+=z*a*b");
        CopyBackV(c0,c);
        c = -a*b;
        v2 = -m*v1;
        Assert(Equal(c,v2,eps),label+" c=-a*b");
        CopyBackV(c0,c);
        c += -a*b;
        v2 = c0 - m*v1;
        Assert(Equal(c,v2,eps2),label+" c+=-a*b");
        CopyBackV(c0,c);
        c -= a*b;
        v2 = c0 - m*v1;
        Assert(Equal(c,v2,eps2),label+" c-=a*b");
        CopyBackV(c0,c);
        c -= -a*b;
        v2 = c0 + m*v1;
        Assert(Equal(c,v2,eps2),label+" c-=-a*b");
        CopyBackV(c0,c);
        c += -x*a*b;
        v2 = c0 - x*m*v1;
        Assert(Equal(c,v2,x*eps2),label+" c+=-x*a*b");
        CopyBackV(c0,c);
        c -= x*a*b;
        v2 = c0 - x*m*v1;
        Assert(Equal(c,v2,x*eps2),label+" c-=x*a*b");
        CopyBackV(c0,c);
        c -= -x*a*b;
        v2 = c0 + x*m*v1;
        Assert(Equal(c,v2,x*eps2),label+" c-=-x*a*b");
        CopyBackV(c0,c);
        c += -z*a*b;
        v2 = c0 - z*m*v1;
        Assert(Equal(c,v2,x*eps2),label+" c+=-z*a*b");
        CopyBackV(c0,c);
        c -= z*a*b;
        v2 = c0 - z*m*v1;
        Assert(Equal(c,v2,x*eps2),label+" c-=z*a*b");
        CopyBackV(c0,c);
        c -= -z*a*b;
        v2 = c0 + z*m*v1;
        Assert(Equal(c,v2,x*eps2),label+" c-=-z*a*b");
        CopyBackV(c0,c);

        c += x*b*Transpose(a);
        v2 = c0 + x*v1*mt;
        Assert(Equal(c,v2,x*eps2),label+" c+=x*b*at");
        CopyBackV(c0,c);
        c += z*b*Transpose(a);
        v2 = c0 + z*v1*mt;
        Assert(Equal(c,v2,x*eps2),label+" c+=z*b*at");
        CopyBackV(c0,c);
        c = -b*Transpose(a);
        v2 = -v1*mt;
        Assert(Equal(c,v2,eps),label+" c=-b*at");
        CopyBackV(c0,c);
        c += -b*Transpose(a);
        v2 = c0 - v1*mt;
        Assert(Equal(c,v2,eps2),label+" c+=-b*at");
        CopyBackV(c0,c);
        c -= b*Transpose(a);
        v2 = c0 - v1*mt;
        Assert(Equal(c,v2,eps2),label+" c-=b*at");
        CopyBackV(c0,c);
        c -= -b*Transpose(a);
        v2 = c0 + v1*mt;
        Assert(Equal(c,v2,eps2),label+" c-=-b*at");
        CopyBackV(c0,c);
        c += -x*b*Transpose(a);
        v2 = c0 - x*v1*mt;
        Assert(Equal(c,v2,x*eps2),label+" c+=-x*b*at");
        CopyBackV(c0,c);
        c -= x*b*Transpose(a);
        v2 = c0 - x*v1*mt;
        Assert(Equal(c,v2,x*eps2),label+" c-=x*b*at");
        CopyBackV(c0,c);
        c -= -x*b*Transpose(a);
        v2 = c0 + x*v1*mt;
        Assert(Equal(c,v2,x*eps2),label+" c-=-x*b*at");
        CopyBackV(c0,c);
        c += -z*b*Transpose(a);
        v2 = c0 - z*v1*mt;
        Assert(Equal(c,v2,x*eps2),label+" c+=-z*b*at");
        CopyBackV(c0,c);
        c -= z*b*Transpose(a);
        v2 = c0 - z*v1*mt;
        Assert(Equal(c,v2,x*eps2),label+" c-=z*b*at");
        CopyBackV(c0,c);
        c -= -z*b*Transpose(a);
        v2 = c0 + z*v1*mt;
        Assert(Equal(c,v2,x*eps2),label+" c-=-z*b*at");
        CopyBackV(c0,c);
    }

    if (showstartdone) std::cout<<"Done MV3a_Full"<<std::endl;
#endif
}

template <class MM, class V1, class V2> 
static void DoTestVM3a_Basic(
    const MM& a, const V1& b, V2& c, std::string label)
{
    typedef typename MM::value_type Ta;
    typedef typename V1::value_type Tb;
    typedef typename V2::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start VM3a"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
        std::cout<<"c = "<<tmv::TMV_Text(c)<<std::endl;
    }

    tmv::Matrix<Ta> m = a;
    tmv::Matrix<Ta> mt = Transpose(a);
    tmv::Vector<Tb> v1 = b;
    tmv::Vector<T> v2 = c;

    RT eps = EPS * (a.colsize() + a.rowsize());
    RT eps2 = eps;
    if (!std::numeric_limits<RT>::is_integer) {
        eps *= Norm(m) * Norm(v1);
        eps2 *= Norm(v1*m) + Norm(m)*Norm(v1);
    }

#ifdef XXD
    if (XXDEBUG6) {
        std::cout<<"a = "<<TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"b = "<<TMV_Text(b)<<"  "<<b.step()<<"  "<<b<<std::endl;
        std::cout<<"c = "<<TMV_Text(c)<<"  "<<c.step()<<"  "<<c<<std::endl;
        std::cout<<"b*a = "<<v1*m<<std::endl;
    }
#endif

    if (CanMultVM(b,a,c)) {
        typename V2::copy_type c0 = c;
        c = b*a;
        v2 = v1*m;
        Assert(Equal(c,v2,eps),label+" c=b*a");
        CopyBackV(c0,c);
        c += b*a;
        v2 = c0 + v1*m;
        Assert(Equal(c,v2,eps2),label+" c+=b*a");
        CopyBackV(c0,c);
        RT x(5);
        T z; SetZ(z);
        c = x*b*a;
        v2 = x*v1*m;
        Assert(Equal(c,v2,x*eps),label+" c=x*b*a");
        CopyBackV(c0,c);
        c = z*b*a;
        v2 = z*v1*m;
        Assert(Equal(c,v2,x*eps),label+" c=z*b*a");
        CopyBackV(c0,c);

        c = Transpose(a)*b;
        v2 = mt*v1;
        Assert(Equal(c,v2,eps),label+" c=at*b");
        CopyBackV(c0,c);
        c += Transpose(a)*b;
        v2 = c0 + mt*v1;
        Assert(Equal(c,v2,eps2),label+" c+=at*b");
        CopyBackV(c0,c);
        c = x*Transpose(a)*b;
        v2 = x*mt*v1;
        Assert(Equal(c,v2,x*eps),label+" c=x*at*b");
        CopyBackV(c0,c);
        c = z*Transpose(a)*b;
        v2 = z*mt*v1;
        Assert(Equal(c,v2,x*eps),label+" c=z*at*b");
        CopyBackV(c0,c);
    }

    if (showstartdone) std::cout<<"Done VM3a"<<std::endl;
}

template <class MM, class V1, class V2> 
static void DoTestVM3a_Full(const MM& a, const V1& b, V2& c, std::string label)
{
    DoTestVM3a_Basic(a,b,c,label);

#if (XTEST & 16)
    typedef typename MM::value_type Ta;
    typedef typename V1::value_type Tb;
    typedef typename V2::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    typedef typename tmv::Traits<T>::complex_type CT;
    if (showstartdone) {
        std::cout<<"Continuing on to VM3a Full"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
        std::cout<<"c = "<<tmv::TMV_Text(c)<<std::endl;
    }
    tmv::Matrix<Ta> m = a;
    tmv::Matrix<Ta> mt = Transpose(a);
    tmv::Vector<Tb> v1 = b;
    tmv::Vector<T> v2 = c;

    RT eps = EPS * (a.colsize() + a.rowsize());
    RT eps2 = eps;
    if (!std::numeric_limits<RT>::is_integer) {
        eps *= Norm(m) * Norm(v1);
        eps2 *= Norm(v1*m) + Norm(m)*Norm(v1);
    }

    if (CanMultVM(b,a,c)) {
        typename V2::copy_type c0 = c;
        RT x(5);
        CT z(3,4);
        c += x*b*a;
        v2 = c0 + x*v1*m;
        Assert(Equal(c,v2,x*eps2),label+" c+=x*b*a");
        CopyBackV(c0,c);
        c += z*b*a;
        v2 = c0 + z*v1*m;
        Assert(Equal(c,v2,x*eps2),label+" c+=z*b*a");
        CopyBackV(c0,c);
        c = -b*a;
        v2 = -v1*m;
        Assert(Equal(c,v2,eps),label+" c=-b*a");
        CopyBackV(c0,c);
        c += -b*a;
        v2 = c0 - v1*m;
        Assert(Equal(c,v2,eps2),label+" c+=-b*a");
        CopyBackV(c0,c);
        c -= b*a;
        v2 = c0 - v1*m;
        Assert(Equal(c,v2,eps2),label+" c-=b*a");
        CopyBackV(c0,c);
        c -= -b*a;
        v2 = c0 + v1*m;
        Assert(Equal(c,v2,eps2),label+" c-=-b*a");
        CopyBackV(c0,c);
        c += -x*b*a;
        v2 = c0 - x*v1*m;
        Assert(Equal(c,v2,x*eps2),label+" c+=-x*b*a");
        CopyBackV(c0,c);
        c -= x*b*a;
        v2 = c0 - x*v1*m;
        Assert(Equal(c,v2,x*eps2),label+" c-=x*b*a");
        CopyBackV(c0,c);
        c -= -x*b*a;
        v2 = c0 + x*v1*m;
        Assert(Equal(c,v2,x*eps2),label+" c-=-x*b*a");
        CopyBackV(c0,c);
        c += -z*b*a;
        v2 = c0 - z*v1*m;
        Assert(Equal(c,v2,x*eps2),label+" c+=-z*b*a");
        CopyBackV(c0,c);
        c -= z*b*a;
        v2 = c0 - z*v1*m;
        Assert(Equal(c,v2,x*eps2),label+" c-=z*b*a");
        CopyBackV(c0,c);
        c -= -z*b*a;
        v2 = c0 + z*v1*m;
        Assert(Equal(c,v2,x*eps2),label+" c-=-z*b*a");
        CopyBackV(c0,c);

        c += x*Transpose(a)*b;
        v2 = c0 + x*mt*v1;
        Assert(Equal(c,v2,x*eps2),label+" c+=x*at*b");
        CopyBackV(c0,c);
        c += z*Transpose(a)*b;
        v2 = c0 + z*mt*v1;
        Assert(Equal(c,v2,x*eps2),label+" c+=z*at*b");
        CopyBackV(c0,c);
        c = -Transpose(a)*b;
        v2 = -mt*v1;
        Assert(Equal(c,v2,eps),label+" c=-at*b");
        CopyBackV(c0,c);
        c += -Transpose(a)*b;
        v2 = c0 - mt*v1;
        Assert(Equal(c,v2,eps2),label+" c+=-at*b");
        CopyBackV(c0,c);
        c -= Transpose(a)*b;
        v2 = c0 - mt*v1;
        Assert(Equal(c,v2,eps2),label+" c-=at*b");
        CopyBackV(c0,c);
        c -= -Transpose(a)*b;
        v2 = c0 + mt*v1;
        Assert(Equal(c,v2,eps2),label+" c-=-at*b");
        CopyBackV(c0,c);
        c += -x*Transpose(a)*b;
        v2 = c0 - x*mt*v1;
        Assert(Equal(c,v2,x*eps2),label+" c+=-x*at*b");
        CopyBackV(c0,c);
        c -= x*Transpose(a)*b;
        v2 = c0 - x*mt*v1;
        Assert(Equal(c,v2,x*eps2),label+" c-=x*at*b");
        CopyBackV(c0,c);
        c -= -x*Transpose(a)*b;
        v2 = c0 + x*mt*v1;
        Assert(Equal(c,v2,x*eps2),label+" c-=-x*at*b");
        CopyBackV(c0,c);
        c += -z*Transpose(a)*b;
        v2 = c0 - z*mt*v1;
        Assert(Equal(c,v2,x*eps2),label+" c+=-z*at*b");
        CopyBackV(c0,c);
        c -= z*Transpose(a)*b;
        v2 = c0 - z*mt*v1;
        Assert(Equal(c,v2,x*eps2),label+" c-=z*at*b");
        CopyBackV(c0,c);
        c -= -z*Transpose(a)*b;
        v2 = c0 + z*mt*v1;
        Assert(Equal(c,v2,x*eps2),label+" c-=-z*at*b");
        CopyBackV(c0,c);
    }

    if (showstartdone) std::cout<<"Done VM3a"<<std::endl;
#endif
}

template <class MM, class V1, class V2> 
static void DoTestMV3R(const MM& a, const V1& b, V2& c, std::string label)
{
    typename V2::reverse_type cr = c.reverse();
    DoTestMV3a_Basic(a,b,c,label);
    DoTestMV3a_Basic(a,b.reverse(),c,label);
    DoTestMV3a_Basic(a,b,cr,label);
    DoTestMV3a_Basic(a,b.reverse(),cr,label);
}

template <class MM, class V1, class V2> 
static void DoTestVM3R(const MM& a, const V1& b, V2& c, std::string label)
{
    typename V2::reverse_type cr = c.reverse();
    DoTestVM3a_Basic(a,b,c,label);
    DoTestVM3a_Basic(a,b.reverse(),c,label);
    DoTestVM3a_Basic(a,b,cr,label);
    DoTestVM3a_Basic(a,b.reverse(),cr,label);
}

template <class MM, class V1, class V2> 
static void DoTestMV3C(const MM& a, const V1& b, V2& c, std::string label)
{
    DoTestMV3a_Full(a,b,c,label);
#if (XTEST & 1)
    typename V2::reverse_type cr = c.reverse();
    DoTestMV3a_Basic(a,b.reverse(),c,label);
    DoTestMV3a_Basic(a,b,cr,label);
    DoTestMV3a_Basic(a,b.reverse(),cr,label);
#endif

    DoTestMV3a_Basic(Conjugate(a),b,c,label+" Conj");
#if (XTEST & 1)
    DoTestMV3a_Basic(Conjugate(a),b.reverse(),c,label+" Conj");
    DoTestMV3a_Basic(Conjugate(a),b,cr,label+" Conj");
    DoTestMV3a_Basic(Conjugate(a),b.reverse(),cr,label+" Conj");
#endif
}

template <class MM, class V1, class V2> 
static void DoTestVM3C(const MM& a, const V1& b, V2& c, std::string label)
{
    DoTestVM3a_Full(a,b,c,label);
#if (XTEST & 1)
    typename V2::reverse_type cr = c.reverse();
    DoTestVM3a_Basic(a,b.reverse(),c,label);
    DoTestVM3a_Basic(a,b,cr,label);
    DoTestVM3a_Basic(a,b.reverse(),cr,label);
#endif

    DoTestVM3a_Basic(Conjugate(a),b,c,label+" Conj");
#if (XTEST & 1)
    DoTestVM3a_Basic(Conjugate(a),b.reverse(),c,label+" Conj");
    DoTestVM3a_Basic(Conjugate(a),b,cr,label+" Conj");
    DoTestVM3a_Basic(Conjugate(a),b.reverse(),cr,label+" Conj");
#endif
}

template <class M1, class M2> 
static void DoTestMM1a_Basic(const M1& a, const M2& b, std::string label)
{
    typedef typename M1::value_type T;
    typedef typename M2::value_type Tb;
    typedef typename tmv::Traits<T>::real_type RT;
    typedef typename tmv::Traits<T>::complex_type CT;
    if (showstartdone) {
        std::cout<<"Start MM1a"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }
#ifdef XXD
    if (XXDEBUG7) {
        std::cout<<"a = "<<tmv::TMV_Text(a)<<" = "<<a<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<" = "<<b<<std::endl;
    }
#endif
    const int M = a.colsize();
    const int N = a.rowsize();

    if (CanAdd(a,b)) {
        tmv::Matrix<T> m1 = a;
        tmv::Matrix<Tb> m2 = b;

        RT eps = EPS * (a.colsize() + a.rowsize());
        if (!std::numeric_limits<RT>::is_integer) eps *= Norm(m1) + Norm(m2);

        tmv::Matrix<CT> m3 = m1+m2;
#ifdef XXD
        if (XXDEBUG7) {
            std::cout<<"m1+m2 = "<<m3<<std::endl;
            std::cout<<"a+b = "<<(a+b)<<std::endl;
        }
#endif
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = m1(i,j) + m2(i,j);
        Assert(Equal(a+b,m3,eps),label+" a+b");
        Assert(Equal(m1+b,m3,eps),label+" m+b");
        Assert(Equal(a+m2,m3,eps),label+" a+m");

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = m1(i,j) - m2(i,j);
#ifdef XXD
        if (XXDEBUG7) {
            std::cout<<"m1-m2 = "<<m3<<std::endl;
            std::cout<<"a-b = "<<(a-b)<<std::endl;
        }
#endif
        Assert(Equal(a-b,m3,eps),label+" a-b");
        Assert(Equal(m1-b,m3,eps),label+" m-b");
        Assert(Equal(a-m2,m3,eps),label+" a-m");

        RT x(5);
        CT z(3,4);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = m1(i,j) + z*m2(i,j);
        Assert(Equal(a+z*b,m3,x*eps),label+" a+z*b");
        Assert(Equal(m1+z*b,m3,x*eps),label+" m+z*b");
        Assert(Equal(a+z*m2,m3,x*eps),label+" a+z*m");

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = m1(i,j) - z*m2(i,j);
        Assert(Equal(a-z*b,m3,x*eps),label+" a-z*b");
        Assert(Equal(m1-z*b,m3,x*eps),label+" m-z*b");
        Assert(Equal(a-z*m2,m3,x*eps),label+" a-z*m");
    }
#ifndef NOELEMMULT
    if (CanElemMultMM(a,b)) {
        tmv::Matrix<T> m1 = a;
        tmv::Matrix<Tb> m2 = b;

        RT eps = EPS * (a.colsize() + a.rowsize());
        if (!std::numeric_limits<RT>::is_integer) eps *= Norm(m1) * Norm(m2);

        tmv::Matrix<CT> m3 = m1;
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = m1(i,j) * m2(i,j);
#ifdef XXD
        if (XXDEBUG7) {
            std::cout<<"ElemProd(m1,m2) = "<<m3<<std::endl;
            std::cout<<"a+b = "<<(a+b)<<std::endl;
        }
#endif
        Assert(Equal(ElemProd(a,b),m3,eps),label+" ElemProd(a,b)");
    }
#endif
 
    if (showstartdone) std::cout<<"Done MM1a"<<std::endl;
}

template <class M1, class M2> 
static void DoTestMM1a_Full(const M1& a, const M2& b, std::string label)
{
    DoTestMM1a_Basic(a,b,label);

#if (XTEST & 16)
    typedef typename M1::value_type T;
    typedef typename M2::value_type Tb;
    typedef typename tmv::Traits<T>::real_type RT;
    typedef typename tmv::Traits<T>::complex_type CT;
    if (showstartdone) {
        std::cout<<"Continuing on to MM1a Full"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }
    if (CanAdd(a,b)) {
        const int M = a.colsize();
        const int N = a.rowsize();

        tmv::Matrix<T> m1 = a;
        tmv::Matrix<Tb> m2 = b;

        RT eps = EPS * (a.colsize() + a.rowsize());
        if (!std::numeric_limits<RT>::is_integer) eps *= Norm(m1) * Norm(m2);

        RT x(5);
        CT z(3,4);
        tmv::Matrix<CT> m3 = m1;
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = m1(i,j) + x*m2(i,j);
        Assert(Equal(a+x*b,m3,x*eps),label+" a+x*b");
        Assert(Equal(m1+x*b,m3,x*eps),label+" m+x*b");
        Assert(Equal(a+x*m2,m3,x*eps),label+" a+x*m");

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = m1(i,j) - x*m2(i,j);
        Assert(Equal(a-x*b,m3,x*eps),label+" a-x*b");
        Assert(Equal(m1-x*b,m3,x*eps),label+" m-x*b");
        Assert(Equal(a-x*m2,m3,x*eps),label+" a-x*m");

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = x*m1(i,j) + m2(i,j);
        Assert(Equal(x*a+b,m3,x*eps),label+" x*a+b");
        Assert(Equal(x*m1+b,m3,x*eps),label+" x*m+b");
        Assert(Equal(x*a+m2,m3,x*eps),label+" x*a+m");

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = x*m1(i,j) - m2(i,j);
        Assert(Equal(x*a-b,m3,x*eps),label+" x*a-b");
        Assert(Equal(x*m1-b,m3,x*eps),label+" x*m-b");
        Assert(Equal(x*a-m2,m3,x*eps),label+" x*a-m");

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = x*m1(i,j) + x*m2(i,j);
        Assert(Equal(x*a+x*b,m3,x*eps),label+" x*a+x*b");
        Assert(Equal(x*m1+x*b,m3,x*eps),label+" x*m+x*b");
        Assert(Equal(x*a+x*m2,m3,x*eps),label+" x*a+x*m");

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = x*m1(i,j) - x*m2(i,j);
        Assert(Equal(x*a-x*b,m3,x*eps),label+" x*a-x*b");
        Assert(Equal(x*m1-x*b,m3,x*eps),label+" x*m-x*b");
        Assert(Equal(x*a-x*m2,m3,x*eps),label+" x*a-x*m");

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = x*m1(i,j) + z*m2(i,j);
        Assert(Equal(x*a+z*b,m3,x*eps),label+" x*a+z*b");
        Assert(Equal(x*m1+z*b,m3,x*eps),label+" x*m+z*b");
        Assert(Equal(x*a+z*m2,m3,x*eps),label+" x*a+z*m");

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = x*m1(i,j) - z*m2(i,j);
        Assert(Equal(x*a-z*b,m3,x*eps),label+" x*a-z*b");
        Assert(Equal(x*m1-z*b,m3,x*eps),label+" x*m-z*b");
        Assert(Equal(x*a-z*m2,m3,x*eps),label+" x*a-z*m");

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = z*m1(i,j) + m2(i,j);
        Assert(Equal(z*a+b,m3,x*eps),label+" z*a+b");
        Assert(Equal(z*m1+b,m3,x*eps),label+" z*m+b");
        Assert(Equal(z*a+m2,m3,x*eps),label+" z*a+m");

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = z*m1(i,j) - m2(i,j);
        Assert(Equal(z*a-b,m3,x*eps),label+" z*a-b");
        Assert(Equal(z*m1-b,m3,x*eps),label+" z*m-b");
        Assert(Equal(z*a-m2,m3,x*eps),label+" z*a-m");

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = z*m1(i,j) + x*m2(i,j);
        Assert(Equal(z*a+x*b,m3,x*eps),label+" z*a+x*b");
        Assert(Equal(z*m1+x*b,m3,x*eps),label+" z*m+x*b");
        Assert(Equal(z*a+x*m2,m3,x*eps),label+" z*a+x*m");

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = z*m1(i,j) - x*m2(i,j);
        Assert(Equal(z*a-x*b,m3,x*eps),label+" z*a-x*b");
        Assert(Equal(z*m1-x*b,m3,x*eps),label+" z*m-x*b");
        Assert(Equal(z*a-x*m2,m3,x*eps),label+" z*a-x*m");

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = z*m1(i,j) + z*m2(i,j);
        Assert(Equal(z*a+z*b,m3,x*eps),label+" z*a+z*b");
        Assert(Equal(z*m1+z*b,m3,x*eps),label+" z*m+z*b");
        Assert(Equal(z*a+z*m2,m3,x*eps),label+" z*a+z*m");

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = z*m1(i,j) - z*m2(i,j);
        Assert(Equal(z*a-z*b,m3,x*eps),label+" z*a-z*b");
        Assert(Equal(z*m1-z*b,m3,x*eps),label+" z*m-z*b");
        Assert(Equal(z*a-z*m2,m3,x*eps),label+" z*a-z*m");
    }
    if (showstartdone) std::cout<<"Done MM1a"<<std::endl;
#endif
}

template <class M1, class M2> 
static void DoTestMM1RR(const M1& a, const M2& b, std::string label)
{
    DoTestMM1a_Basic(a,b,label);

#if (XTEST & 2)
    DoTestMM1a_Basic(Transpose(b),Transpose(a), label+" TransB TransA");
#endif
}

template <class M1, class M2> 
static void DoTestMM1RC(const M1& a, const M2& b, std::string label)
{
    DoTestMM1a_Basic(a,b,label);

#if (XTEST & 2)
    DoTestMM1a_Basic(a,Conjugate(b), label+" ConjB");
    DoTestMM1a_Basic(Transpose(b),Transpose(a), label+" TransB TransA");
    DoTestMM1a_Basic(Adjoint(b),Transpose(a), label+" AdjB TransA");
#endif
}

template <class M1, class M2> 
static void DoTestMM1CR(const M1& a, const M2& b, std::string label)
{
    DoTestMM1a_Basic(a,b,label);

#if (XTEST & 2)
    DoTestMM1a_Basic(Conjugate(a),b, label+" ConjA");
    DoTestMM1a_Basic(Transpose(b),Transpose(a), label+" TransB TransA");
    DoTestMM1a_Basic(Transpose(b),Adjoint(a), label+" TransB AdjA");
#endif
}

template <class M1, class M2> 
static void DoTestMM1CC(const M1& a, const M2& b, std::string label)
{
    DoTestMM1a_Full(a,b,label);

#if (XTEST & 2)
    DoTestMM1a_Basic(Conjugate(a),b, label+" ConjA");
    DoTestMM1a_Basic(a,Conjugate(b), label+" ConjB");
    DoTestMM1a_Basic(Conjugate(a),Conjugate(b), label+" ConjA ConjB");
    DoTestMM1a_Basic(Transpose(b), Transpose(a),label+" TransB TransA");
    DoTestMM1a_Basic(Transpose(b),Adjoint(a), label+" TransB AdjA");
    DoTestMM1a_Basic(Adjoint(b),Transpose(a), label+" AdjB TransA");
    DoTestMM1a_Basic(Adjoint(b),Adjoint(a), label+" AdjA ConjB");
#endif
} 

template <class M1, class M2> 
static void DoTestMM2a_Basic(M1& a, const M2& b, std::string label)
{
    typedef typename M1::value_type T;
    typedef typename M2::value_type Tb;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start MM2a"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<" "<<a<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<" "<<b<<std::endl;
    }

#ifdef XXD
    if (XXDEBUG8) {
        std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
    }
#endif
    const tmv::Matrix<T> m1 = a;
    const tmv::Matrix<Tb> m2 = b;
#ifdef XXD
    if (XXDEBUG8) {
        std::cout<<"m1 = "<<tmv::TMV_Text(m1)<<"  "<<m1<<std::endl;
        std::cout<<"m2 = "<<tmv::TMV_Text(m2)<<"  "<<m2<<std::endl;
    }
#endif

    typename M1::copy_type a0 = a;
    tmv::Matrix<T> m3 = a;

    RT eps = EPS * (a.colsize() + a.rowsize());
    if (!std::numeric_limits<RT>::is_integer) 
        eps *= Norm(m1) * Norm(m2);

#ifdef XXD
    if (XXDEBUG8) {
        std::cout<<"eps = "<<eps<<std::endl;
    }
#endif

#ifndef NOADDEQ
    if (CanAddEq(a,b)) {
        const int M = a.colsize();
        const int N = a.rowsize();

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = m1(i,j) + m2(i,j);
        a += b;
#ifdef XXD
        if (XXDEBUG8) {
            std::cout<<"a += b = "<<a<<std::endl;
            std::cout<<"m3 = "<<m3<<std::endl;
            std::cout<<"diff = "<<(a-m3)<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(a-m3)<<std::endl;
            std::cout<<"cf. eps = "<<eps<<std::endl;
        }
#endif
        Assert(Equal(a,m3,eps),label+" a += b");
        CopyBackM(a0,a);
#ifdef ALIASOK
        a = a+b;
        if (XXDEBUG8) {
            std::cout<<"a.ptr = "<<a.cptr()<<std::endl;
            std::cout<<"b.ptr = "<<b.cptr()<<std::endl;
            std::cout<<"a = a+b = "<<a<<std::endl;
            std::cout<<"m3 = "<<m3<<std::endl;
            std::cout<<"diff = "<<(a-m3)<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(a-m3)<<std::endl;
            std::cout<<"cf. eps = "<<eps<<std::endl;
        }
        Assert(Equal(a,m3,eps),label+" a = a+b");
        CopyBackM(a0,a);
#endif
    }
#endif // NOADDEQ

#ifndef NOELEMMULT
#ifdef ALIASOK
    if (CanElemMultMM(a,b)) {
        const int M = a.colsize();
        const int N = a.rowsize();

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = m1(i,j) * m2(i,j);
        a = ElemProd(a,b);
#ifdef XXD
        if (XXDEBUG8) {
            std::cout<<"a = ElemProd(a,b) = "<<a<<std::endl;
            std::cout<<"m3 = "<<m3<<std::endl;
            std::cout<<"diff = "<<(a-m3)<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(a-m3)<<std::endl;
            std::cout<<"cf. eps = "<<eps<<std::endl;
        }
#endif
        Assert(Equal(a,m3,eps),label+" a = ElemProd(a,b)");
        CopyBackM(a0,a);
        a = ElemProd(b,a);
        Assert(Equal(a,m3,eps),label+" a = ElemProd(b,a)");
        CopyBackM(a0,a);
    }

    if (CanAddElemMultMM(a,b,a)) {
        const int M = a.colsize();
        const int N = a.rowsize();

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = m1(i,j) + m1(i,j) * m2(i,j);
        a += ElemProd(a,b);
        Assert(Equal(a,m3,eps),label+" a += ElemProd(a,b)");
        CopyBackM(a0,a);
        a += ElemProd(b,a);
        Assert(Equal(a,m3,eps),label+" a += ElemProd(b,a)");
        CopyBackM(a0,a);

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = m1(i,j) - m1(i,j) * m2(i,j);
        a -= ElemProd(a,b);
        Assert(Equal(a,m3,eps),label+" a -= ElemProd(a,b)");
        CopyBackM(a0,a);
        a -= ElemProd(b,a);
        Assert(Equal(a,m3,eps),label+" a -= ElemProd(b,a)");
        CopyBackM(a0,a);

        RT x = 5;
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = x * m1(i,j) * m2(i,j);
        a = x*ElemProd(a,b);
        Assert(Equal(a,m3,eps),label+" a = x*ElemProd(a,b)");
        CopyBackM(a0,a);
        a = x*ElemProd(b,a);
        Assert(Equal(a,m3,eps),label+" a = x*ElemProd(b,a)");
        CopyBackM(a0,a);
        a = ElemProd(x*a,b);
        Assert(Equal(a,m3,eps),label+" a = ElemProd(x*a,b)");
        CopyBackM(a0,a);
        a = ElemProd(b,x*a);
        Assert(Equal(a,m3,eps),label+" a = ElemProd(b,x*a)");
        CopyBackM(a0,a);
        a = ElemProd(a,x*b);
        Assert(Equal(a,m3,eps),label+" a = ElemProd(a,x*b)");
        CopyBackM(a0,a);
        a = ElemProd(x*b,a);
        Assert(Equal(a,m3,eps),label+" a = ElemProd(x*b,a)");
        CopyBackM(a0,a);

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = m1(i,j) + x * m1(i,j) * m2(i,j);
        a += x*ElemProd(a,b);
        Assert(Equal(a,m3,eps),label+" a += x*ElemProd(a,b)");
        CopyBackM(a0,a);
        a += x*ElemProd(b,a);
        Assert(Equal(a,m3,eps),label+" a += x*ElemProd(b,a)");
        CopyBackM(a0,a);
        a += ElemProd(x*a,b);
        Assert(Equal(a,m3,eps),label+" a += ElemProd(x*a,b)");
        CopyBackM(a0,a);
        a += ElemProd(b,x*a);
        Assert(Equal(a,m3,eps),label+" a += ElemProd(b,x*a)");
        CopyBackM(a0,a);
        a += ElemProd(a,x*b);
        Assert(Equal(a,m3,eps),label+" a += ElemProd(a,x*b)");
        CopyBackM(a0,a);
        a += ElemProd(x*b,a);
        Assert(Equal(a,m3,eps),label+" a += ElemProd(x*b,a)");
        CopyBackM(a0,a);

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m3(i,j) = m1(i,j) - x * m1(i,j) * m2(i,j);
        a -= x*ElemProd(a,b);
        Assert(Equal(a,m3,eps),label+" a -= x*ElemProd(a,b)");
        CopyBackM(a0,a);
        a -= x*ElemProd(b,a);
        Assert(Equal(a,m3,eps),label+" a -= x*ElemProd(b,a)");
        CopyBackM(a0,a);
        a -= ElemProd(x*a,b);
        Assert(Equal(a,m3,eps),label+" a -= ElemProd(x*a,b)");
        CopyBackM(a0,a);
        a -= ElemProd(b,x*a);
        Assert(Equal(a,m3,eps),label+" a -= ElemProd(b,x*a)");
        CopyBackM(a0,a);
        a -= ElemProd(a,x*b);
        Assert(Equal(a,m3,eps),label+" a -= ElemProd(a,x*b)");
        CopyBackM(a0,a);
        a -= ElemProd(x*b,a);
        Assert(Equal(a,m3,eps),label+" a -= ElemProd(x*b,a)");
        CopyBackM(a0,a);
    }
#endif
#endif

    if (showstartdone) std::cout<<"Done MM2a"<<std::endl;
}

template <class M1, class M2> 
static void DoTestMM2a_Full(M1& a, const M2& b, std::string label)
{
    DoTestMM2a_Basic(a,b,label);

#if (XTEST & 16)
    typedef typename M1::value_type T;
    typedef typename M2::value_type Tb;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Continuing on to MM2a Full"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }

#ifndef NOADDEQ
    const tmv::Matrix<T> m1 = a;
    const tmv::Matrix<Tb> m2 = b;

    RT eps = EPS * (a.colsize() + a.rowsize());
    if (!std::numeric_limits<RT>::is_integer) 
        eps *= Norm(m1) * Norm(m2);

    if (CanAddEq(a,b)) {
        typename M1::copy_type a0 = a;
        tmv::Matrix<T> m3 = a;
        a += -b;
        m3 = m1-m2;
        Assert(Equal(a,m3,eps),label+" a += -b");
        CopyBackM(a0,a);
        a -= b;
        m3 = m1-m2;
        Assert(Equal(a,m3,eps),label+" a -= b");
        CopyBackM(a0,a);
        a -= -b;
        m3 = m1+m2;
        Assert(Equal(a,m3,eps),label+" a -= -b");
        CopyBackM(a0,a);
#ifdef ALIASOK
        a = a-b;
        m3 = m1-m2;
        Assert(Equal(a,m3,eps),label+" a = a-b");
        CopyBackM(a0,a);
        a = b+a; 
        m3 = m2+m1;
        Assert(Equal(a,m3,eps),label+" a = b+a");
        CopyBackM(a0,a);
        if (XXDEBUG8) {
            std::cout<<"a = "<<a<<std::endl;
            std::cout<<"b = "<<b<<std::endl;
        }
        a = b-a;
        m3 = m2-m1;
        if (XXDEBUG8) {
            std::cout<<"a = b-a = "<<a<<std::endl;
            std::cout<<"m3 = m2-m1 = "<<m3<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(a-m3)<<std::endl;
        }
        Assert(Equal(a,m3,eps),label+" a = b-a");
        CopyBackM(a0,a);
#endif
    }
#endif // NOADDEQ

    if (showstartdone) std::cout<<"Done MM2a_Full"<<std::endl;
#endif
}

template <class M1, class M2> 
static void DoTestMM2RR(M1& a, const M2& b, std::string label)
{
    DoTestMM2a_Basic(a,b,label);
}

template <class M1, class M2> 
static void DoTestMM2RC(M1& a, const M2& b, std::string label)
{
    DoTestMM2a_Basic(a,b,label);
#if (XTEST & 2)
    DoTestMM2a_Basic(a,Conjugate(b),label+" ConjB");
#endif
}

template <class M1, class M2> 
static void DoTestMM2CR(M1& a, const M2& b, std::string label)
{
    typename M1::conjugate_type ac = a.conjugate();
    DoTestMM2a_Basic(a,b,label);
    DoTestMM2a_Basic(ac,b,label+" ConjA");
}

template <class M1, class M2> 
static void DoTestMM2CC(M1& a, const M2& b, std::string label)
{
    typename M1::conjugate_type ac = a.conjugate();
    DoTestMM2a_Full(a,b,label);
    DoTestMM2a_Basic(ac,b,label+" ConjA");

#if (XTEST & 2)
    DoTestMM2a_Basic(a,Conjugate(b),label+" ConjB");
    DoTestMM2a_Basic(ac,Conjugate(b),label+" ConjA ConjB");
#endif
}

template <class M1, class M2> 
static void DoTestMM3a_Basic(const M1& a, const M2& b, std::string label)
{
    typedef typename M1::value_type T;
    typedef typename M2::value_type Tb;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start MM3a"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }

#ifdef XXD
    if (XXDEBUG7) {
        std::cout<<"a = "<<a<<std::endl;
        std::cout<<"b = "<<b<<std::endl;
    }
#endif

    if (CanMultMM(a,b)) {
        tmv::Matrix<T> m1 = a;
        tmv::Matrix<Tb> m2 = b;
        tmv::Matrix<ProductType(T,Tb)> mm = m1*m2;

        RT eps = EPS * (a.colsize() + a.rowsize());
        if (!std::numeric_limits<RT>::is_integer) 
            eps *= Norm(m1) * Norm(m2);

#ifdef XXD
        if (XXDEBUG7) {
            std::cout<<"CanMult("<<tmv::TMV_Text(a)<<","<<tmv::TMV_Text(b)<<")\n";
            std::cout<<"m1*m2 = "<<mm<<std::endl;
            std::cout<<"m1*b = "<<(m1*b)<<std::endl;
            std::cout<<"a*m2 = "<<(a*m2)<<std::endl;
            std::cout<<"a*b = "<<(a*b)<<std::endl;
        }
#endif
        Assert(Equal(m1*b,mm,eps),label+" m*b");
        Assert(Equal(a*m2,mm,eps),label+" a*m");
        Assert(Equal(a*b,mm,eps),label+" a*b");
    }

    if (showstartdone) std::cout<<"Done MM3a"<<std::endl;
}

template <class M1, class M2> 
static void DoTestMM3a_Full(const M1& a, const M2& b, std::string label)
{
    DoTestMM3a_Basic(a,b,label);

#if (XTEST & 16)
    typedef typename M1::value_type T;
    typedef typename M2::value_type Tb;
    typedef typename tmv::Traits<T>::real_type RT;
    typedef typename tmv::Traits2<T,Tb>::type PT;
    if (showstartdone) {
        std::cout<<"Continuing on to MM3a Full"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }
    if (CanMultMM(a,b)) {
        tmv::Matrix<T> m1 = a;
        tmv::Matrix<Tb> m2 = b;
        tmv::Matrix<PT> mm = m1*m2;

        RT eps = EPS * (a.colsize() + a.rowsize());
        if (!std::numeric_limits<RT>::is_integer) eps *= Norm(m1) * Norm(m2);

        RT x(5);
        std::complex<RT> z(3,4);

        Assert(Equal(((x*a)*b),(x*mm),x*eps),label+" (x*a)*b");
        Assert(Equal((x*(a*b)),(x*mm),x*eps),label+" x*(a*b)");
        Assert(Equal((a*(x*b)),(x*mm),x*eps),label+" a*(x*b)");
        Assert(Equal(((z*a)*b),(z*mm),x*eps),label+" (z*a)*b");
        Assert(Equal((z*(a*b)),(z*mm),x*eps),label+" z*(a*b)");
        Assert(Equal((a*(z*b)),(z*mm),x*eps),label+" a*(z*b)");
    }
#endif
    if (showstartdone) std::cout<<"Done MM3a_Full"<<std::endl;
}

template <class M1, class M2> 
static void DoTestMM3RR(const M1& a, const M2& b, std::string label)
{
    DoTestMM3a_Basic(a,b,label);

#if (XTEST & 2)
    DoTestMM3a_Basic(Transpose(b),Transpose(a), label+" TransB TransA");
#endif
}

template <class M1, class M2> 
static void DoTestMM3RC(const M1& a, const M2& b, std::string label)
{
    DoTestMM3a_Basic(a,b,label);

#if (XTEST & 2)
    DoTestMM3a_Basic(a,Conjugate(b),label+" ConjB");
    DoTestMM3a_Basic(Transpose(b),Transpose(a), label+" TransB TransA");
    DoTestMM3a_Basic(Adjoint(b),Transpose(a), label+" AdjB TransA");
#endif
}

template <class M1, class M2> 
static void DoTestMM3CR(const M1& a, const M2& b, std::string label)
{
    DoTestMM3a_Basic(a,b,label);

#if (XTEST & 2)
    DoTestMM3a_Basic(Conjugate(a),b,label+" ConjA");
    DoTestMM3a_Basic(Transpose(b),Transpose(a), label+" TransB TransA");
    DoTestMM3a_Basic(Transpose(b),Adjoint(a), label+" TransB AdjA");
#endif
}

template <class M1, class M2> 
static void DoTestMM3CC(const M1& a, const M2& b, std::string label)
{
    DoTestMM3a_Full(a,b,label);

#if (XTEST & 2)
    DoTestMM3a_Basic(Conjugate(a),b,label+" ConjA");
    DoTestMM3a_Basic(a,Conjugate(b),label+" ConjB");
    DoTestMM3a_Basic(Conjugate(a),Conjugate(b), label+" ConjA ConjB");
    DoTestMM3a_Basic(Transpose(b),Transpose(a), label+" TransB TransA");
    DoTestMM3a_Basic(Transpose(b),Adjoint(a), label+" TransB AdjA");
    DoTestMM3a_Basic(Adjoint(b),Transpose(a), label+" AdjB TransA");
    DoTestMM3a_Basic(Adjoint(b),Adjoint(a), label+" AdjA ConjB");
#endif
}

template <class M1, class M2> 
static void DoTestMM4a_Basic(M1& a, const M2& b, std::string label)
{
    typedef typename M1::value_type T;
    typedef typename M2::value_type Tb;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start MM4a"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }

    if (CanMultMM(a,b)) {
#ifdef XXD
        if (XXDEBUG8) {
            std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
            std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
        }
#endif
        const tmv::Matrix<T> m1 = a;
        const tmv::Matrix<Tb> m2 = b;
#ifdef XXD
        if (XXDEBUG8) {
            std::cout<<"m1 = "<<tmv::TMV_Text(m1)<<"  "<<m1<<std::endl;
            std::cout<<"m2 = "<<tmv::TMV_Text(m2)<<"  "<<m2<<std::endl;
        }
#endif

        RT eps = EPS * (a.colsize() + a.rowsize());
        RT eps2 = eps;
        if (!std::numeric_limits<RT>::is_integer) {
            eps *= Norm(m1) * Norm(m2);
            eps2 *= Norm(m1*m2) + Norm(m1) * Norm(m2);
        }

        tmv::Matrix<ProductType(T,Tb)> mm = m1*m2;
        tmv::Matrix<ProductType(T,Tb)> m3 = mm;
        tmv::Matrix<ProductType(T,Tb)> m4 = mm;
        tmv::Matrix<ProductType(T,Tb)> m0 = mm;
        m0 /= T(2);
#ifdef XXD
        if (XXDEBUG8) {
            std::cout<<"CanMult("<<tmv::TMV_Text(a)<<","<<tmv::TMV_Text(b)<<"\n";
            std::cout<<"m0 = "<<m0<<std::endl;
            std::cout<<"m1 * m2 = "<<mm<<std::endl;
            std::cout<<"a * b = "<<(m3 = a*b)<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(m3 - mm)<<std::endl;
        }
#endif
        RT x(5);
        T z; SetZ(z);
        m3 = a*b;
        m4 = mm;
        Assert(Equal(m3,m4,eps),label+" m = a*b");
        m3 = m0;
        m3 += a*b;
        m4 = m0 + mm;
#ifdef XXD
        if (XXDEBUG8) {
            std::cout<<"m0 + m1 * m2 = "<<m4<<std::endl;
            std::cout<<"m3 += a * b = "<<m3<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(m3 - m4)<<std::endl;
        }
#endif
        Assert(Equal(m3,m4,eps2),label+" m += a*b");
        m3 = m0;
        m3 = x*a*b;
        m4 = x*mm;
        Assert(Equal(m3,m4,x*eps2),label+" m = x*a*b");
        m3 = m0;
        m3 = z*a*b;
        m4 = z*mm;
        Assert(Equal(m3,m4,x*eps2),label+" m = z*a*b");
        m3 = m0;

#ifndef NOMULTEQ
        if (CanMultMM(a,b,a)) {
            typename M1::copy_type a0 = a;
#ifdef XXD
            if (XXDEBUG8) {
                std::cout<<"CanMult("<<tmv::TMV_Text(a)<<","<<tmv::TMV_Text(b)
                    <<","<<tmv::TMV_Text(a)<<")\n";
            }
#endif
            a *= b;
            m4 = mm;
#ifdef XXD
            if (XXDEBUG8) {
                std::cout<<"a *= b = "<<a<<std::endl;
                std::cout<<"m4 = "<<m4<<std::endl;
            }
#endif
            Assert(Equal(a,m4,eps),label+" a *= b");
            CopyBackM(a0,a);
#ifdef ALIASOK
            a = a*b;
            m4 = mm;
            Assert(Equal(a,m4,eps),label+" a = a*b");
            CopyBackM(a0,a);
#endif
        }
#endif // NOMULTEQ
    }

    if (showstartdone) std::cout<<"Done MM4a"<<std::endl;
}

template <class M1, class M2> 
static void DoTestMM4a_Full(M1& a, const M2& b, std::string label)
{
    DoTestMM4a_Basic(a,b,label);

#if (XTEST & 16)
    typedef typename M1::value_type T;
    typedef typename M2::value_type Tb;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Continuing on to MM4a Full"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
    }
    if (CanMultMM(a,b)) {
        const tmv::Matrix<T> m1 = a;
        const tmv::Matrix<Tb> m2 = b;

        RT eps = EPS * (a.colsize() + a.rowsize() + b.rowsize());
        RT eps2 = eps;
        if (!std::numeric_limits<RT>::is_integer) {
            eps *= Norm(m1) * Norm(m2);
            eps2 *= Norm(m1*m2) + Norm(m1) * Norm(m2);
        }

        tmv::Matrix<ProductType(T,Tb)> mm = m1*m2;
        tmv::Matrix<ProductType(T,Tb)> m3 = mm;
        tmv::Matrix<ProductType(T,Tb)> m4 = mm;
        tmv::Matrix<ProductType(T,Tb)> m0 = mm;
        m0 /= T(2);
        RT x(5);
        T z; SetZ(z);
        m3 = m0;
        m3 += x*a*b;
        m4 = m0 + x*mm;
        Assert(Equal(m3,m4,x*eps2),label+" m += x*a*b");
        m3 = m0;
        m3 += z*a*b;
        m4 = m0 + z*mm;
        Assert(Equal(m3,m4,x*eps2),label+" m += z*a*b");
        m3 = m0;
        m3 = -a*b;
        m4 = -mm;
        Assert(Equal(m3,m4,eps),label+" m = -a*b");
        m3 = m0;
        m3 += -a*b;
        m4 = m0 - mm;
        Assert(Equal(m3,m4,eps2),label+" m += -a*b");
        m3 = m0;
        m3 -= a*b;
        m4 = m0 - mm;
        Assert(Equal(m3,m4,eps2),label+" m -= a*b");
        m3 = m0;
        m3 -= -a*b;
        m4 = m0 + mm;
        Assert(Equal(m3,m4,eps2),label+" m -= -a*b");
        m3 = m0;
        m3 = -x*a*b;
        m4 = -x*mm;
        Assert(Equal(m3,m4,x*eps2),label+" m = -x*a*b");
        m3 = m0;
        m3 += -x*a*b;
        m4 = m0 - x*mm;
        Assert(Equal(m3,m4,x*eps2),label+" m += -x*a*b");
        m3 = m0;
        m3 -= x*a*b;
        m4 = m0 - x*mm;
        Assert(Equal(m3,m4,x*eps2),label+" m -= x*a*b");
        m3 = m0;
        m3 -= -x*a*b;
        m4 = m0 + x*mm;
        Assert(Equal(m3,m4,x*eps2),label+" m -= -x*a*b");
        m3 = m0;
        m3 = -z*a*b;
        m4 = -z*mm;
        Assert(Equal(m3,m4,x*eps2),label+" m = -z*a*b");
        m3 = m0;
        m3 += -z*a*b;
        m4 = m0 - z*mm;
        Assert(Equal(m3,m4,x*eps2),label+" m += -z*a*b");
        m3 = m0;
        m3 -= z*a*b;
        m4 = m0 - z*mm;
        Assert(Equal(m3,m4,x*eps2),label+" m -= z*a*b");
        m3 = m0;
        m3 -= -z*a*b;
        m4 = m0 + z*mm;
        Assert(Equal(m3,m4,x*eps2),label+" m -= -z*a*b");
        m3 = m0;

#ifndef BASIC_MULTMM_ONLY
#ifndef NOMULTEQ
        if (CanMultXMM(a,b,a)) {
            typename M1::copy_type a0 = a;
#ifdef XXD
            if (XXDEBUG8) {
                std::cout<<"CanMultXM("<<tmv::TMV_Text(a)<<","<<tmv::TMV_Text(b)
                    <<","<<tmv::TMV_Text(a)<<")\n";
            }
#endif
            m0 = a0;
            m4 = m0;
            a *= -b;
            m4 = -mm;
            Assert(Equal(a,m4,eps),label+" a *= -b");
            CopyBackM(a0,a);
            a *= x*b;
            m4 = x*mm;
            Assert(Equal(a,m4,x*eps),label+" a *= x*b");
            CopyBackM(a0,a);
            a *= z*b;
            m4 = z*mm;
            Assert(Equal(a,m4,x*eps),label+" a *= z*b");
            CopyBackM(a0,a);
            a *= -x*b;
            m4 = -x*mm;
            Assert(Equal(a,m4,x*eps),label+" a *= -x*b");
            CopyBackM(a0,a);
            a *= -z*b;
            m4 = -z*mm;
            Assert(Equal(a,m4,x*eps),label+" a *= -z*b");
            CopyBackM(a0,a);
#ifdef ALIASOK
            a = -a*b;
            m4 = -mm;
            Assert(Equal(a,m4,eps),label+" a = -a*b");
            CopyBackM(a0,a);
            a += a*b;
            m4 = m0 + mm;
            Assert(Equal(a,m4,eps2),label+" a += a*b");
            CopyBackM(a0,a);
            a += x*a*b;
            m4 = m0 + x*mm;
            Assert(Equal(a,m4,x*eps2),label+" a += x*a*b");
            CopyBackM(a0,a);
            a += z*a*b;
            m4 = m0 + z*mm;
            Assert(Equal(a,m4,x*eps2),label+" a += z*a*b");
            CopyBackM(a0,a);
            a -= a*b;
            m4 = m0 - mm;
            Assert(Equal(a,m4,eps2),label+" a -= a*b");
            CopyBackM(a0,a);
            a -= x*a*b;
            m4 = m0 - x*mm;
            Assert(Equal(a,m4,x*eps2),label+" a -= x*a*b");
            CopyBackM(a0,a);
            a -= z*a*b;
            m4 = m0 - z*mm;
            Assert(Equal(a,m4,x*eps2),label+" a -= z*a*b");
            CopyBackM(a0,a);
#endif
        }
#endif // NOMULTEQ
#endif
    }

    if (showstartdone) std::cout<<"Done MM4a_Full"<<std::endl;
#endif
}

template <class M1, class M2> 
static void DoTestMM4RR(M1& a, const M2& b, std::string label)
{
    DoTestMM4a_Basic(a,b,label);
}

template <class M1, class M2> 
static void DoTestMM4RC(M1& a, const M2& b, std::string label)
{
    DoTestMM4a_Basic(a,b,label);
#if (XTEST & 2)
    DoTestMM4a_Basic(a,Conjugate(b),label+" ConjB");
#endif
}

template <class M1, class M2> 
static void DoTestMM4CR(M1& a, const M2& b, std::string label)
{
    typename M1::conjugate_type ac = a.conjugate();
    DoTestMM4a_Basic(a,b,label);
    DoTestMM4a_Basic(ac,b,label+" ConjA");
}

template <class M1, class M2> 
static void DoTestMM4CC(M1& a, const M2& b, std::string label)
{
    typename M1::conjugate_type ac = a.conjugate();
    DoTestMM4a_Full(a,b,label);
    DoTestMM4a_Basic(ac,b,label+" ConjA");

#if (XTEST & 2)
    DoTestMM4a_Basic(a,Conjugate(b),label+" ConjB");
    DoTestMM4a_Basic(ac,Conjugate(b),label+" ConjA ConjB");
#endif
}

template <class M1, class M2, class M3> 
static void DoTestMM5a_Basic(
    const M1& a, const M2& b, M3& c, std::string label)
{
    typedef typename M1::value_type Ta;
    typedef typename M2::value_type Tb;
    typedef typename M3::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start MM5a"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
        std::cout<<"c = "<<tmv::TMV_Text(c)<<std::endl;
    }

    const tmv::Matrix<Ta> m1 = a;
    const tmv::Matrix<Tb> m2 = b;
    const tmv::Matrix<T> m3 = c;
    const typename M3::copy_type c0 = c;
    tmv::Matrix<T> m4 = c;

    const int M = c.colsize();
    const int N = c.rowsize();
    const int K = a.rowsize();
    RT eps = EPS * (M + N + K);
    RT eps2 = eps;
    if (!std::numeric_limits<RT>::is_integer) {
        eps *= Norm(m1) * Norm(m2);
        eps2 *= Norm(m1*m2) + Norm(m1) * Norm(m2);
    }

#ifdef XXD
    if (XXDEBUG9) {
        std::cout<<"a = "<<TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"b = "<<TMV_Text(b)<<"  "<<b<<std::endl;
        std::cout<<"c = "<<TMV_Text(c)<<"  "<<c<<std::endl;
    }
#endif

    if (CanMultMM(a,b,c)) {
        tmv::Matrix<T> mm = m1*m2;
        c = a*b;
        m4 = mm;
#ifdef XXD
        if (XXDEBUG9) {
            std::cout<<"c = a*b = "<<c<<std::endl;
            std::cout<<"m1*m2 = "<<m4<<std::endl;
        }
#endif
        Assert(Equal(c,m4,eps),label+" c=a*b");
        CopyBackM(c0,c);
#ifndef BASIC_MULTMM_ONLY
        c += a*b;
        m4 = m3 + mm;
#ifdef XXD
        if (XXDEBUG9) {
            std::cout<<"c += a*b = "<<c<<std::endl;
            std::cout<<"m3 + m1*m2 = "<<m4<<std::endl;
        }
#endif
        Assert(Equal(c,m4,eps2),label+" c+=a*b");
        CopyBackM(c0,c);
        RT x(5);
        T z; SetZ(z);
        c = x*a*b;
        m4 = x*mm;
        Assert(Equal(c,m4,x*eps),label+" c=x*a*b");
        CopyBackM(c0,c);
        c = z*a*b;
        m4 = z*mm;
#ifdef XXD
        if (XXDEBUG9) {
            std::cout<<"c = z*a*b = "<<c<<std::endl;
            std::cout<<"m4 = z*m1*m2 = "<<m4<<std::endl;
            std::cout<<"c-m4 = "<<(c-m4)<<std::endl;
            std::cout<<"Norm(c-m4) = "<<Norm(c-m4)<<std::endl;
        }
#endif
        Assert(Equal(c,m4,x*eps),label+" c=z*a*b");
        CopyBackM(c0,c);
#endif
    }
#ifndef NOELEMMULT
    if (CanElemMultMM(a,b)) {
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m4(i,j) = m1(i,j) * m2(i,j);
        c = ElemProd(a,b);
#ifdef XXD
        if (XXDEBUG8) {
            std::cout<<"c = ElemProd(a,b) = "<<c<<std::endl;
            std::cout<<"m4 = "<<m4<<std::endl;
            std::cout<<"diff = "<<(c-m4)<<std::endl;
            std::cout<<"Norm(diff) = "<<Norm(c-m4)<<std::endl;
            std::cout<<"cf. eps = "<<eps<<std::endl;
        }
#endif
        Assert(Equal(c,m4,eps),label+" c = ElemProd(a,b)");
        CopyBackM(c0,c);
        c = ElemProd(b,a);
        Assert(Equal(c,m4,eps),label+" c = ElemProd(b,a)");
        CopyBackM(c0,c);
    }
    if (CanAddElemMultMM(a,b,c)) {
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m4(i,j) = m3(i,j) + m1(i,j) * m2(i,j);
        c += ElemProd(a,b);
        Assert(Equal(c,m4,eps),label+" c += ElemProd(a,b)");
        CopyBackM(c0,c);
        c += ElemProd(b,a);
        Assert(Equal(c,m4,eps),label+" c += ElemProd(b,a)");
        CopyBackM(c0,c);

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m4(i,j) = m3(i,j) - m1(i,j) * m2(i,j);
        c -= ElemProd(a,b);
        Assert(Equal(c,m4,eps),label+" c -= ElemProd(a,b)");
        CopyBackM(c0,c);
        c -= ElemProd(b,a);
        Assert(Equal(c,m4,eps),label+" c -= ElemProd(b,a)");
        CopyBackM(c0,c);

        RT x = 5;
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m4(i,j) = x * m1(i,j) * m2(i,j);
        c = x*ElemProd(a,b);
        Assert(Equal(c,m4,eps),label+" c = x*ElemProd(a,b)");
        CopyBackM(c0,c);
        c = x*ElemProd(b,a);
        Assert(Equal(c,m4,eps),label+" c = x*ElemProd(b,a)");
        CopyBackM(c0,c);
        c = ElemProd(x*a,b);
        Assert(Equal(c,m4,eps),label+" c = ElemProd(x*a,b)");
        CopyBackM(c0,c);
        c = ElemProd(b,x*a);
        Assert(Equal(c,m4,eps),label+" c = ElemProd(b,x*a)");
        CopyBackM(c0,c);
        c = ElemProd(a,x*b);
        Assert(Equal(c,m4,eps),label+" c = ElemProd(a,x*b)");
        CopyBackM(c0,c);
        c = ElemProd(x*b,a);
        Assert(Equal(c,m4,eps),label+" c = ElemProd(x*b,a)");
        CopyBackM(c0,c);

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m4(i,j) = m3(i,j) + x * m1(i,j) * m2(i,j);
        c += x*ElemProd(a,b);
        Assert(Equal(c,m4,eps),label+" c += x*ElemProd(a,b)");
        CopyBackM(c0,c);
        c += x*ElemProd(b,a);
        Assert(Equal(c,m4,eps),label+" c += x*ElemProd(b,a)");
        CopyBackM(c0,c);
        c += ElemProd(x*a,b);
        Assert(Equal(c,m4,eps),label+" c += ElemProd(x*a,b)");
        CopyBackM(c0,c);
        c += ElemProd(b,x*a);
        Assert(Equal(c,m4,eps),label+" c += ElemProd(b,x*a)");
        CopyBackM(c0,c);
        c += ElemProd(a,x*b);
        Assert(Equal(c,m4,eps),label+" c += ElemProd(a,x*b)");
        CopyBackM(c0,c);
        c += ElemProd(x*b,a);
        Assert(Equal(c,m4,eps),label+" c += ElemProd(x*b,a)");
        CopyBackM(c0,c);

        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            m4(i,j) = m3(i,j) - x * m1(i,j) * m2(i,j);
        c -= x*ElemProd(a,b);
        Assert(Equal(c,m4,eps),label+" c -= x*ElemProd(a,b)");
        CopyBackM(c0,c);
        c -= x*ElemProd(b,a);
        Assert(Equal(c,m4,eps),label+" c -= x*ElemProd(b,a)");
        CopyBackM(c0,c);
        c -= ElemProd(x*a,b);
        Assert(Equal(c,m4,eps),label+" c -= ElemProd(x*a,b)");
        CopyBackM(c0,c);
        c -= ElemProd(b,x*a);
        Assert(Equal(c,m4,eps),label+" c -= ElemProd(b,x*a)");
        CopyBackM(c0,c);
        c -= ElemProd(a,x*b);
        Assert(Equal(c,m4,eps),label+" c -= ElemProd(a,x*b)");
        CopyBackM(c0,c);
        c -= ElemProd(x*b,a);
        Assert(Equal(c,m4,eps),label+" c -= ElemProd(x*b,a)");
        CopyBackM(c0,c);
    }
#endif

    if (showstartdone) std::cout<<"Done MM5a"<<std::endl;
}

template <class M1, class M2, class M3> 
static void DoTestMM5a_Full(const M1& a, const M2& b, M3& c, std::string label)
{
    DoTestMM5a_Basic(a,b,c,label);

#ifndef BASIC_MULTMM_ONLY
#if (XTEST & 16)
    typedef typename M1::value_type Ta;
    typedef typename M2::value_type Tb;
    typedef typename M3::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Continuing on to MM5a Full"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<std::endl;
        std::cout<<"c = "<<tmv::TMV_Text(c)<<std::endl;
    }
    tmv::Matrix<Ta> m = a;
    tmv::Matrix<Tb> m2 = b;
    tmv::Matrix<T> m3 = c;

    RT eps = EPS * (a.colsize() + a.rowsize() + b.rowsize());
    RT eps2 = eps;
    if (!std::numeric_limits<RT>::is_integer) {
        eps *= Norm(m) * Norm(m2);
        eps2 *= Norm(m*m2) + Norm(m) * Norm(m2);
    }

    if (CanMultMM(a,b,c)) {
        typename M3::copy_type c0 = c;
        tmv::Matrix<T> mm = m*m2;
        RT x(5);
        T z; SetZ(z);
        c += x*a*b;
        m3 = c0 + x*mm;
        Assert(Equal(c,m3,x*eps2),label+" c+=x*a*b");
        CopyBackM(c0,c);
        c += z*a*b;
        m3 = c0 + z*mm;
        Assert(Equal(c,m3,x*eps2),label+" c+=z*a*b");
        CopyBackM(c0,c);
        c = -a*b;
        m3 = -mm;
        Assert(Equal(c,m3,eps),label+" c=-a*b");
        CopyBackM(c0,c);
        c += -a*b;
        m3 = c0 - mm;
        Assert(Equal(c,m3,eps2),label+" c+=-a*b");
        CopyBackM(c0,c);
        c -= a*b;
        m3 = c0 - mm;
        Assert(Equal(c,m3,eps2),label+" c-=a*b");
        CopyBackM(c0,c);
        c -= -a*b;
        m3 = c0 + mm;
        Assert(Equal(c,m3,eps2),label+" c-=-a*b");
        CopyBackM(c0,c);
        c += -x*a*b;
        m3 = c0 - x*mm;
        Assert(Equal(c,m3,x*eps2),label+" c+=-x*a*b");
        CopyBackM(c0,c);
        c -= x*a*b;
        m3 = c0 - x*mm;
        Assert(Equal(c,m3,x*eps2),label+" c-=x*a*b");
        CopyBackM(c0,c);
        c -= -x*a*b;
        m3 = c0 + x*mm;
        Assert(Equal(c,m3,x*eps2),label+" c-=-x*a*b");
        CopyBackM(c0,c);
        c += -z*a*b;
        m3 = c0 - z*mm;
        Assert(Equal(c,m3,x*eps2),label+" c+=-z*a*b");
        CopyBackM(c0,c);
        c -= z*a*b;
        m3 = c0 - z*mm;
        Assert(Equal(c,m3,x*eps2),label+" c-=z*a*b");
        CopyBackM(c0,c);
        c -= -z*a*b;
        m3 = c0 + z*mm;
        Assert(Equal(c,m3,x*eps2),label+" c-=-z*a*b");
        CopyBackM(c0,c);
    }

    if (showstartdone) std::cout<<"Done MM5a_Full"<<std::endl;
#endif
#endif
}

template <class M1, class M2, class M3> 
static void DoTestMM5R(const M1& a, const M2& b, M3& c, std::string label)
{
    DoTestMM5a_Basic(a,b,c,label);
}

template <class M1, class M2, class M3> 
static void DoTestMM5C(const M1& a, const M2& b, M3& c, std::string label)
{
    DoTestMM5a_Full(a,b,c,label);
    DoTestMM5a_Basic(Conjugate(a),b,c,label+" ConjA");
#if (XTEST & 2)
    typename M3::conjugate_type cc = c.conjugate();
    DoTestMM5a_Basic(a,Conjugate(b),c,label+" ConjB");
    DoTestMM5a_Basic(Conjugate(a),Conjugate(b),c,label+" ConjAB");
    DoTestMM5a_Basic(a,b,cc,label+" ConjC");
    DoTestMM5a_Basic(a,Conjugate(b),cc,label+" ConjBC");
    DoTestMM5a_Basic(Conjugate(a),b,cc,label+" ConjAC");
    DoTestMM5a_Basic(Conjugate(a),Conjugate(b),cc,label+" ConjABC");
#endif
}

template <class M, class V1, class V2> 
static void DoTestOProda_Basic(
    M& a, const V1& v1, const V2& v2, std::string label)
{
    typedef typename M::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Start OProd"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"v1 = "<<tmv::TMV_Text(v1)<<std::endl;
        std::cout<<"v2 = "<<tmv::TMV_Text(v2)<<std::endl;
    }
#ifdef XXD
    if (XXDEBUG9) {
        std::cout<<"a = "<<a<<std::endl;
        std::cout<<"v1 = "<<v1<<std::endl;
        std::cout<<"v2 = "<<v2<<std::endl;
    }
#endif

    tmv::Matrix<T> vv = tmv::Vector<T>(v1)^tmv::Vector<T>(v2);

    RT eps = EPS * (a.colsize() + a.rowsize());
    RT eps2 = eps;
    if (!std::numeric_limits<RT>::is_integer) {
        eps *= Norm(v1) * Norm(v2);
        eps2 *= Norm(vv) + Norm(v1) * Norm(v2);
    }

    a = v1^v2;
#ifdef XXD
    if (XXDEBUG9) {
        std::cout<<"v1 ^ v2 = "<<vv<<std::endl;
        std::cout<<"a = v1 ^ v2 = "<<a<<std::endl;
        std::cout<<"Norm(a-v^v) = "<<Norm(a-vv)<<"  cf "<<eps<<std::endl;
    }
#endif
    Assert(Equal(a,vv,eps),label+" a = v1^v2");
    typename M::copy_type a0 = a;
    a += v1^v2;
    Assert(Equal(a,(a0+vv),eps2),label+" a += v1^v2");
    CopyBackM(a0,a);
    a -= v1^v2;
    Assert(Equal(a,(a0-vv),eps2),label+" a -= v1^v2");
    CopyBackM(a0,a);
    RT x(5);
    T z; SetZ(z);
    a = x * (v1^v2);
    Assert(Equal(a,x*vv,x*eps),label+" a = x * (v1^v2)");
    CopyBackM(a0,a);
#ifdef SYMOPROD
    if (a.issym()) {
#endif
        a = z * (v1^v2);
#ifdef XXD
        if (XXDEBUG9) {
            std::cout<<"z * v1 ^ v2 = "<<z*vv<<std::endl;
            std::cout<<"a = z * v1 ^ v2 = "<<a<<std::endl;
            std::cout<<"Norm(a-v^v) = "<<Norm(a-z*vv)
                <<"  cf "<<x*eps<<std::endl;
        }
#endif
        Assert(Equal(a,z*vv,x*eps),label+" a = z * (v1^v2)");
        CopyBackM(a0,a);
#ifdef SYMOPROD
    }
#endif

    if (showstartdone) std::cout<<"Done OProd"<<std::endl;
}

template <class M, class V1, class V2> 
static void DoTestOProda_Full(
    M& a, const V1& v1, const V2& v2, std::string label)
{
    DoTestOProda_Basic(a,v1,v2,label);

#if (XTEST & 16)
    typedef typename M::value_type T;
    typedef typename tmv::Traits<T>::real_type RT;
    if (showstartdone) {
        std::cout<<"Continuing on to OProda Full"<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<std::endl;
        std::cout<<"v1 = "<<tmv::TMV_Text(v1)<<std::endl;
        std::cout<<"v2 = "<<tmv::TMV_Text(v2)<<std::endl;
    }
    tmv::Matrix<T> vv = tmv::Vector<T>(v1)^tmv::Vector<T>(v2);
    typename M::copy_type a0 = vv;

    RT eps = EPS * (a.colsize() + a.rowsize());
    RT eps2 = eps;
    if (!std::numeric_limits<RT>::is_integer) {
        eps *= Norm(v1) * Norm(v2);
        eps2 *= Norm(vv) + Norm(v1) * Norm(v2);
    }

    RT x(5);
    T z; SetZ(z);
    a = (x * v1)^v2;
    Assert(Equal(a,x*vv,x*eps),label+" a = (x*v1) ^ v2)");
    CopyBackM(a0,a);
    a = v1 ^ (x * v2);
    Assert(Equal(a,x*vv,x*eps),label+" a = v1 ^ (x*v2)");
    CopyBackM(a0,a);
    a += x * (v1^v2);
    Assert(Equal(a,(a0+x*vv),x*eps2),label+" a += x * (v1^v2)");
    CopyBackM(a0,a);
    a += (x * v1)^v2;
    Assert(Equal(a,(a0+x*vv),x*eps2),label+" a += (x*v1) ^ v2)");
    CopyBackM(a0,a);
    a += v1 ^ (x * v2);
    Assert(Equal(a,(a0+x*vv),x*eps2),label+" a += v1 ^ (x*v2)");
    CopyBackM(a0,a);
    a -= x * (v1^v2);
    Assert(Equal(a,(a0-x*vv),x*eps2),label+" a -= x * (v1^v2)");
    CopyBackM(a0,a);
    a -= (x * v1)^v2;
    Assert(Equal(a,(a0-x*vv),x*eps2),label+" a -= (x*v1) ^ v2)");
    CopyBackM(a0,a);
    a -= v1 ^ (x * v2);
    Assert(Equal(a,(a0-x*vv),x*eps2),label+" a -= v1 ^ (x*v2)");
    CopyBackM(a0,a);

#ifdef SYMOPROD
    if (a.issym()) {
#endif
        a = (z * v1)^v2;
        Assert(Equal(a,z*vv,x*eps),label+" a = (z*v1) ^ v2)");
        CopyBackM(a0,a);
        a = v1 ^ (z * v2);
        Assert(Equal(a,z*vv,x*eps),label+" a = v1 ^ (z*v2)");
        CopyBackM(a0,a);
        a += z * (v1^v2);
        Assert(Equal(a,(a0+z*vv),x*eps2),label+" a += z * (v1^v2)");
        CopyBackM(a0,a);
        a += (z * v1)^v2;
        Assert(Equal(a,(a0+z*vv),x*eps2),label+" a += (z*v1) ^ v2)");
        CopyBackM(a0,a);
        a += v1 ^ (z * v2);
        Assert(Equal(a,(a0+z*vv),x*eps2),label+" a += v1 ^ (z*v2)");
        CopyBackM(a0,a);
        a -= z * (v1^v2);
        Assert(Equal(a,(a0-z*vv),x*eps2),label+" a -= z * (v1^v2)");
        CopyBackM(a0,a);
        a -= (z * v1)^v2;
        Assert(Equal(a,(a0-z*vv),x*eps2),label+" a -= (z*v1) ^ v2)");
        CopyBackM(a0,a);
        a -= v1 ^ (z * v2);
        Assert(Equal(a,(a0-z*vv),x*eps2),label+" a -= v1 ^ (z*v2)");
        CopyBackM(a0,a);
        a = (x * v1)^(x * v2);
        Assert(Equal(a,x*x*vv,x*x*eps),label+" a = (x*v1) ^ (x*v2))");
        CopyBackM(a0,a);
        a += (x * v1)^(x * v2);
        Assert(Equal(a,(a0+x*x*vv),x*x*eps2),label+" a += (x*v1) ^ (x*v2)");
        CopyBackM(a0,a);
        a -= (x * v1)^(x * v2);
        Assert(Equal(a,(a0-x*x*vv),x*x*eps2),label+" a -= (x*v1) ^ (x*v2)");
        CopyBackM(a0,a);
        a = (x * v1)^(z * v2);
        Assert(Equal(a,x*z*vv,x*x*eps),label+" a = (x*v1) ^ (z*v2))");
        CopyBackM(a0,a);
        a += (x * v1)^(z * v2);
        Assert(Equal(a,(a0+x*z*vv),x*x*eps2),label+" a += (x*v1) ^ (z*v2)");
        CopyBackM(a0,a);
        a -= (x * v1)^(z * v2);
        Assert(Equal(a,(a0-x*z*vv),x*x*eps2),label+" a -= (x*v1) ^ (z*v2)");
        CopyBackM(a0,a);
        a = (z * v1)^(x * v2);
        Assert(Equal(a,z*x*vv,x*x*eps),label+" a = (z*v1) ^ (x*v2))");
        CopyBackM(a0,a);
        a += (z * v1)^(x * v2);
        Assert(Equal(a,(a0+z*x*vv),x*x*eps2),label+" a += (z*v1) ^ (x*v2)");
        CopyBackM(a0,a);
        a -= (z * v1)^(x * v2);
        Assert(Equal(a,(a0-z*x*vv),x*x*eps2),label+" a -= (z*v1) ^ (x*v2)");
        CopyBackM(a0,a);
        a = (z * v1)^(z * v2);
        Assert(Equal(a,z*z*vv,x*x*eps),label+" a = (z*v1) ^ (z*v2))");
        CopyBackM(a0,a);
        a += (z * v1)^(z * v2);
        Assert(Equal(a,(a0+z*z*vv),x*x*eps2),label+" a += (z*v1) ^ (z*v2)");
        CopyBackM(a0,a);
        a -= (z * v1)^(z * v2);
        Assert(Equal(a,(a0-z*z*vv),x*x*eps2),label+" a -= (z*v1) ^ (z*v2)");
        CopyBackM(a0,a);
#ifdef SYMOPROD
    }
#endif

    if (showstartdone) std::cout<<"Done OProda_Full"<<std::endl;
#endif
} 

template <class M, class V1, class V2> 
static void DoTestOProdR(M& a, const V1& v1, const V2& v2, std::string label)
{
    DoTestOProda_Basic(a,v1,v2,label);
    DoTestOProda_Basic(a,v1.reverse(),v2.reverse(),label+" RevBC");
#ifndef SYMOPROD
    DoTestOProda_Basic(a,v1.reverse(),v2,label+" RevB");
    DoTestOProda_Basic(a,v1,v2.reverse(),label+" RevC");
#endif
}

template <class M, class V1, class V2> 
static void DoTestOProdC(M& a, const V1& v1, const V2& v2, std::string label)
{
    typename M::conjugate_type ac = a.conjugate();

    DoTestOProda_Full(a,v1,v2,label);
    DoTestOProda_Basic(a,v1.reverse(),v2.reverse(),label+" RevBC");
    DoTestOProda_Basic(ac,v1,v2,label+" ConjA");
    DoTestOProda_Basic(ac,v1.reverse(),v2.reverse(),
                          label+" ConjA RevBC");
#ifndef SYMOPROD
    DoTestOProda_Basic(a,v1.reverse(),v2,label+" RevB");
    DoTestOProda_Basic(a,v1,v2.reverse(),label+" RevC");
    DoTestOProda_Basic(ac,v1.reverse(),v2,label+" ConjA RevB");
    DoTestOProda_Basic(ac,v1,v2.reverse(),label+" ConjA RevC");
#endif
}

template <class M, class CM> 
static void TestMatrixArith1(M& a, CM& ca, std::string label)
{
    typedef typename M::value_type T;
    typedef std::complex<T> CT;
    if (showstartdone) {
        std::cout<<"Start TestMatrixArith1 "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
    }

    CT z(9,-2);
    T x = 12;

#if 1
#ifndef NO_REAL_ARITH
    DoTestMR(a,label+" R");
#endif
#ifndef NO_COMPLEX_ARITH
    DoTestMC(ca,label+" C");
#endif
#endif

#if 1
#ifndef NO_REAL_ARITH
    DoTestMX1R(a,x,label+" R,R");
#endif
#ifndef NO_COMPLEX_ARITH
    DoTestMX1C(ca,z,label+" C,C");
#if (XTEST & 4)
#ifndef NO_REAL_ARITH
    DoTestMX1R(a,z,label+" R,C");
    DoTestMX1C(ca,x,label+" C,R");
#endif
#endif
#endif
#endif

#ifndef NOASSIGN

#if 1
#ifndef NO_REAL_ARITH
    DoTestMX2R(a,x,label+" R,R");
#endif
#ifndef NO_COMPLEX_ARITH
    DoTestMX2C(ca,z,label+" C,C");
#if (XTEST & 4)
#ifndef NO_REAL_ARITH
    DoTestMX2C(ca,x,label+" C,R");
#endif
#endif
#endif
#endif

#if 1
#ifdef ALIASOK
    DoTestMM2RR(a,a,label+" self_arith");
    DoTestMM4RR(a,a,label+" self_arith");
#ifndef NO_COMPLEX_ARITH
    DoTestMM2CC(ca,ca,label+" self_arith");
    DoTestMM4CC(ca,ca,label+" self_arith");
#endif
#endif
#endif

#endif

    if (showstartdone) std::cout<<"Done Test1"<<std::endl;
}

template <class M, class CM, class V1, class CV1, class V2, class CV2> 
static void TestMatrixArith2a(
    const M& a, const CM& ca, V1& b, CV1& cb, 
    V2& c, CV2& cc, std::string label)
{
    if (showstartdone) {
        std::cout<<"Start TestMatrixArith2a "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
        std::cout<<"cb = "<<tmv::TMV_Text(cb)<<"  "<<cb<<std::endl;
        std::cout<<"c = "<<tmv::TMV_Text(c)<<"  "<<c<<std::endl;
        std::cout<<"cc = "<<tmv::TMV_Text(cc)<<"  "<<cc<<std::endl;
    }
#ifndef NO_REAL_ARITH
    DoTestMV1R(a,b,label+" R,R");
#endif
#ifndef NO_COMPLEX_ARITH
    DoTestMV1C(ca,cb,label+" C,C");
#if (XTEST & 4)
#ifndef NO_REAL_ARITH
    DoTestMV1R(a,cb,label+" R,C");
    DoTestMV1C(ca,b,label+" C,R");
#endif
#endif
#endif

#ifndef NONSQUARE
#ifndef NO_REAL_ARITH
    DoTestMV2R(a,b,label+" R,R");
#endif
#ifndef NO_COMPLEX_ARITH
    DoTestMV2C(ca,cb,label+" C,C");
#if (XTEST & 4)
#ifndef NO_REAL_ARITH
    DoTestMV2R(a,cb,label+" R,C");
#endif
#endif
#endif
#endif

    if (showstartdone) std::cout<<"Done Test2a"<<std::endl;
}

template <class M, class CM, class V1, class CV1, class V2, class CV2> 
static void TestMatrixArith2b(
    const M& a, const CM& ca, V1& b, CV1& cb, 
    V2& c, CV2& cc, std::string label)
{
    if (showstartdone) {
        std::cout<<"Start TestMatrixArith2b "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
        std::cout<<"cb = "<<tmv::TMV_Text(cb)<<"  "<<cb<<std::endl;
        std::cout<<"c = "<<tmv::TMV_Text(c)<<"  "<<c<<std::endl;
        std::cout<<"cc = "<<tmv::TMV_Text(cc)<<"  "<<cc<<std::endl;
    }

#ifndef NO_REAL_ARITH
    DoTestVM1R(a,c,label+" R,R");
#endif
#ifndef NO_COMPLEX_ARITH
    DoTestVM1C(ca,cc,label+" C,C");
#if (XTEST & 4)
#ifndef NO_REAL_ARITH
    DoTestVM1R(a,cc,label+" R,C");
    DoTestVM1C(ca,c,label+" C,R");
#endif
#endif
#endif

#ifndef NONSQUARE
#ifndef NO_REAL_ARITH
    DoTestVM2R(a,c,label+" R,R");
#endif
#ifndef NO_COMPLEX_ARITH
    DoTestVM2C(ca,cc,label+" C,C");
#if (XTEST & 4)
#ifndef NO_REAL_ARITH
    DoTestVM2R(a,cc,label+" R,C");
#endif
#endif
#endif
#endif

    if (showstartdone) std::cout<<"Done Test2b"<<std::endl;
}

template <class M, class CM> 
static void TestMatrixArith2(M& a, CM& ca, std::string label)
{
    typedef typename M::value_type T;
    typedef std::complex<T> CT;
    tmv::Vector<T> v(a.rowsize());
    for(int i=0;i<int(a.rowsize());i++) v(i) = T(i+3);
    tmv::Vector<CT> cv = CT(4,5) * v;

    tmv::Vector<T> v5(5*a.rowsize(),T(0));
    tmv::Vector<CT> cv5(5*a.rowsize(),CT(0));
    tmv::VectorView<T> vs = v5.cSubVector(0,5*a.rowsize(),5);
    tmv::VectorView<CT> cvs = cv5.cSubVector(0,5*a.rowsize(),5);
    vs = v;
    cvs = cv;

    tmv::Vector<T> w(a.colsize());
    for(int i=0;i<int(a.colsize());i++) w(i) = T(2*i-6);
    tmv::Vector<CT> cw = CT(-1,2) * w;

    tmv::Vector<T> w5(5*a.colsize(),T(0));
    tmv::Vector<CT> cw5(5*a.colsize(),CT(0));
    tmv::VectorView<T> ws = w5.cSubVector(0,5*a.colsize(),5);
    tmv::VectorView<CT> cws = cw5.cSubVector(0,5*a.colsize(),5);
    ws = w;
    cws = cw;

    TestMatrixArith2a(a,ca,v,cv,w,cw,label);
    TestMatrixArith2a(a,ca,vs,cvs,w,cw,label);
    TestMatrixArith2b(a,ca,v,cv,w,cw,label);
    TestMatrixArith2b(a,ca,v,cv,ws,cws,label);
}

template <class M, class CM, class V1, class CV1, class V2, class CV2> 
static void TestMatrixArith3a(
    const M& a, const CM& ca, V1& b, CV1& cb,
    V2& c, CV2& cc, std::string label)
{
    if (showstartdone) {
        std::cout<<"Start TestMatrixArith3a "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
        std::cout<<"cb = "<<tmv::TMV_Text(cb)<<"  "<<cb<<std::endl;
        std::cout<<"c = "<<tmv::TMV_Text(c)<<"  "<<c<<std::endl;
        std::cout<<"cc = "<<tmv::TMV_Text(cc)<<"  "<<cc<<std::endl;
    }

#ifndef NO_REAL_ARITH
    DoTestMV3R(a,b,c,label+" R,R,R");
#endif
#ifndef NO_COMPLEX_ARITH
    DoTestMV3C(ca,cb,cc,label+" C,C,C");
#if (XTEST & 4)
#ifndef NO_REAL_ARITH
    DoTestMV3R(a,b,cc,label+" C,R,R");
    DoTestMV3R(a,cb,cc,label+" C,R,C");
    DoTestMV3C(ca,b,cc,label+" C,C,R");
#endif
#endif
#endif

    if (showstartdone) std::cout<<"Done Test3a"<<std::endl;
}

template <class M, class CM, class V1, class CV1, class V2, class CV2> 
static void TestMatrixArith3b(
    const M& a, const CM& ca, V1& b, CV1& cb,
    V2& c, CV2& cc, std::string label)
{
    if (showstartdone) {
        std::cout<<"Start TestMatrixArith3b "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
        std::cout<<"cb = "<<tmv::TMV_Text(cb)<<"  "<<cb<<std::endl;
        std::cout<<"c = "<<tmv::TMV_Text(c)<<"  "<<c<<std::endl;
        std::cout<<"cc = "<<tmv::TMV_Text(cc)<<"  "<<cc<<std::endl;
    }

#ifndef NO_REAL_ARITH
    DoTestVM3R(a,c,b,label+" R,R,R");
#endif
#ifndef NO_COMPLEX_ARITH
    DoTestVM3C(ca,cc,cb,label+" C,C,C");
#if (XTEST & 4)
#ifndef NO_REAL_ARITH
    DoTestVM3R(a,c,cb,label+" C,R,R");
    DoTestVM3R(a,cc,cb,label+" C,R,C");
    DoTestVM3C(ca,c,cb,label+" C,C,R");
#endif
#endif
#endif

    if (showstartdone) std::cout<<"Done Test3b"<<std::endl;
}

template <class M, class CM> 
static void TestMatrixArith3(M& a, CM& ca, std::string label)
{
    typedef typename M::value_type T;
    typedef std::complex<T> CT;
    tmv::Vector<T> vx(a.rowsize());
    for(int i=0;i<int(a.rowsize());i++) vx(i) = T(i+3);
    tmv::Vector<CT> cvx = CT(4,5) * vx;
    tmv::VectorView<T> v = vx.view();
    tmv::VectorView<CT> cv = cvx.view();

    tmv::Vector<T> v5(5*a.rowsize(),T(0));
    tmv::Vector<CT> cv5(5*a.rowsize(),CT(0));
    tmv::VectorView<T> vs = v5.cSubVector(0,5*a.rowsize(),5);
    tmv::VectorView<CT> cvs = cv5.cSubVector(0,5*a.rowsize(),5);
    vs = v;
    cvs = cv;

    tmv::Vector<T> wx(a.colsize());
    for(int i=0;i<int(a.colsize());i++) wx(i) = T(2*i-6);
    tmv::Vector<CT> cwx = CT(-1,2) * wx;
    tmv::VectorView<T> w = wx.view();
    tmv::VectorView<CT> cw = cwx.view();

    tmv::Vector<T> w5(5*a.colsize(),T(0));
    tmv::Vector<CT> cw5(5*a.colsize(),CT(0));
    tmv::VectorView<T> ws = w5.cSubVector(0,5*a.colsize(),5);
    tmv::VectorView<CT> cws = cw5.cSubVector(0,5*a.colsize(),5);
    ws = w;
    cws = cw;

    TestMatrixArith3a(a,ca,v,cv,w,cw,label);
    TestMatrixArith3a(a,ca,v,cv,ws,cws,label);
    TestMatrixArith3a(a,ca,vs,cvs,w,cw,label);
    TestMatrixArith3a(a,ca,vs,cvs,ws,cws,label); 

    TestMatrixArith3b(a,ca,v,cv,w,cw,label);
    TestMatrixArith3b(a,ca,v,cv,ws,cws,label);
    TestMatrixArith3b(a,ca,vs,cvs,w,cw,label);
    TestMatrixArith3b(a,ca,vs,cvs,ws,cws,label); 
}

template <class M1, class CM1, class M2, class CM2> 
static void TestMatrixArith4(
    M1& a, CM1& ca, const M2& b, const CM2& cb, std::string label)
{
    if (showstartdone) {
        std::cout<<"Start TestMatrixArith4 "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
        std::cout<<"cb = "<<tmv::TMV_Text(cb)<<"  "<<cb<<std::endl;
    }

#ifndef NO_REAL_ARITH
    DoTestMM1RR(a,b,label+" R,R");
#ifndef NOASSIGN
    DoTestMM2RR(a,b,label+" R,R");
#endif
#endif
#ifndef NO_COMPLEX_ARITH
    DoTestMM1CC(ca,cb,label+" C,C");
#ifndef NOASSIGN
    DoTestMM2CC(ca,cb,label+" C,C");
#endif
#if (XTEST & 4)
#ifndef NO_REAL_ARITH
    DoTestMM1RC(a,cb,label+" R,C");
    DoTestMM1CR(ca,b,label+" C,R");
#ifndef NOASSIGN
    DoTestMM2CR(ca,b,label+" C,R");
#endif
#endif
#endif
#endif

    if (showstartdone) std::cout<<"Done Test4"<<std::endl;
}

template <class M1, class CM1, class M2, class CM2>
static void TestMatrixArith5(
    M1& a, CM1& ca, const M2& b, const CM2& cb, std::string label)
{
    if (showstartdone) {
        std::cout<<"Start TestMatrixArith5 "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
        std::cout<<"cb = "<<tmv::TMV_Text(cb)<<"  "<<cb<<std::endl;
    }

#ifndef NO_REAL_ARITH
    DoTestMM3RR(a,b,label+" R,R");
#ifndef NOASSIGN
    DoTestMM4RR(a,b,label+" R,R");
#endif
#endif
#ifndef NO_COMPLEX_ARITH
    DoTestMM3CC(ca,cb,label+" C,C");
#ifndef NOASSIGN
    DoTestMM4CC(ca,cb,label+" C,C");
#endif
#if (XTEST & 4)
#ifndef NO_REAL_ARITH
    DoTestMM3RC(a,cb,label+" R,C");
    DoTestMM3CR(ca,b,label+" C,R");
#ifndef NOASSIGN
    DoTestMM4CR(ca,b,label+" C,R");
#endif
#endif
#endif
#endif

    if (showstartdone) std::cout<<"Done Test5"<<std::endl;
}

template <class M1, class CM1, class M2, class CM2, class M3, class CM3> 
static void TestMatrixArith6(
    const M1& a, const CM1& ca,
    const M2& b, const CM2& cb, M3& c, CM3& cc, std::string label)
{
    if (showstartdone) {
        std::cout<<"Start TestMatrixArith6 "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
        std::cout<<"b = "<<tmv::TMV_Text(b)<<"  "<<b<<std::endl;
        std::cout<<"cb = "<<tmv::TMV_Text(cb)<<"  "<<cb<<std::endl;
        std::cout<<"c = "<<tmv::TMV_Text(c)<<"  "<<c<<std::endl;
        std::cout<<"cc = "<<tmv::TMV_Text(cc)<<"  "<<cc<<std::endl;
    }

#ifndef NO_REAL_ARITH
    DoTestMM5R(a,b,c,label+" R,R,R");
#endif
#ifndef NO_COMPLEX_ARITH
    DoTestMM5C(ca,cb,cc,label+" C,C,C");
#if (XTEST & 4)
#ifndef NO_REAL_ARITH
    DoTestMM5R(a,b,cc,label+" C,R,R");
    DoTestMM5R(a,cb,cc,label+" C,R,C");
    DoTestMM5R(ca,b,cc,label+" C,C,R");
#endif
#endif
#endif

    if (showstartdone) std::cout<<"Done Test6"<<std::endl;
}

template <class M1, class CM1, class M2, class CM2> 
static void TestMatrixArith6x(
    const M1& a, const CM1& ca,
    const M2& b, const CM2& cb, std::string label)
{
    typedef typename M1::value_type T;
    typedef std::complex<T> CT;
    if (CanMultMM(a,b)) {
        tmv::Matrix<T,tmv::ColMajor> c1(a*b);
        tmv::Matrix<CT,tmv::ColMajor> cc1(ca*cb);
        TestMatrixArith6(a,ca,b,cb,c1,cc1,label);
#if (XTEST & 2)
        tmv::Matrix<T,tmv::RowMajor> c2(c1);
        tmv::Matrix<CT,tmv::RowMajor> cc2(cc1);
        TestMatrixArith6(a,ca,b,cb,c2,cc2,label);
#endif

#if (XTEST & 1)
        tmv::Matrix<T> c3x(4*c1.colsize(),5*c1.rowsize(),T(0));
        tmv::Matrix<CT> cc3x(4*c1.colsize(),5*c1.rowsize(),CT(0));
        tmv::MatrixView<T> c3 = c3x.cSubMatrix(0,c3x.colsize(),0,c3x.rowsize(),4,5);
        tmv::MatrixView<CT> cc3 = cc3x.cSubMatrix(0,c3x.colsize(),0,c3x.rowsize(),4,5);
        c3 = c1;
        cc3 = cc1;
        TestMatrixArith6(a,ca,b,cb,c3,cc3,label);
#endif
    }
}

template <class M, class CM, class V1, class CV1, class V2, class CV2> 
static void TestMatrixArith7(
    M& a, CM& ca, const V1& v1, const CV1& cv1,
    const V2& v2, const CV2& cv2, std::string label)
{
    if (showstartdone) {
        std::cout<<"Start TestMatrixArith6 "<<label<<std::endl;
        std::cout<<"a = "<<tmv::TMV_Text(a)<<"  "<<a<<std::endl;
        std::cout<<"ca = "<<tmv::TMV_Text(ca)<<"  "<<ca<<std::endl;
        std::cout<<"v1 = "<<tmv::TMV_Text(v1)<<"  "<<v1<<std::endl;
        std::cout<<"cv1 = "<<tmv::TMV_Text(cv1)<<"  "<<cv1<<std::endl;
        std::cout<<"v2 = "<<tmv::TMV_Text(v2)<<"  "<<v2<<std::endl;
        std::cout<<"cv2 = "<<tmv::TMV_Text(cv2)<<"  "<<cv2<<std::endl;
    }

#ifndef NO_REAL_ARITH
    DoTestOProdR(a,v1,v2,label+" R,R,R");
#endif
#ifndef NO_COMPLEX_ARITH
    DoTestOProdC(ca,cv1,cv2,label+" C,C,C");
#if (XTEST & 4)
#ifndef NO_REAL_ARITH
    DoTestOProdC(ca,v1,v2,label+" C,R,R");
#ifndef SYMOPROD
    DoTestOProdC(ca,cv1,v2,label+" C,C,R");
    DoTestOProdC(ca,v1,cv2,label+" C,R,C");
#endif
#endif
#endif
#endif

    if (showstartdone) std::cout<<"Done Test7"<<std::endl;
}

#undef ProductType
