
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV_Vec.h"
#include <fstream>
#include <cstdio>

#include "TMV_TestVectorArith.h"

template <class T> static void TestVectorReal()
{
    if (showstartdone) {
        std::cout<<"Start Test Real Vector"<<std::endl;
    }
    const int N = 100;

    tmv::Vector<T> v(N);

    for (int i=0; i<N; ++i) v(i) = T(i);

    for (int i=0; i<N; ++i) Assert(v(i) == T(i),"Setting Vector");

    tmv::VectorView<T,2> v2 = v.subVector(0,N,2);
    for (int i=0; i<N/2; ++i) Assert(v2(i) == T(2*i),
                                     "Reading Vector with stride = 2");

    for (int i=0; i<N/2; ++i) v2[i] = T(i + 1000);
    for (int i=0; i<N/2; ++i) Assert(v(2*i) == T(i+1000),
                                     "Writing Vector with stride = 2");

    tmv::Vector<T> v3 = v2;
    for (int i=0; i<N/2; ++i) Assert(v3[i] == v2[i],
                                     "Copying Vector with stride = 2");

    for (int i=0; i<N; ++i) v[i] = T(i);
    v.swap(2,5);
    Assert(v(2) == T(5) && v(5) == T(2),"Swapping elements of Vector");
    v.swap(2,5);
    Assert(v(2) == T(2) && v(5) == T(5),"Swapping elements of Vector");

    T sum = N*(N-1)/2;
    if (showacc) {
        std::cout<<"SumElements = "<<SumElements(v)<<std::endl;
        std::cout<<"expected "<<sum<<std::endl;
        std::cout<<"== "<<(SumElements(v)==sum)<<std::endl;
        std::cout<<"diff = "<<(SumElements(v)-sum)<<std::endl;
    }
    Assert(SumElements(v) == sum,"Vector SumElements(v)");

    v.reverseSelf();
    for (int i=0; i<N; ++i) Assert(v(i) == T(N-i-1),"Reversing Vector");

    for (int i=0; i<N; ++i) v(i) = T(i+10);
    v(23) = T(10*N);
    v(42) = T(0.25);
    v(15) = T(-20*N);
    int imax,imin;
    if (showacc) {
        std::cout<<"v = "<<v<<std::endl;
        std::cout<<"v.MaxAbs = "<<v.maxAbsElement(&imax)<<std::endl;
        std::cout<<"imax = "<<imax<<std::endl;
        std::cout<<"v.MinAbs = "<<v.minAbsElement(&imin)<<std::endl;
        std::cout<<"imin = "<<imin<<std::endl;
    }
    Assert(v.maxAbsElement(&imax) == T(20*N),
           "MaxAbsElement of Vector did not return correct value");
    Assert(imax == 15,
           "MaxAbsElement of Vector did not return correct index");
    Assert(v.minAbsElement(&imin) == T(0.25),
           "MinAbsElement of Vector did not return correct value");
    Assert(imin == 42,
           "MinAbsElement of Vector did not return correct index");
    Assert(v.maxElement(&imax) == T(10*N),
           "MaxElement of Vector did not return correct value");
    Assert(imax == 23,
           "MaxElement of Vector did not return correct index");
    Assert(v.minElement(&imin) == T(-20*N),
           "MinElement of Vector did not return correct value");
    Assert(imin == 15,
           "MinElement of Vector did not return correct index");

    tmv::Vector<T> a(N);
    tmv::Vector<T> b(N);
    for (int i=0; i<N; ++i) a(i) = T(3+i);

    b = a;
    for (int i=0; i<N; ++i) Assert(a(i) == b(i),"Vector1 = Vector2");

    Assert(a == b,"Testing Equality of Vectors");

    b(4) = T(0);
    Assert(a != b,"Vector = Vector copied address, not values");

    tmv::VectorF<T> af(N);
    for (int i=1; i<=N; ++i) af(i) = T(3+i-1);
    for (int i=1; i<=N; ++i) Assert(af(i) == a(i-1),"FortranStyle Vector access");
    tmv::ConstVectorViewF<T> afcv = af.view();
    for (int i=1; i<=N; ++i) Assert(afcv(i) == a(i-1),"FortranStyle Vector CV access");
    tmv::VectorViewF<T> afv = af.view();
    for (int i=1; i<=N; ++i) Assert(afv(i) == a(i-1),"FortranStyle Vector V access");
    Assert(a == af,"FortransStyle Vector = CStyle Vector");
    tmv::ConstVectorView<T> afcv_c = afcv;
    Assert(afcv_c == a,"CStyle View of FortransStyle Vector = CStyle Vector");
    Assert(afcv == a,"FortranStyle View of Vector == CStyle Vector");

    for (int i=0; i<N; ++i) b(i) = T(5+2*i);

    v = a+b;
    for (int i=0; i<N; ++i) Assert(v(i) == T(8+3*i),"Adding Vectors");

    v = a-b;
    if (showacc) {
        std::cout<<"a = "<<a<<std::endl;
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"a-b = "<<v<<std::endl;
    }
    for (int i=0; i<N; ++i) Assert(v(i) == T(-2-i),"Subtracting Vectors");

    // b(i) = 5+2i
    // a(i) = 3+i
    // a(i)*b(i) = 15+11i+2i^2
    // Sum = 15N + 11N(N-1)/2 + 2*N*(N-1)*(2N-1)/6
    T prod = 15*N + 11*N*(N-1)/2 + 2*N*(N-1)*(2*N-1)/6;
    if (showacc) {
        std::cout<<"a = "<<a<<std::endl;
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"a*b = "<<a*b<<std::endl;
        std::cout<<"prod = "<<prod<<std::endl;
    }
    Assert(a*b == prod,"Multiplying Vectors");

    tmv::Vector<T> c(5);
    c = v.subVector(10,70,12);
    for (int i=0; i<5; ++i) Assert(c(i) == v(10+12*i),"SubVector");

    if (tmv::Epsilon<T>() == T(0)) return;

    for(int i=0;i<N;++i) a(i) = T(i+10);
    for(int i=0;i<N;++i) b(i) = T(-3*i+191);

    prod = 2900;
    T normsum = tmv::TMV_SQRT(T(1373700));
    T normdiff = tmv::TMV_SQRT(T(1362100));
    if (showacc) {
        std::cout<<"a*b = "<<a*b<<std::endl;
        std::cout<<"expected prod = "<<prod<<std::endl;
        std::cout<<"abs(diff) = "<<tmv::TMV_ABS(a*b-prod)<<std::endl;
        std::cout<<"eps = "<<EPS*Norm(a)*Norm(b)<<std::endl;
        std::cout<<"a+b = "<<a+b<<std::endl;
        std::cout<<"Norm(a+b) = "<<Norm(a+b)<<std::endl;
        std::cout<<"expected normsum = "<<normsum<<std::endl;
        std::cout<<"abs(diff) = "<<tmv::TMV_ABS(Norm(a+b)-normsum)<<std::endl;
        std::cout<<"eps = "<<EPS*tmv::TMV_ABS(Norm1(a)+Norm1(b))<<std::endl;
        std::cout<<"Norm1(a) = "<<Norm1(a)<<std::endl;
        std::cout<<"Norm1(b) = "<<Norm1(b)<<std::endl;
    }
    Assert(tmv::TMV_ABS(a*b - prod) < EPS*Norm(a)*Norm(b),"Inner Product");
    Assert(tmv::TMV_ABS(Norm(a+b) - normsum) < EPS*tmv::TMV_ABS(Norm1(a)+Norm1(b)),"Vector Sum");
    Assert(tmv::TMV_ABS(Norm(a-b) - normdiff) < EPS*tmv::TMV_ABS(Norm1(a)+Norm1(b)),"Vector Diff");

    const int NN=20;
    tmv::Vector<T> w(NN);
    w << 3.3,1.2,5.4,-1.2,4.3,-9.4,0,-2,4,-11.5,
      -12,14,33,1,-9.3,-3.9,4.9,10,-31,1.e-33;

    tmv::Vector<T> origw = w;
    tmv::Vector<T> w2 = w;
    int perm[NN];
    if (showacc)
        std::cout<<"unsorted w = "<<w<<std::endl;

    w.sort(perm);
    if (showacc)
        std::cout<<"sorted w = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(w(i-1) <= w(i),"Sort real Vector");
    }
    w2.permute(perm);
    Assert(w2 == w,"Sort real Vector -- perm");
    w = origw;
    w.sort();
    Assert(w2 == w,"Sort real Vector -- without perm");
    w.reversePermute(perm);
    if (showacc)
        std::cout<<"reverse permute sorted Vector = "<<w<<std::endl;
    Assert(w == origw,"Reverse permute sorted Vector = orig");
    w2 = origw;

    w.sort(perm,tmv::Ascend,tmv::AbsComp);
    if (showacc)
        std::cout<<"sorted w abs = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(tmv::TMV_ABS(w(i-1)) <= tmv::TMV_ABS(w(i)),"Sort real Vector abs");
    }
    w2.permute(perm);
    if (showacc)
        std::cout<<"permuted w2 = "<<w<<std::endl;
    Assert(w2 == w,"Sort real Vector abs -- perm");
    w = origw;
    w.sort(tmv::Ascend,tmv::AbsComp);
    if (showacc)
        std::cout<<"sorted w abs (without perm) = "<<w<<std::endl;
    Assert(w2 == w,"Sort real Vector abs -- without perm");
    w = w2 = origw;

    w.sort(perm,tmv::Descend);
    if (showacc)
        std::cout<<"sorted w desc = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(w(i-1) >= w(i),"Sort real Vector desc");
    }
    w2.permute(perm);
    Assert(w2 == w,"Sort real Vector desc -- perm");
    w = origw;
    w.sort(tmv::Descend);
    Assert(w2 == w,"Sort real Vector desc -- without perm");
    w = w2 = origw;

    w.sort(perm,tmv::Descend,tmv::AbsComp);
    if (showacc)
        std::cout<<"sorted w desc abs = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(tmv::TMV_ABS(w(i-1)) >= tmv::TMV_ABS(w(i)),"Sort real Vector desc abs");
    }
    w2.permute(perm);
    Assert(w2 == w,"Sort real Vector desc abs -- perm");
    w = origw;
    w.sort(tmv::Descend,tmv::AbsComp);
    Assert(w2 == w,"Sort real Vector desc abs -- without perm");
    w = w2 = origw;

    if (showstartdone) {
        std::cout<<"Done Test Real Vector"<<std::endl;
    }
}

template <class T> static void TestVectorComplex()
{
    if (showstartdone) {
        std::cout<<"Start Test Complex Vector"<<std::endl;
    }
    const int N = 100;

    tmv::Vector<std::complex<T> > v(N);
    for (int i=0; i<N; ++i) v(i) = std::complex<T>(T(i),T(i+1234));

    for (int i=0; i<N; ++i) Assert(v(i).real() == T(i),
                                   "CVector set");
    for (int i=0; i<N; ++i) Assert(v(i).imag() == T(i+1234),
                                   "CVector set");

    tmv::VectorView<std::complex<T>,2> v1(v.subVector(0,N,2));
    for (int i=0; i<N/2; ++i) Assert(v1(i) == std::complex<T>(T(2*i),T(2*i+1234)),
                                     "CVector stride=2");

    for (int i=0; i<N/2; ++i) v1[i] = std::complex<T>(T(i),T(i+1234));
    for (int i=0; i<N/2; ++i) Assert(v[2*i] == std::complex<T>(T(i),T(i+1234)),
                                     "setting CVector with stride = 2");

    for (int i=0; i<N; ++i) v(i) = std::complex<T>(T(i),T(i+1234));

    v.swap(2,5);
    Assert(v[2] == std::complex<T>(5,5+1234),"Swap in CVector");
    Assert(v[5] == std::complex<T>(2,2+1234),"Swap in CVector");
    v.swap(2,5);

    tmv::Vector<std::complex<T> > v2 = v.conjugate();

    for (int i=0; i<N; ++i) 
        Assert(v2(i) == std::complex<T>(T(i),T(-i-1234)), "Conjugate CVector");
    Assert(v2 == v.conjugate(),"Conjugate == CVector");

    if (tmv::Epsilon<T>() == T(0)) return;

    std::complex<T> prod_act(0);
    for (int i=0; i<N; ++i) prod_act += v[i] * v2[i];
    std::complex<T> prod = v*v2;
    Assert(tmv::TMV_ABS(prod-prod_act) < EPS*tmv::TMV_ABS(prod_act),
           "CVector * CVector");
    prod = v*v.conjugate();
    Assert(tmv::TMV_ABS(prod.imag()) < EPS,"prod is real");
    Assert(tmv::TMV_ABS(prod-prod_act) < EPS*tmv::TMV_ABS(prod_act),
           "CVector * conj(CVector)");
    T norm1 = tmv::TMV_SQRT(prod.real());

    T norm2 = Norm(v);
    if (showacc) {
        std::cout<<"v = "<<v<<std::endl;
        std::cout<<"v2 = "<<v2<<std::endl;
        std::cout<<"v*v2 = "<<v*v2<<std::endl;
        std::cout<<"norm1 = "<<norm1<<std::endl;
        std::cout<<"norm2 = "<<norm2<<std::endl;
    }
    Assert(tmv::TMV_ABS(norm1 - norm2) < EPS*norm1,"Norm CVector");

    std::complex<T> sum_act(0);
    for (int i=0; i<N; ++i) sum_act += v[i];
    std::complex<T> sumel = v.sumElements();
    if (showacc) {
        std::cout<<"sumel = "<<sumel<<std::endl;
        std::cout<<"sumact = "<<sum_act<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS(sumel-sum_act)<<std::endl;
    }
    Assert(tmv::TMV_ABS(sumel-sum_act) < EPS*tmv::TMV_ABS(sum_act),"CVector SumElements");

    T sumabs_act(0);
    for (int i=0; i<N; ++i) sumabs_act += tmv::TMV_ABS(v[i]);
    T sumabsel = v.sumAbsElements();
    Assert(tmv::TMV_ABS(sumabsel-sumabs_act) < EPS*tmv::TMV_ABS(sumabs_act),
           "CVector SumAbsElements");

    v.conjugateSelf();
    Assert(v == v2,"ConjugateSelf CVector");
    v = v.conjugate();
    Assert(v == v2.conjugate(),"v = v.conjugate() CVector");

    tmv::Vector<T> a(N);
    for(int i=0;i<N;++i) a(i) = T(i+10);
    tmv::Vector<T> b(N);
    for(int i=0;i<N;++i) b(i) = T(-3*i+191);

    tmv::Vector<std::complex<T> > ca = a;
    Assert(Norm(ca-a) < EPS*Norm(a),"Copy real V -> complex V");

    ca *= std::complex<T>(3,4)/T(5);
    tmv::Vector<std::complex<T> > cb = b*std::complex<T>(3,4)/T(5);

    prod = T(29)*std::complex<T>(-28,96);
    T normsum = tmv::TMV_SQRT(T(1373700));
    T normdiff = tmv::TMV_SQRT(T(1362100));
    if (showacc) {
        std::cout<<"ca*cb = "<<ca*cb<<std::endl;
        std::cout<<"expected prod = "<<prod<<std::endl;
        std::cout<<"abs(diff) = "<<tmv::TMV_ABS(ca*cb-prod)<<std::endl;
        std::cout<<"eps = "<<EPS*Norm(ca)*Norm(cb)<<std::endl;
    }
    Assert(tmv::TMV_ABS(ca*cb - prod) < EPS*Norm(ca)*Norm(cb),"CInner Product");
    Assert(tmv::TMV_ABS(Norm(ca+cb) - normsum) < EPS*tmv::TMV_ABS(Norm(ca)+Norm(cb)),"CVector Sum");
    Assert(tmv::TMV_ABS(Norm(ca-cb) - normdiff) < EPS*tmv::TMV_ABS(Norm(ca)+Norm(cb)),"CVector Diff");

    const int NN=20;
    tmv::Vector<std::complex<T> > w(NN);
    w << 3.3,1.2,5.4,-1.2,4.3,-9.4,0,-2,4,-11.5,
      -12,14,33,1,-9.3,-3.9,4.9,10,-31,1.e-33;

    tmv::Vector<T> iw(NN);
    iw << 1.4,9.8,-0.2,-8.6,3.0,-4.4,3,9,-1.9,-11.4,
       11.1,-140,-23,11,5.2,-3.9,4.8,99,-71,-0.5;
    w.imagPart() = iw;

    tmv::Vector<std::complex<T> > origw = w;
    tmv::Vector<std::complex<T> > w2 = w;
    int perm[NN];
    if (showacc)
        std::cout<<"unsorted w = "<<w<<std::endl;

    w.sort(perm);
    if (showacc)
        std::cout<<"sorted w = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(real(w(i-1)) <= real(w(i)),"Sort complex Vector");
    }
    w2.permute(perm);
    Assert(w2 == w,"Sort complex Vector -- perm");
    w = origw;
    w.sort();
    Assert(w2 == w,"Sort complex Vector -- without perm");
    w.reversePermute(perm);
    if (showacc)
        std::cout<<"reverse permute sorted Vector = "<<w<<std::endl;
    Assert(w == origw,"Reverse permute sorted Vector = orig");
    w2 = origw;

    w.sort(perm,tmv::Ascend,tmv::AbsComp);
    if (showacc)
        std::cout<<"sorted w abs = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(tmv::TMV_ABS(w(i-1)) <= tmv::TMV_ABS(w(i)),"Sort complex Vector abs");
    }
    w2.permute(perm);
    Assert(w2 == w,"Sort complex Vector abs -- perm");
    w = origw;
    w.sort(tmv::Ascend,tmv::AbsComp);
    Assert(w2 == w,"Sort complex Vector abs -- without perm");
    w = w2 = origw;

    w.sort(perm,tmv::Ascend,tmv::ArgComp);
    if (showacc)
        std::cout<<"sorted w arg = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(tmv::TMV_ARG(w(i-1)) <= tmv::TMV_ARG(w(i)),"Sort complex Vector arg");
    }
    w2.permute(perm);
    Assert(w2 == w,"Sort complex Vector arg -- perm");
    w = origw;
    w.sort(tmv::Ascend,tmv::ArgComp);
    Assert(w2 == w,"Sort complex Vector arg -- without perm");
    w = w2 = origw;

    w.sort(perm,tmv::Ascend,tmv::ImagComp);
    if (showacc)
        std::cout<<"sorted w imag = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(imag(w(i-1)) <= imag(w(i)),"Sort complex Vector imag");
    }
    w2.permute(perm);
    Assert(w2 == w,"Sort complex Vector imag -- perm");
    w = origw;
    w.sort(tmv::Ascend,tmv::ImagComp);
    Assert(w2 == w,"Sort complex Vector imag -- without perm");
    w = w2 = origw;

    w.sort(perm,tmv::Descend);
    if (showacc)
        std::cout<<"sorted w desc = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(real(w(i-1)) >= real(w(i)),"Sort complex Vector desc");
    }
    w2.permute(perm);
    Assert(w2 == w,"Sort complex Vector desc -- perm");
    w = origw;
    w.sort(tmv::Descend);
    Assert(w2 == w,"Sort complex Vector desc -- without perm");
    w = w2 = origw;

    w.sort(perm,tmv::Descend,tmv::AbsComp);
    if (showacc)
        std::cout<<"sorted w desc abs = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(tmv::TMV_ABS(w(i-1)) >= tmv::TMV_ABS(w(i)),"Sort complex Vector desc abs");
    }
    w2.permute(perm);
    Assert(w2 == w,"Sort complex Vector desc abs -- perm");
    w = origw;
    w.sort(tmv::Descend,tmv::AbsComp);
    Assert(w2 == w,"Sort complex Vector desc abs -- without perm");
    w = w2 = origw;

    w.sort(perm,tmv::Descend,tmv::ArgComp);
    if (showacc)
        std::cout<<"sorted w desc arg = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(tmv::TMV_ARG(w(i-1)) >= tmv::TMV_ARG(w(i)),"Sort complex Vector desc arg");
    }
    w2.permute(perm);
    Assert(w2 == w,"Sort complex Vector desc arg -- perm");
    w = origw;
    w.sort(tmv::Descend,tmv::ArgComp);
    Assert(w2 == w,"Sort complex Vector desc arg -- without perm");
    w = w2 = origw;

    w.sort(perm,tmv::Descend,tmv::ImagComp);
    if (showacc)
        std::cout<<"sorted w desc imag = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(imag(w(i-1)) >= imag(w(i)),"Sort complex Vector desc imag");
    }
    w2.permute(perm);
    Assert(w2 == w,"Sort complex Vector desc imag -- perm");
    w = origw;
    w.sort(tmv::Descend,tmv::ImagComp);
    Assert(w2 == w,"Sort complex Vector desc imag -- without perm");
    w = w2 = origw;

    if (showstartdone) {
        std::cout<<"Done Test Complex Vector"<<std::endl;
    }
}

template <class T> static void TestVectorArith()
{
    typedef tmv::VectorView<T> V;
    typedef tmv::VectorView<std::complex<T> > CV;
    typedef tmv::VectorViewF<T> VF;
    typedef tmv::VectorViewF<std::complex<T> > CVF;

    if (showstartdone) {
        std::cout<<"Start Test Vector Arith"<<std::endl;
    }
    const int N = 100;

    tmv::Vector<T> a(N);
    for(int i=0;i<N;++i) a(i) = T(i+10);
    tmv::Vector<std::complex<T> > ca = a*std::complex<T>(2,-1);;
    V aa = a.view();
    CV caa = ca.view();
    TestVectorArith1<T>(aa,caa,"Vector C");

    tmv::Vector<T> b(N);
    for(int i=0;i<N;++i) b(i) = T(-3*i+2);
    tmv::Vector<std::complex<T> > cb = b*std::complex<T>(-5,1);
    V bb = b.view();
    CV cbb = cb.view();
    TestVectorArith2<T>(aa,caa,b,cbb,"Vector CC");

    tmv::Vector<T> a10(10*N);
    V as = a10.subVector(0,10*N,10);
    as = a;
    tmv::Vector<std::complex<T> > ca10(10*N);
    CV cas = ca10.subVector(0,10*N,10);
    cas = ca;
    TestVectorArith1<T>(as,cas,"Vector C Step");

    tmv::Vector<T> b10(10*N);
    V bs = b10.subVector(0,10*N,10);
    bs = b;
    tmv::Vector<std::complex<T> > cb10(10*N);
    CV cbs = cb10.subVector(0,10*N,10);
    cbs = cb;
    TestVectorArith2<T>(as,cas,bb,cbb,"Vector C StepA");
    TestVectorArith2<T>(aa,caa,bs,cbs,"Vector C StepB");
    TestVectorArith2<T>(as,cas,bs,cbs,"Vector C StepAB");

    V ar = aa.reverse();
    CV car = caa.reverse();
    TestVectorArith1<T>(ar,car,"Vector C Rev");

    V br = bb.reverse();
    CV cbr = cbb.reverse();
    TestVectorArith2<T>(ar,car,bb,cbb,"Vector C RevA");
    TestVectorArith2<T>(aa,caa,br,cbr,"Vector C RevB");
    TestVectorArith2<T>(ar,car,br,cbr,"Vector C RevAB");

#ifdef XTEST
#if (XTEST & 32)
    VF af = a.fView();
    CVF caf = ca.fView();
    TestVectorArith1<T>(af,caf,"Vector F");

#if 1
    VF bf = b.fView();
    CVF cbf = cb.fView();
    TestVectorArith2<T>(af,caf,b,cbb,"Vector FC");
    TestVectorArith2<T>(aa,caa,bf,cbf,"Vector CF");
    TestVectorArith2<T>(af,caf,bf,cbf,"Vector FF");

    VF asf = as;
    CVF casf = cas;
    TestVectorArith1<T>(asf,casf,"Vector F Step");

    VF bsf = bs;
    CVF cbsf = cbs;
    TestVectorArith2<T>(asf,casf,bf,cbf,"Vector F StepA");
    TestVectorArith2<T>(af,caf,bsf,cbsf,"Vector F StepB");
    TestVectorArith2<T>(asf,casf,bsf,cbsf,"Vector F StepAB");

    VF arf = ar;
    VF brf = br;
    CVF carf = car;
    CVF cbrf = cbr;
    TestVectorArith1<T>(arf,carf,"Vector F Rev");
    TestVectorArith2<T>(arf,carf,bf,cbf,"Vector F RevA");
    TestVectorArith2<T>(af,caf,brf,cbrf,"Vector F RevB");
    TestVectorArith2<T>(arf,carf,brf,cbrf,"Vector F RevAB");
#endif
#endif
#endif

    if (showstartdone) {
        std::cout<<"Done Test Vector Arith"<<std::endl;
    }
}

template <class T> static void TestVectorIO()
{
    if (showstartdone) {
        std::cout<<"Start Test Vector I/O"<<std::endl;
    }
    const int N = 20;
    tmv::Vector<T> v(N);
    tmv::Vector<std::complex<T> > cv(N);
    for (int i=0; i<N; ++i) {
        v(i) = T(i+34);
        cv(i) = std::complex<T>(T(i),T(N-i));
    }

    std::ofstream fout("tmvtest_vector_io.dat");
    if (!fout) throw std::runtime_error(
        "Couldn't open tmvtest_vector_io.dat for output");
    fout << v << std::endl << cv << std::endl;
    fout.close();

    tmv::Vector<T> xv1(N);
    tmv::Vector<std::complex<T> > xcv1(N);
    std::ifstream fin("tmvtest_vector_io.dat");
    if (!fin) throw std::runtime_error(
        "Couldn't open tmvtest_vector_io.dat for input");
    fin >> xv1 >> xcv1;
    fin.close();
    Assert(v == xv1,"Vector I/O check #1");
    Assert(cv == xcv1,"CVector I/O check #1");

    std::auto_ptr<tmv::Vector<T> > xv2;
    std::auto_ptr<tmv::Vector<std::complex<T> > > xcv2;
    fin.open("tmvtest_vector_io.dat");
    if (!fin) throw std::runtime_error(
        "Couldn't open tmvtest_vector_io.dat for input");
    fin >> xv2 >> xcv2;
    fin.close();
    Assert(v == *xv2,"Vector I/O check #2");
    Assert(cv == *xcv2,"CVector I/O check #2");

#ifndef XTEST
    std::remove("tmvtest_vector_io.dat");
#endif

    if (showstartdone) {
        std::cout<<"Done Test Vector I/O"<<std::endl;
    }
}

template <class T> void TestAllVector()
{
#if 1
    TestVectorReal<T>();
    TestVectorComplex<T>();
    TestVectorIO<T>();
#endif

#if 1
    TestVectorArith<T>();
#endif

    std::cout<<"Vector<"<<tmv::TMV_Text(T())<<"> passed all tests\n";
}

#ifdef TEST_DOUBLE
template void TestAllVector<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllVector<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllVector<long double>();
#endif
#ifdef TEST_INT
template void TestAllVector<int>();
#endif
