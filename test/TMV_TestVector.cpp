
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"
#include <fstream>
#include <cstdio>

#include "TMV_TestVectorArith.h"

template <class T> 
static void TestVectorReal()
{
    typedef typename tmv::Traits<T>::float_type FT;
    if (showstartdone) {
        std::cout<<"Start Test Real Vector"<<std::endl;
    }
    const int N = 100;

    tmv::Vector<T> v(N);
    Assert(int(v.size()) == N,"Creating Vector");
    for (int i=0; i<N; ++i) v(i) = T(i);
    for (int i=0; i<N; ++i) Assert(v(i) == T(i),"Setting Vector");

    tmv::VectorView<T> v2 = v.subVector(0,N,2);
    for (int i=0; i<N/2; ++i) 
        Assert(v2(i) == T(2*i),"Reading Vector with stride = 2");

    for (int i=0; i<N/2; ++i) v2[i] = T(i + 1000);
    for (int i=0; i<N/2; ++i) 
        Assert(v(2*i) == T(i+1000),"Writing Vector with stride = 2");

    tmv::Vector<T> v3 = v2;
    for (int i=0; i<N/2; ++i) 
        Assert(v3[i] == v2[i],"Copying Vector with stride = 2");

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
    Assert(Equal2(v.maxAbsElement(&imax),T(20*N),EPS),
           "MaxAbsElement of Vector did not return correct value");
    Assert(imax == 15,
           "MaxAbsElement of Vector did not return correct index");
    Assert(Equal2(v.minAbsElement(&imin),T(0.25),EPS),
           "MinAbsElement of Vector did not return correct value");
    Assert(imin == 42,
           "MinAbsElement of Vector did not return correct index");
    Assert(Equal2(v.maxElement(&imax),T(10*N),EPS),
           "MaxElement of Vector did not return correct value");
    Assert(imax == 23,
           "MaxElement of Vector did not return correct index");
    Assert(Equal2(v.minElement(&imin),T(-20*N),EPS),
           "MinElement of Vector did not return correct value");
    Assert(imin == 15,
           "MinElement of Vector did not return correct index");
    Assert(Equal2(v.maxAbs2Element(&imax),T(20*N),EPS),
           "MaxAbs2Element of Vector did not return correct value");
    Assert(imax == 15,
           "MaxAbs2Element of Vector did not return correct index");
    Assert(Equal2(v.minAbs2Element(&imin),T(0.25),EPS),
           "MinAbs2Element of Vector did not return correct value");
    Assert(imin == 42,
           "MinAbs2Element of Vector did not return correct index");

    tmv::Vector<T> a(N);
    tmv::Vector<T> b(N);
    for (int i=0; i<N; ++i) a(i) = T(3+i);

    b = a;
    for (int i=0; i<N; ++i) Assert(a(i) == b(i),"Vector1 = Vector2");

    Assert(a == b,"Testing Equality of Vectors");

    b(4) = T(0);
    Assert(a != b,"Vector = Vector copied address, not values");

    tmv::Vector<T,tmv::FortranStyle> af(N);
    for (int i=1; i<=N; ++i) af(i) = T(3+i-1);
    for (int i=1; i<=N; ++i) 
        Assert(af(i) == a(i-1),"FortranStyle Vector access");
    tmv::ConstVectorViewF<T> afcv = af.view();
    for (int i=1; i<=N; ++i) 
        Assert(afcv(i) == a(i-1),"FortranStyle Vector CV access");
    tmv::VectorViewF<T> afv = af.view();
    for (int i=1; i<=N; ++i) 
        Assert(afv(i) == a(i-1),"FortranStyle Vector V access");
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

    for(int i=0;i<N;++i) a(i) = T(i+10);
    for(int i=0;i<N;++i) b(i) = T(-3*i+191);

    prod = 2900;
    T normsqsum = 1373700;
    T normsqdiff = 1362100;
    if (showacc) {
        std::cout<<"a*b = "<<a*b<<std::endl;
        std::cout<<"expected prod = "<<prod<<std::endl;
        std::cout<<"abs(diff) = "<<tmv::TMV_ABS(a*b-prod)<<std::endl;
        std::cout<<"eps = "<<EPS<<std::endl;
        std::cout<<"a+b = "<<a+b<<std::endl;
        std::cout<<"NormSq(a+b) = "<<NormSq(a+b)<<std::endl;
        std::cout<<"expected normsqsum = "<<normsqsum<<std::endl;
        std::cout<<"abs(diff) = "<<tmv::TMV_ABS(NormSq(a+b)-normsqsum)<<std::endl;
        std::cout<<"eps = "<<EPS*tmv::TMV_ABS(Norm1(a)+Norm1(b))<<std::endl;
        std::cout<<"Norm1(a) = "<<Norm1(a)<<std::endl;
        std::cout<<"Norm1(b) = "<<Norm1(b)<<std::endl;
    }
    FT eps = EPS;
    if (!std::numeric_limits<T>::is_integer) eps *= Norm(a) * Norm(b);
    Assert(Equal2(a*b,prod,eps),"Inner Product");
    Assert(Equal2(NormSq(a+b),normsqsum,EPS*normsqsum),"Vector Sum");
    Assert(Equal2(NormSq(a-b),normsqdiff,EPS*normsqdiff),"Vector Diff");

    const int NN=20;
    tmv::Vector<T> w(NN);
    w << 33,12,54,-12,43,-94,0,-20,40,-115,
      -120,140,330,10,-93,-39,49,100,-310,1;

    tmv::Vector<T> origw = w;
    tmv::Vector<T> w2 = w;
    tmv::Permutation P(NN);
    if (showacc) std::cout<<"unsorted w = "<<w<<std::endl;

    w.sort(P);
    if (showacc) std::cout<<"sorted w = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(w(i-1) <= w(i),"Sort real Vector");
    }
    w2 = P * w2;
    Assert(w2 == w,"Sort real Vector -- perm");
    w = origw;
    w.sort();
    Assert(w2 == w,"Sort real Vector -- without perm");
    w = P.inverse() * w;
    if (showacc) std::cout<<"reverse permute sorted Vector = "<<w<<std::endl;
    Assert(w == origw,"Reverse permute sorted Vector = orig");
    w2 = origw;

    w.sort(P,tmv::Ascend,tmv::AbsComp);
    if (showacc) std::cout<<"sorted w abs = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(tmv::TMV_ABS(w(i-1)) <= tmv::TMV_ABS(w(i)),
               "Sort real Vector abs");
    }
    w2 = P * w2;
    if (showacc) std::cout<<"permuted w2 = "<<w<<std::endl;
    Assert(w2 == w,"Sort real Vector abs -- perm");
    w = origw;
    w.sort(tmv::Ascend,tmv::AbsComp);
    if (showacc)
        std::cout<<"sorted w abs (without perm) = "<<w<<std::endl;
    Assert(w2 == w,"Sort real Vector abs -- without perm");
    w = w2 = origw;

    w.sort(P,tmv::Descend);
    if (showacc) std::cout<<"sorted w desc = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(w(i-1) >= w(i),"Sort real Vector desc");
    }
    w2 = P * w2;
    Assert(w2 == w,"Sort real Vector desc -- perm");
    w = origw;
    w.sort(tmv::Descend);
    Assert(w2 == w,"Sort real Vector desc -- without perm");
    w = w2 = origw;

    w.sort(P,tmv::Descend,tmv::AbsComp);
    if (showacc) std::cout<<"sorted w desc abs = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(tmv::TMV_ABS(w(i-1)) >= tmv::TMV_ABS(w(i)),
               "Sort real Vector desc abs");
    }
    w2 = P * w2;
    Assert(w2 == w,"Sort real Vector desc abs -- perm");
    w = origw;
    w.sort(tmv::Descend,tmv::AbsComp);
    Assert(w2 == w,"Sort real Vector desc abs -- without perm");
    w = w2 = origw;

    v.resize(2);
    Assert(v.size() == 2,"v.resize(2)");
    for (int i=0; i<2; ++i) v(i) = T(i);
    for (int i=0; i<2; ++i) Assert(v(i) == T(i),"Setting resized Vector");

    v.resize(2*N);
    Assert(int(v.size()) == 2*N,"v.resize(2*N)");
    for (int i=0; i<2*N; ++i) v(i) = T(i);
    for (int i=0; i<2*N; ++i) Assert(v(i) == T(i),"Setting resized Vector");

    if (showstartdone) {
        std::cout<<"Done Test Real Vector"<<std::endl;
    }
}

template <class T> 
static void TestVectorComplex()
{
    typedef typename tmv::Traits<T>::float_type FT;
    if (showstartdone) {
        std::cout<<"Start Test Complex Vector"<<std::endl;
    }
    const int N = 100;

    tmv::Vector<std::complex<T> > v(N);
    for (int i=0; i<N; ++i) v(i) = std::complex<T>(T(i),T(i+1234));

    for (int i=0; i<N; ++i) 
        Assert(v(i).real() == T(i), "CVector set");
    for (int i=0; i<N; ++i) 
        Assert(v(i).imag() == T(i+1234), "CVector set");

    tmv::VectorView<std::complex<T> > v1(v.subVector(0,N,2));
    for (int i=0; i<N/2; ++i) 
        Assert(v1(i) == std::complex<T>(T(2*i),T(2*i+1234)),
               "CVector stride=2");

    for (int i=0; i<N/2; ++i) v1[i] = std::complex<T>(T(i),T(i+1234));
    for (int i=0; i<N/2; ++i) 
        Assert(v[2*i] == std::complex<T>(T(i),T(i+1234)),
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

    tmv::Vector<std::complex<T> > v3(N);
    for (int i=0; i<N; ++i) v3(i) = std::complex<T>(i+10,2*i);
    v3(23) = std::complex<T>(40*N,9*N);
    v3(42) = std::complex<T>(0,1);
    v3(15) = std::complex<T>(-32*N,24*N);
    int imax,imin;
    if (showacc) {
        std::cout<<"v = "<<v3<<std::endl;
        std::cout<<"v.MaxAbs = "<<v3.maxAbsElement(&imax)<<std::endl;
        std::cout<<"imax = "<<imax<<std::endl;
        std::cout<<"v.MinAbs = "<<v3.minAbsElement(&imin)<<std::endl;
        std::cout<<"imin = "<<imin<<std::endl;
    }

    Assert(Equal2(v3.maxAbsElement(&imax),FT(41*N),EPS),
           "MaxAbsElement of complex Vector did not return correct value");
    Assert(imax == 23,
           "MaxAbsElement of complex Vector did not return correct index");
    Assert(Equal2(v3.minAbsElement(&imin),FT(1),EPS),
           "MinAbsElement of complex Vector did not return correct value");
    Assert(imin == 42,
           "MinAbsElement of complex Vector did not return correct index");
    Assert(Equal2(v3.maxAbs2Element(&imax),T(56*N),EPS),
           "MaxAbs2Element of complex Vector did not return correct value");
    Assert(imax == 15,
           "MaxAbs2Element of complex Vector did not return correct index");
    Assert(Equal2(v3.minAbs2Element(&imin),T(1),EPS),
           "MinAbs2Element of complex Vector did not return correct value");
    Assert(imin == 42,
           "MinAbs2Element of complex Vector did not return correct index");

    std::complex<T> prod_act(0);
    for (int i=0; i<N; ++i) prod_act += v[i] * v2[i];
    std::complex<T> prod = v*v2;
    Assert(Equal2(prod,prod_act,EPS*tmv::TMV_ABS2(prod_act)),
           "CVector * CVector");
    prod = v*v.conjugate();
    prod_act = T(0);
    for (int i=0; i<N; ++i) prod_act += v[i] * std::conj(v[i]);
    Assert(Equal2(prod.imag(),T(0),EPS),"prod is real");
    Assert(Equal2(prod,prod_act,EPS*tmv::TMV_ABS2(prod_act)),
           "CVector * conj(CVector)");

    FT norm1 = tmv::TMV_SQRT(prod.real());
    FT norm2 = Norm(v);
    if (showacc) {
        std::cout<<"v = "<<v<<std::endl;
        std::cout<<"v2 = "<<v2<<std::endl;
        std::cout<<"v*v2 = "<<v*v2<<std::endl;
        std::cout<<"norm1 = "<<norm1<<std::endl;
        std::cout<<"norm2 = "<<norm2<<std::endl;
    }
    Assert(Equal2(norm1,norm2,EPS*norm1),"Norm CVector");

    std::complex<T> sum_act(0);
    for (int i=0; i<N; ++i) sum_act += v[i];
    std::complex<T> sumel = v.sumElements();
    if (showacc) {
        std::cout<<"sumel = "<<sumel<<std::endl;
        std::cout<<"sumact = "<<sum_act<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS(sumel-sum_act)<<std::endl;
    }
    Assert(Equal2(sumel,sum_act,EPS*tmv::TMV_ABS2(sum_act)),
           "CVector SumElements");

    FT sumabs_act(0);
    for (int i=0; i<N; ++i) sumabs_act += tmv::TMV_ABS(v[i]);
    FT sumabsel = v.sumAbsElements();
    if (showacc) {
        std::cout<<"sumabsel = "<<sumabsel<<std::endl;
        std::cout<<"sumabs_act = "<<sumabs_act<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS(sumabsel-sumabs_act)<<std::endl;
    }
    Assert(Equal2(sumabsel,sumabs_act,EPS*tmv::TMV_ABS2(sumabs_act)),
           "CVector SumAbsElements");
    T sumabs2_act(0);
    for (int i=0; i<N; ++i) sumabs2_act += tmv::TMV_ABS2(v[i]);
    T sumabs2el = v.sumAbs2Elements();
    Assert(Equal2(sumabs2el,sumabs2_act,EPS*tmv::TMV_ABS2(sumabs2_act)),
           "CVector SumAbs2Elements");

    v.conjugateSelf();
    Assert(v == v2,"ConjugateSelf CVector");
    v = v.conjugate();
    Assert(v == v2.conjugate(),"v = v.conjugate() CVector");

    tmv::Vector<T> a(N);
    for(int i=0;i<N;++i) a(i) = T(i+10);
    tmv::Vector<T> b(N);
    for(int i=0;i<N;++i) b(i) = T(-3*i+191);

    tmv::Vector<std::complex<T> > ca = a;
    Assert(Equal(ca,a,EPS),"Copy real V -> complex V");

    ca *= std::complex<T>(3,4);
    tmv::Vector<std::complex<T> > cb = b*std::complex<T>(3,4);

    prod = T(29)*T(25)*std::complex<T>(-28,96);
    T normsqsum = 34342500;
    T normsqdiff = 34052500;
    if (showacc) {
        std::cout<<"ca*cb = "<<ca*cb<<std::endl;
        std::cout<<"expected prod = "<<prod<<std::endl;
        std::cout<<"abs(diff) = "<<tmv::TMV_ABS(ca*cb-prod)<<std::endl;
        std::cout<<"eps = "<<EPS*Norm(ca)*Norm(cb)<<std::endl;
        std::cout<<"NormSq(ca+cb) = "<<NormSq(ca+cb)<<std::endl;
        std::cout<<"expected sum = "<<normsqsum<<std::endl;
        std::cout<<"abs(diff) = "<<tmv::TMV_ABS2(NormSq(ca+cb)-normsqsum)<<std::endl;
        std::cout<<"eps = "<<EPS*normsqsum<<std::endl;
        std::cout<<"NormSq(ca-cb) = "<<NormSq(ca-cb)<<std::endl;
        std::cout<<"expected diff = "<<normsqdiff<<std::endl;
        std::cout<<"abs(diff) = "<<tmv::TMV_ABS2(NormSq(ca+cb)-normsqdiff)<<std::endl;
        std::cout<<"eps = "<<EPS*normsqdiff<<std::endl;
    }
    FT eps = EPS;
    if (!std::numeric_limits<T>::is_integer) eps *= Norm(ca) * Norm(cb);
    Assert(Equal2(ca*cb,prod,eps),"CInner Product");
    Assert(Equal2(NormSq(ca+cb),normsqsum,EPS*normsqsum),"CVector Sum");
    Assert(Equal2(NormSq(ca-cb),normsqdiff,EPS*normsqdiff),"CVector Diff");

    const int NN=20;
    tmv::Vector<std::complex<T> > w(NN);
    w << 33,12,54,-12,43,-94,0,-20,40,-115,
      -120,140,330,10,-93,-39,49,100,-310,1;

    tmv::Vector<T> iw(NN);
    iw << 14,98,-2,-86,30,-44,30,90,-19,-114,
       111,-1400,-230,110,52,-39,48,990,-710,-5;
    w.imagPart() = iw;

    tmv::Vector<std::complex<T> > origw = w;
    tmv::Vector<std::complex<T> > w2 = w;
    tmv::Permutation P(NN);
    if (showacc) std::cout<<"unsorted w = "<<w<<std::endl;

    w.sort(P);
    if (showacc) std::cout<<"sorted w = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(real(w(i-1)) <= real(w(i)),"Sort complex Vector");
    }
    w2 = P * w2;
    Assert(w2 == w,"Sort complex Vector -- perm");
    w = origw;
    w.sort();
    Assert(w2 == w,"Sort complex Vector -- without perm");
    w = P.inverse() * w;
    if (showacc) std::cout<<"reverse permute sorted Vector = "<<w<<std::endl;
    Assert(w == origw,"Reverse permute sorted Vector = orig");
    w2 = origw;

    w.sort(P,tmv::Ascend,tmv::AbsComp);
    if (showacc) std::cout<<"sorted w abs = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(tmv::TMV_ABS(w(i-1)) <= tmv::TMV_ABS(w(i)),
               "Sort complex Vector abs");
    }
    w2 = P * w2;
    Assert(w2 == w,"Sort complex Vector abs -- perm");
    w = origw;
    w.sort(tmv::Ascend,tmv::AbsComp);
    Assert(w2 == w,"Sort complex Vector abs -- without perm");
    w = w2 = origw;

    w.sort(P,tmv::Ascend,tmv::ArgComp);
    if (showacc) std::cout<<"sorted w arg = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(tmv::TMV_ARG(w(i-1)) <= tmv::TMV_ARG(w(i)),
               "Sort complex Vector arg");
    }
    w2 = P * w2;
    Assert(w2 == w,"Sort complex Vector arg -- perm");
    w = origw;
    w.sort(tmv::Ascend,tmv::ArgComp);
    Assert(w2 == w,"Sort complex Vector arg -- without perm");
    w = w2 = origw;

    w.sort(P,tmv::Ascend,tmv::ImagComp);
    if (showacc) std::cout<<"sorted w imag = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(imag(w(i-1)) <= imag(w(i)),"Sort complex Vector imag");
    }
    w2 = P * w2;
    Assert(w2 == w,"Sort complex Vector imag -- perm");
    w = origw;
    w.sort(tmv::Ascend,tmv::ImagComp);
    Assert(w2 == w,"Sort complex Vector imag -- without perm");
    w = w2 = origw;

    w.sort(P,tmv::Descend);
    if (showacc) std::cout<<"sorted w desc = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(real(w(i-1)) >= real(w(i)),"Sort complex Vector desc");
    }
    w2 = P * w2;
    Assert(w2 == w,"Sort complex Vector desc -- perm");
    w = origw;
    w.sort(tmv::Descend);
    Assert(w2 == w,"Sort complex Vector desc -- without perm");
    w = w2 = origw;

    w.sort(P,tmv::Descend,tmv::AbsComp);
    if (showacc) std::cout<<"sorted w desc abs = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(tmv::TMV_ABS(w(i-1)) >= tmv::TMV_ABS(w(i)),
               "Sort complex Vector desc abs");
    }
    w2 = P * w2;
    Assert(w2 == w,"Sort complex Vector desc abs -- perm");
    w = origw;
    w.sort(tmv::Descend,tmv::AbsComp);
    Assert(w2 == w,"Sort complex Vector desc abs -- without perm");
    w = w2 = origw;

    w.sort(P,tmv::Descend,tmv::ArgComp);
    if (showacc) std::cout<<"sorted w desc arg = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(tmv::TMV_ARG(w(i-1)) >= tmv::TMV_ARG(w(i)),
               "Sort complex Vector desc arg");
    }
    w2 = P * w2;
    Assert(w2 == w,"Sort complex Vector desc arg -- perm");
    w = origw;
    w.sort(tmv::Descend,tmv::ArgComp);
    Assert(w2 == w,"Sort complex Vector desc arg -- without perm");
    w = w2 = origw;

    w.sort(P,tmv::Descend,tmv::ImagComp);
    if (showacc)
        std::cout<<"sorted w desc imag = "<<w<<std::endl;
    for(int i=1;i<NN;++i) {
        Assert(imag(w(i-1)) >= imag(w(i)),"Sort complex Vector desc imag");
    }
    w2 = P * w2;
    Assert(w2 == w,"Sort complex Vector desc imag -- perm");
    w = origw;
    w.sort(tmv::Descend,tmv::ImagComp);
    Assert(w2 == w,"Sort complex Vector desc imag -- without perm");
    w = w2 = origw;

    if (showstartdone) {
        std::cout<<"Done Test Complex Vector"<<std::endl;
    }
}

template <class T> 
static void TestVectorArith()
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

#if (XTEST & 32)
    VF af = a.fView();
    CVF caf = ca.fView();
    TestVectorArith1<T>(af,caf,"Vector F");

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

    if (showstartdone) {
        std::cout<<"Done Test Vector Arith"<<std::endl;
    }
}

template <class T> 
static void TestVectorIO()
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

    std::remove("tmvtest_vector_io.dat");

    if (showstartdone) {
        std::cout<<"Done Test Vector I/O"<<std::endl;
    }
}

template <class T> void TestVector()
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
template void TestVector<double>();
#endif
#ifdef TEST_FLOAT
template void TestVector<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestVector<long double>();
#endif
#ifdef TEST_INT
template void TestVector<int>();
#endif
