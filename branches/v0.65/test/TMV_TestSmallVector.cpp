
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_SmallVectorArith.h"
#include <fstream>
#include <cstdio>
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV_TestVectorArith.h"

#define N 100
#define NN 20

template <class T> 
static void TestSmallVectorReal()
{
    tmv::SmallVector<T,N> v;

    for (int i=0; i<N; ++i) v(i) = T(i);

    for (int i=0; i<N; ++i) Assert(v(i) == T(i),"Setting SmallVector");

    tmv::VectorView<T> v1 = v.subVector(0,N,2);
    for (int i=0; i<N/2; ++i) Assert(v1(i) == T(2*i), "SmallVector stride=2");

    for (int i=0; i<N/2; ++i) v1[i] = T(i+1234);
    for (int i=0; i<N/2; ++i) Assert(v[2*i] == T(i+1234),
                                     "setting SmallVector with stride = 2");
    for (int i=0; i<N; ++i) v(i) = T(i);

    v.swap(2,5);
    Assert(v(2) == T(5) && v(5) == T(2),"Swapping elements of SmallVector");
    v.swap(2,5);
    Assert(v(2) == T(2) && v(5) == T(5),"Swapping elements of SmallVector");

    T sum = N*(N-1)/2;
    Assert(SumElements(v) == sum,"SmallVector SumElements(v)");

    v.reverseSelf();
    for (int i=0; i<N; ++i) Assert(v(i) == T(N-i-1),"Reversing SmallVector");

    for (int i=0; i<N; ++i) v(i) = T(i+10);
    v(23) = T(10*N);
    v(42) = T(1)/T(4);
    v(15) = T(-20*N);
    int imax,imin;
    Assert(v.maxAbsElement(&imax) == T(20*N),
           "MaxAbsElement of SmallVector did not return correct value");
    Assert(imax == 15,
           "MaxAbsElement of SmallVector did not return correct index");
    Assert(v.minAbsElement(&imin) == T(1)/T(4),
           "MinAbsElement of SmallVector did not return correct value");
    Assert(imin == 42,
           "MinAbsElement of SmallVector did not return correct index");
    Assert(v.maxElement(&imax) == T(10*N),
           "MaxElement of SmallVector did not return correct value");
    Assert(imax == 23,
           "MaxElement of SmallVector did not return correct index");
    Assert(v.minElement(&imin) == T(-20*N),
           "MinElement of SmallVector did not return correct value");
    Assert(imin == 15,
           "MinElement of SmallVector did not return correct index");

    tmv::SmallVector<T,N> a;
    tmv::SmallVector<T,N> b;
    for (int i=0; i<N; ++i) a(i) = T(3+i);

    b = a;
    for (int i=0; i<N; ++i) Assert(a(i) == b(i),"SmallVector1 = SmallVector2");

    Assert(a == b,"Testing Equality of SmallVectors");

    b(4) = 0;
    Assert(a != b,"SmallVector = SmallVector copied address, not values");

    tmv::SmallVector<T,N,tmv::FortranStyle> af;
    for (int i=1; i<=N; ++i) af(i) = T(3+i-1);
    for (int i=1; i<=N; ++i) Assert(af(i) == a(i-1),
                                    "FortranStyle SmallVector access");
    Assert(a == af,"FortransStyle SmallVector = CStyle SmallVector");

    for (int i=0; i<N; ++i) b(i) = T(5+2*i);

    v = a+b;
    for (int i=0; i<N; ++i) Assert(v(i) == T(8+3*i),"Adding SmallVectors");

    v = a-b;
    tmv::Vector<T> vv = a-b;
    for (int i=0; i<N; ++i) Assert(v(i) == T(-2-i),"Subtracting SmallVectors");

    a(0) = 1;
    b(0) = 1024;
    for (int i=1; i<10; ++i) {
        a(i) = a(i-1)*T(2);
        b(i) = b(i-1)/T(2);
    }
    for (int i=10; i<N; ++i) a(i) = 0;

    if (showacc) {
        std::cout<<"a = "<<a<<std::endl;
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"a*b = "<<a*b<<std::endl;
        std::cout<<"a*b-10240 = "<<(a*b-T(10240))<<std::endl;
        //ProductType(T,T) res(0);
        //for(int i=0;i<N;++i) res += a[i]*a[i]; 
        //std::cout<<"res = "<<res<<std::endl;
        //std::cout<<"MultVV = "<<tmv::MultVV<N>(a.cptr(),b.cptr())<<std::endl;
    }
    Assert(a*b == T(10240),"Multiplying SmallVectors");

    tmv::SmallVector<T,5> c = v.subVector(10,70,12);
    for (int i=0; i<5; ++i) Assert(c(i) == v(10+12*i),"SubSmallVector");

    for(int i=0;i<N;++i) a(i) = T(i+10);
    for(int i=0;i<N;++i) b(i) = T(-3*i+191);

    T prod = 2900;
    T normsqsum = 1373700;
    T normsqdiff = 1362100;
    T eps = EPS;
    if (!std::numeric_limits<T>::is_integer) eps *= Norm(a) * Norm(b);
    Assert(Equal2(a*b,prod,eps),"Inner Product");
    tmv::SmallVector<T,N> temp;
    T eps2 = EPS * tmv::TMV_ABS2(Norm1(a)+Norm1(b));
    Assert(Equal2(NormSq(temp=a+b),normsqsum,eps2),"SmallVector Sum");
    Assert(Equal2(NormSq(temp=a-b),normsqdiff,eps2),"SmallVector Diff");

    tmv::SmallVector<T,NN> w;
    w <<
        33,12,54,-12,43,-94,0,-20,40,-115,
        -120,140,330,10,-93,-39,49,100,-310,1;

    tmv::SmallVector<T,NN> origw = w;
    int perm[NN];

    if (showacc)
        std::cout<<"unsorted w = "<<w<<std::endl;
    w.sort(perm);
    for(int i=1;i<NN;++i) {
        Assert(w(i-1) <= w(i),"Sort real SmallVector");
    }
    if (showacc)
        std::cout<<"sorted w = "<<w<<std::endl;

    w.sort(0);
    w.reversePermute(perm);
    if (showacc)
        std::cout<<"Reverse permuted w = "<<w<<std::endl;
    Assert(w==origw,"Reverse permute sorted SmallVector = orig");
    w.sort(0);
    origw.permute(perm);
    if (showacc)
        std::cout<<"Sort permuted w = "<<origw<<std::endl;
    Assert(w==origw,"Permute SmallVector = sorted SmallVector");
}

template <class T> 
static void TestSmallVectorComplex()
{
    tmv::SmallVector<std::complex<T>,N> v;
    for (int i=0; i<N; ++i) v(i) = std::complex<T>(T(i),T(i+1234));

    for (int i=0; i<N; ++i) 
        Assert(real(v(i)) == T(i),"CSmallVector set");
    for (int i=0; i<N; ++i) 
        Assert(imag(v(i)) == T(i+1234),"CSmallVector set");

    tmv::VectorView<std::complex<T> > v1 = v.subVector(0,N,2);
    for (int i=0; i<N/2; ++i) 
        Assert(v1(i)==std::complex<T>(T(2*i),T(2*i+1234)),
               "CSmallVector stride=2");

    for (int i=0; i<N/2; ++i) v1[i] = std::complex<T>(T(i),T(i+9876));
    for (int i=0; i<N/2; ++i) 
        Assert(v[2*i]==std::complex<T>(T(i),T(i+9876)),
               "setting CSmallVector with stride = 2");

    for (int i=0; i<N; ++i) v(i) = std::complex<T>(T(i),T(i+1234));

    v.swap(2,5);
    Assert(v[2] == std::complex<T>(5,5+1234),"Swap in CSmallVector");
    Assert(v[5] == std::complex<T>(2,2+1234),"Swap in CSmallVector");
    v.swap(2,5);

    tmv::SmallVector<std::complex<T>,N> v2 = v.conjugate();

    for (int i=0; i<N; ++i) 
        Assert(v2(i) == std::complex<T>(T(i),T(-i-1234)),
               "Conjugate CSmallVector");
    Assert(v2 == v.conjugate(),"Conjugate == CSmallVector");

    tmv::SmallVector<std::complex<T>,N> v3;
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

    if (!std::numeric_limits<T>::is_integer) {
        Assert(Equal2(v3.maxAbsElement(&imax),T(41*N),EPS),
               "MaxAbsElement of Vector did not return correct value");
        Assert(imax == 23,
               "MaxAbsElement of Vector did not return correct index");
        Assert(Equal2(v3.minAbsElement(&imin),T(1),EPS),
               "MinAbsElement of Vector did not return correct value");
        Assert(imin == 42,
               "MinAbsElement of Vector did not return correct index");
    }
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

    if (!std::numeric_limits<T>::is_integer) {
        T norm1 = tmv::TMV_SQRT(prod.real());
        T norm2 = Norm(v);
        if (showacc) {
            std::cout<<"v = "<<v<<std::endl;
            std::cout<<"v2 = "<<v2<<std::endl;
            std::cout<<"v*v2 = "<<v*v2<<std::endl;
            std::cout<<"norm1 = "<<norm1<<std::endl;
            std::cout<<"norm2 = "<<norm2<<std::endl;
        }
        Assert(Equal2(norm1,norm2,EPS*norm1),"Norm CVector");
    }


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

    if (!std::numeric_limits<T>::is_integer) {
        T sumabs_act(0);
        for (int i=0; i<N; ++i) sumabs_act += tmv::TMV_ABS(v[i]);
        T sumabsel = v.sumAbsElements();
        Assert(Equal2(sumabsel,sumabs_act,EPS*tmv::TMV_ABS2(sumabs_act)),
               "CVector SumAbsElements");
    }
    T sumabs2_act(0);
    for (int i=0; i<N; ++i) sumabs2_act += tmv::TMV_ABS2(v[i]);
    T sumabs2el = v.sumAbs2Elements();
    Assert(Equal2(sumabs2el,sumabs2_act,EPS*tmv::TMV_ABS2(sumabs2_act)),
           "CVector SumAbs2Elements");

    v.conjugateSelf();
    Assert(v == v2,"ConjugateSelf CVector");
    v = v.conjugate();
    Assert(v == v2.conjugate(),"v = v.conjugate() CVector");

    tmv::SmallVector<T,N> a;
    for(int i=0;i<N;++i) a(i) = T(i+10);
    tmv::SmallVector<T,N> b;
    for(int i=0;i<N;++i) b(i) = T(-3*i+191);

    tmv::SmallVector<std::complex<T>,N> ca = a;
    Assert(Equal(ca,a,EPS),"Copy real V -> complex V");

    ca *= std::complex<T>(3,4);
    tmv::SmallVector<std::complex<T>,N> cb = b*std::complex<T>(3,4);
    prod = T(29)*T(25)*std::complex<T>(-28,96);
    T normsqsum = 34342500;
    T normsqdiff = 34052500;
    if (showacc) {
        std::cout<<"ca*cb = "<<ca*cb<<std::endl;
        std::cout<<"expected prod = "<<prod<<std::endl;
        std::cout<<"abs(diff) = "<<tmv::TMV_ABS(ca*cb-prod)<<std::endl;
        std::cout<<"eps = "<<EPS*Norm(ca)*Norm(cb)<<std::endl;
    }
    T eps = EPS;
    if (!std::numeric_limits<T>::is_integer) eps *= Norm(ca) * Norm(cb);
    Assert(Equal2(ca*cb,prod,eps),"CInner Product");
    T eps2 = EPS;
    if (!std::numeric_limits<T>::is_integer) eps2 *= Norm1(ca) * Norm1(cb);
    Assert(Equal2(NormSq(ca+cb),normsqsum,eps2),"CVector Sum");
    Assert(Equal2(NormSq(ca-cb),normsqdiff,eps2),"CVector Diff");

    tmv::SmallVector<std::complex<T>,NN> w;
    w <<
        33,12,54,-12,43,-94,0,-20,40,-115,
        -120,140,330,10,-93,-39,49,100,-310,1;
    tmv::SmallVector<T,NN> iw;
    iw <<
        14,98,-2,-86,30,-44,30,90,-19,-114,
        111,-1400,-230,110,52,-38,49,990,-710,-5;
    w.imagPart() = iw;

    tmv::SmallVector<std::complex<T>,NN> origw = w;
    int perm[NN];

    if (showacc)
        std::cout<<"unsorted w = "<<w<<std::endl;
    w.sort(perm);
    for(int i=1;i<NN;++i) {
        Assert(real(w(i-1)) <= real(w(i)),"Sort complex SmallVector");
    }
    if (showacc)
        std::cout<<"sorted w = "<<w<<std::endl;

    w.reversePermute(perm);
    Assert(w==origw,"Reverse permute sorted SmallVector = orig");
    w.sort(0);
    origw.permute(perm);
    Assert(w==origw,"Permute SmallVector = sorted SmallVector");
}

template <class T> 
static void TestSmallVectorArith()
{
    tmv::SmallVector<T,NN> a;
    for(int i=0;i<NN;++i) a(i) = T(i+10);
    tmv::SmallVector<T,NN> b;
    for(int i=0;i<NN;++i) b(i) = T(-3*i+2);

    tmv::SmallVector<std::complex<T>,NN> ca = a*std::complex<T>(2,-1);;
    tmv::SmallVector<std::complex<T>,NN> cb = b*std::complex<T>(-5,1);

    TestVectorArith1<T>(a,ca,"SmallVector C");
    TestVectorArith2<T>(a,ca,b,cb,"SmallVector CC");

    tmv::SmallVector<T,NN,tmv::FortranStyle> af = a;
    tmv::SmallVector<T,NN,tmv::FortranStyle> bf = b;
    tmv::SmallVector<std::complex<T>,NN,tmv::FortranStyle> caf = ca;
    tmv::SmallVector<std::complex<T>,NN,tmv::FortranStyle> cbf = cb;

    TestVectorArith1<T>(af,caf,"SmallVector F");
    TestVectorArith2<T>(af,caf,bf,cbf,"SmallVector FF");
    TestVectorArith2<T>(a,ca,bf,cbf,"SmallVector CF");
    TestVectorArith2<T>(af,caf,b,cb,"SmallVector FC");

#ifndef NOMIX_SMALL
    tmv::VectorView<T> av = a.view();
    tmv::VectorView<std::complex<T> > cav = ca.view();
    tmv::VectorView<T> bv = b.view();
    tmv::VectorView<std::complex<T> > cbv = cb.view();

    TestVectorArith2<T>(av,cav,b,cb,"SmallVector/Vector");
    TestVectorArith2<T>(a,ca,bv,cbv,"Vector/SmallVector");
#endif
}

template <class T> 
static void TestSmallVectorIO()
{
    tmv::SmallVector<T,NN> v;
    tmv::SmallVector<std::complex<T>,NN> cv;
    for (int i=0; i<NN; ++i) {
        v(i) = T(i+34);
        cv(i) = std::complex<T>(T(i),T(N-i));
    }

    std::ofstream fout("tmvtest_smallvector_io.dat");
    if (!fout) {
#ifdef NOTHROW
        std::cerr<<"Couldn't open tmvtest_smallvector_io.dat for output\n"; 
        exit(1); 
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_smallvector_io.dat for output");
#endif
    }
    fout << v << std::endl << cv << std::endl;
    fout.close();

    tmv::SmallVector<T,NN> xv1;
    tmv::SmallVector<std::complex<T>,NN> xcv1;
    std::ifstream fin("tmvtest_smallvector_io.dat");
    if (!fin) {
#ifdef NOTHROW
        std::cerr<<"Couldn't open tmvtest_smallvector_io.dat for input\n"; 
        exit(1); 
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_smallvector_io.dat for input");
#endif
    }
    fin >> xv1 >> xcv1;
    fin.close();
    Assert(v == xv1,"SmallVector I/O check #1");
    Assert(cv == xcv1,"CSmallVector I/O check #1");

    std::auto_ptr<tmv::Vector<T> > xv2;
    std::auto_ptr<tmv::Vector<std::complex<T> > > xcv2;
    fin.open("tmvtest_smallvector_io.dat");
    if (!fin) {
#ifdef NOTHROW
        std::cerr<<"Couldn't open tmvtest_smallvector_io.dat for input\n"; 
        exit(1); 
#else
        throw std::runtime_error(
            "Couldn't open tmvtest_smallvector_io.dat for input");
#endif
    }
    fin >> xv2 >> xcv2;
    fin.close();
    Assert(v == *xv2,"SmallVector I/O check #2");
    Assert(cv == *xcv2,"CSmallVector I/O check #2");

    std::remove("tmvtest_smallvector_io.dat");
}

template <class T> 
void TestAllSmallVector()
{
    TestSmallVectorReal<T>();
    TestSmallVectorComplex<T>();
    if (tmv::TMV_Epsilon<T>() > T(0)) {
        TestSmallVectorArith<T>();
    }
    TestSmallVectorIO<T>();

    std::cout<<"SmallVector<"<<tmv::TMV_Text(T())<<"> passed all tests\n";
}

#ifdef TEST_DOUBLE
template void TestAllSmallVector<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllSmallVector<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllSmallVector<long double>();
#endif
#ifdef TEST_INT
template void TestAllSmallVector<int>();
#endif
