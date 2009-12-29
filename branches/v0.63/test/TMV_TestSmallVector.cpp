
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_SmallVectorArith.h"
#include <fstream>
#include <cstdio>

#include "TMV_TestVectorArith.h"

#define N 100
#define NN 20

template <class T> 
static void TestSmallVectorReal()
{
    tmv::SmallVector<T,N> v;

    for (int i=0; i<N; ++i) v(i) = T(i);

    for (int i=0; i<N; ++i) Assert(v(i) == T(i),"Setting SmallVector");

    tmv::VectorView<T> v1 = v.SubVector(0,N,2);
    for (int i=0; i<N/2; ++i) Assert(v1(i) == T(2*i), "SmallVector stride=2");

    for (int i=0; i<N/2; ++i) v1[i] = T(i+1234);
    for (int i=0; i<N/2; ++i) Assert(v[2*i] == T(i+1234),
                                     "setting SmallVector with stride = 2");
    for (int i=0; i<N; ++i) v(i) = T(i);

    v.Swap(2,5);
    Assert(v(2) == T(5) && v(5) == T(2),"Swapping elements of SmallVector");
    v.Swap(2,5);
    Assert(v(2) == T(2) && v(5) == T(5),"Swapping elements of SmallVector");

    T sum = N*(N-1)/2;
    Assert(SumElements(v) == sum,"SmallVector SumElements(v)");

    v.ReverseSelf();
    for (int i=0; i<N; ++i) Assert(v(i) == T(N-i-1),"Reversing SmallVector");

    for (int i=0; i<N; ++i) v(i) = T(i+10);
    v(23) = T(10*N);
    v(42) = T(1)/T(4);
    v(15) = T(-20*N);
    int imax,imin;
    Assert(v.MaxAbsElement(&imax) == T(20*N),
           "MaxAbsElement of SmallVector did not return correct value");
    Assert(imax == 15,
           "MaxAbsElement of SmallVector did not return correct index");
    Assert(v.MinAbsElement(&imin) == T(1)/T(4),
           "MinAbsElement of SmallVector did not return correct value");
    Assert(imin == 42,
           "MinAbsElement of SmallVector did not return correct index");
    Assert(v.MaxElement(&imax) == T(10*N),
           "MaxElement of SmallVector did not return correct value");
    Assert(imax == 23,
           "MaxElement of SmallVector did not return correct index");
    Assert(v.MinElement(&imin) == T(-20*N),
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

    tmv::SmallVector<T,5> c = v.SubVector(10,70,12);
    for (int i=0; i<5; ++i) Assert(c(i) == v(10+12*i),"SubSmallVector");

    for(int i=0;i<N;++i) a(i) = T(i+10);
    for(int i=0;i<N;++i) b(i) = T(-3*i+191);

    T prod = 2900;
    T normsum = tmv::TMV_SQRT(T(1373700));
    T normdiff = tmv::TMV_SQRT(T(1362100));
    Assert(std::abs(a*b - prod) <= EPS*Norm(a)*Norm(b),"Inner Product");
    tmv::SmallVector<T,N> temp;
    Assert(std::abs(Norm(temp=a+b) - normsum) <= EPS*std::abs(Norm1(a)+Norm1(b)),"SmallVector Sum");
    Assert(std::abs(Norm(temp=a-b) - normdiff) <= EPS*std::abs(Norm1(a)+Norm1(b)),"SmallVector Diff");

    tmv::SmallVector<T,NN> w;
    w <<
        3.3,1.2,5.4,-1.2,4.3,-9.4,0,-2,4,-11.5,
        -12,14,33,1,-9.3,-3.9,4.9,10,-31,1.e-33;

    tmv::SmallVector<T,NN> origw = w;
    int perm[NN];

    if (showacc)
        std::cout<<"unsorted w = "<<w<<std::endl;
    w.Sort(perm);
    for(int i=1;i<NN;++i) {
        Assert(w(i-1) <= w(i),"Sort real SmallVector");
    }
    if (showacc)
        std::cout<<"sorted w = "<<w<<std::endl;

    w.Sort(0);
    w.ReversePermute(perm);
    if (showacc)
        std::cout<<"Reverse permuted w = "<<w<<std::endl;
    Assert(w==origw,"Reverse permute sorted SmallVector = orig");
    w.Sort(0);
    origw.Permute(perm);
    if (showacc)
        std::cout<<"Sort permuted w = "<<origw<<std::endl;
    Assert(w==origw,"Permute SmallVector = sorted SmallVector");
}

template <class T> 
static void TestSmallVectorComplex()
{
    tmv::SmallVector<std::complex<T>,N> v;
    for (int i=0; i<N; ++i) v(i) = std::complex<T>(T(i),T(i+1234));

    for (int i=0; i<N; ++i) Assert(v(i).real() == T(i),
                                   "CSmallVector set");
    for (int i=0; i<N; ++i) Assert(v(i).imag() == T(i+1234),
                                   "CSmallVector set");

    tmv::VectorView<std::complex<T> > v1 = v.SubVector(0,N,2);
    for (int i=0; i<N/2; ++i) Assert(v1(i)==std::complex<T>(T(2*i),T(2*i+1234)),
                                     "CSmallVector stride=2");

    for (int i=0; i<N/2; ++i) v1[i] = std::complex<T>(T(i),T(i+9876));
    for (int i=0; i<N/2; ++i) Assert(v[2*i]==std::complex<T>(T(i),T(i+9876)),
                                     "setting CSmallVector with stride = 2");

    for (int i=0; i<N; ++i) v(i) = std::complex<T>(T(i),T(i+1234));

    v.Swap(2,5);
    Assert(v[2] == std::complex<T>(5,5+1234),"Swap in CSmallVector");
    Assert(v[5] == std::complex<T>(2,2+1234),"Swap in CSmallVector");
    v.Swap(2,5);

    tmv::SmallVector<std::complex<T>,N> v2 = v.Conjugate();

    for (int i=0; i<N; ++i) Assert(v2(i) == std::complex<T>(T(i),T(-i-1234)),
                                   "Conjugate CSmallVector");
    Assert(v2 == v.Conjugate(),"Conjugate == CSmallVector");

    Assert(std::abs((v*v2).imag()) <= EPS,"CSmallVector * CSmallVector");
    T norm1 = tmv::TMV_SQRT((v*v2).real());
    T norm2 = Norm(v);
    if (showacc) {
        std::cout<<"v = "<<v<<std::endl;
        std::cout<<"v2 = "<<v2<<std::endl;
        std::cout<<"v*v2 = "<<v*v2<<std::endl;
        std::cout<<"norm1 = "<<norm1<<std::endl;
        std::cout<<"norm2 = "<<norm2<<std::endl;
    }
    Assert(std::abs(norm1 - norm2) <= EPS*norm1,"Norm CSmallVector");

    Assert(v2 == v.ConjugateSelf(),"ConjugateSelf CSmallVector");

    tmv::SmallVector<T,N> a;
    for(int i=0;i<N;++i) a(i) = T(i+10);
    tmv::SmallVector<T,N> b;
    for(int i=0;i<N;++i) b(i) = T(-3*i+191);

    tmv::SmallVector<std::complex<T>,N> ca = a;
    tmv::SmallVector<std::complex<T>,N> temp;
    Assert(Norm(temp=ca-a) <= EPS*Norm(a),"Copy real V -> complex V");
    ca *= std::complex<T>(3,4);
    tmv::SmallVector<std::complex<T>,N> cb = b*std::complex<T>(3,4);

    std::complex<T> prod = T(29)*std::complex<T>(-28,96)*T(25);
    T normsum = tmv::TMV_SQRT(T(1373700)*T(25));
    T normdiff = tmv::TMV_SQRT(T(1362100)*T(25));
    Assert(std::abs(ca*cb - prod) <= EPS*Norm(ca)*Norm(cb),"CInner Product");
    Assert(std::abs(Norm(temp=ca+cb) - normsum) <= EPS*std::abs(Norm(ca)+Norm(cb)),"CSmallVector Sum");
    Assert(std::abs(Norm(temp=ca-cb) - normdiff) <= EPS*std::abs(Norm(ca)+Norm(cb)),"CSmallVector Diff");

    tmv::SmallVector<std::complex<T>,NN> w;
    w <<
        3.3,1.2,5.4,-1.2,4.3,-9.4,0,-2,4,-11.5,
        -12,14,33,1,-9.3,-3.9,4.9,10,-31,1.e-33;
    tmv::SmallVector<T,NN> iw;
    iw <<
        1.4,9.8,-0.2,-8.6,3.0,-4.4,3,9,-1.9,-11.4,
        11.1,-140,-23,11,5.2,-3.8,4.9,99,-71,-0.5;
    w.Imag() = iw;

    tmv::SmallVector<std::complex<T>,NN> origw = w;
    int perm[NN];

    if (showacc)
        std::cout<<"unsorted w = "<<w<<std::endl;
    w.Sort(perm);
    for(int i=1;i<NN;++i) {
        Assert(w(i-1).real() <= w(i).real(),"Sort complex SmallVector");
    }
    if (showacc)
        std::cout<<"sorted w = "<<w<<std::endl;

    //w.Sort(0);
    w.ReversePermute(perm);
    Assert(w==origw,"Reverse permute sorted SmallVector = orig");
    w.Sort(0);
    origw.Permute(perm);
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

    tmv::SmallVector<T,NN> a0;
    tmv::SmallVector<std::complex<T>,NN> ca0;

    TestVectorArith1<T>(a0,ca0,a,ca,"SmallVector C");
    TestVectorArith2<T>(a0,ca0,a,ca,b,cb,"SmallVector CC");

    tmv::SmallVector<T,NN,tmv::FortranStyle> af = a;
    tmv::SmallVector<T,NN,tmv::FortranStyle> bf = b;
    tmv::SmallVector<std::complex<T>,NN,tmv::FortranStyle> caf = ca;
    tmv::SmallVector<std::complex<T>,NN,tmv::FortranStyle> cbf = cb;

    TestVectorArith1<T>(a0,ca0,af,caf,"SmallVector F");
    TestVectorArith2<T>(a0,ca0,af,caf,bf,cbf,"SmallVector FF");
    TestVectorArith2<T>(a0,ca0,a,ca,bf,cbf,"SmallVector CF");
    TestVectorArith2<T>(a0,ca0,af,caf,b,cb,"SmallVector FC");

#ifdef XTEST
    // These tests lead to segmentation faults with ATLAS BLAS.
    tmv::VectorView<T> av = a.View();
    tmv::VectorView<std::complex<T> > cav = ca.View();
    tmv::VectorView<T> bv = b.View();
    tmv::VectorView<std::complex<T> > cbv = cb.View();

    TestVectorArith2<T>(a0,ca0,av,cav,b,cb,"SmallVector/Vector");
    TestVectorArith2<T>(a0,ca0,a,ca,bv,cbv,"Vector/SmallVector");
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

#ifndef XTEST
    std::remove("tmvtest_smallvector_io.dat");
#endif
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

#ifdef INST_DOUBLE
template void TestAllSmallVector<double>();
#endif
#ifdef INST_FLOAT
template void TestAllSmallVector<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestAllSmallVector<long double>();
#endif
#ifdef INST_INT
template void TestAllSmallVector<int>();
#endif
