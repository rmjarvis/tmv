
#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV_Vec.h"
#include <fstream>
#include <cstdio>

template <class T> static void TestSmallVectorReal()
{
    const int N = 100;

    tmv::SmallVector<T,N> v;

    for (int i=0; i<N; ++i) v(i) = T(i);

    for (int i=0; i<N; ++i) Assert(v(i) == T(i),"Setting SmallVector");

    if (N % 2 == 0) {
        tmv::SmallVectorView<T,N/2,2> v1 = v.subVector(0,N,2);
        for (int i=0; i<N/2; ++i) Assert(v1(i) == T(2*i), "SmallVector stride=2");

        for (int i=0; i<N/2; ++i) v1[i] = T(i+1234);
        for (int i=0; i<N/2; ++i) Assert(v[2*i] == T(i+1234),
                                         "setting SmallVector with stride = 2");
        for (int i=0; i<N; ++i) v(i) = T(i);
    }

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
           "maxAbsElement of SmallVector did not return correct value");
    Assert(imax == 15,
           "maxAbsElement of SmallVector did not return correct index");
    Assert(v.minAbsElement(&imin) == T(1)/T(4),
           "minAbsElement of SmallVector did not return correct value");
    Assert(imin == 42,
           "minAbsElement of SmallVector did not return correct index");
    Assert(v.maxElement(&imax) == T(10*N),
           "maxElement of SmallVector did not return correct value");
    Assert(imax == 23,
           "maxElement of SmallVector did not return correct index");
    Assert(v.minElement(&imin) == T(-20*N),
           "minElement of SmallVector did not return correct value");
    Assert(imin == 15,
           "minElement of SmallVector did not return correct index");

    tmv::SmallVector<T,N> a;
    tmv::SmallVector<T,N> b;
    for (int i=0; i<N; ++i) a(i) = T(3+i);

    b = a;
    for (int i=0; i<N; ++i) Assert(a(i) == b(i),"SmallVector1 = SmallVector2");

    Assert(a == b,"Testing Equality of SmallVectors");

    b(4) = 0;
    Assert(a != b,"Testing Inequality of SmallVectors");

    tmv::SmallVector<T,N,tmv::FortranStyle> af;
    for (int i=1; i<=N; ++i) af(i) = T(3+i-1);
    for (int i=1; i<=N; ++i) 
        Assert(af(i) == a(i-1), "FortranStyle SmallVector access");
    Assert(a == af,"FortransStyle SmallVector = CStyle SmallVector");

    for (int i=0; i<N; ++i) b(i) = T(5+2*i);

    v = a+b;
    for (int i=0; i<N; ++i) 
        Assert(v(i) == T(8+3*i),"Adding SmallVectors");

    v = a-b;
    for (int i=0; i<N; ++i) 
        Assert(v(i) == T(-2-i),"Subtracting SmallVectors");

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

    T prod = 0, normsum = 0, normdiff = 0;
    for(int i=0;i<N;++i) {
        prod += a[i] * b[i];
        normsum += (a[i]+b[i])*(a[i]+b[i]);
        normdiff += (a[i]-b[i])*(a[i]-b[i]);
    }
    normsum = tmv::TMV_SQRT(normsum);
    normdiff = tmv::TMV_SQRT(normdiff);
    Assert(std::abs(a*b - prod) <= EPS*Norm(a)*Norm(b),"Inner Product");
    tmv::SmallVector<T,N> temp;
    Assert(std::abs(Norm(temp=a+b) - normsum) <= 
           EPS*std::abs(Norm1(a)+Norm1(b)),"SmallVector Sum");
    Assert(std::abs(Norm(temp=a-b) - normdiff) <= 
           EPS*std::abs(Norm1(a)+Norm1(b)),"SmallVector Diff");

    tmv::SmallVector<T,20> w;
    w << 3.3,1.2,5.4,-1.2,4.3,-9.4,0,-2,4,-11.5,
      -12,14,33,1,-9.3,-3.9,4.9,10,-31,1.e-33;

    tmv::SmallVector<T,20> origw = w;
    int perm[20];

    if (showacc)
        std::cout<<"unsorted w = "<<w<<std::endl;
    w.sort(perm);
    for(int i=1;i<20;++i) {
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

template <class T> static void TestSmallVectorComplex()
{
    const int N = 100;

    tmv::SmallVector<std::complex<T>,N> v;
    for (int i=0; i<N; ++i) v(i) = std::complex<T>(T(i),T(i+1234));

    for (int i=0; i<N; ++i) 
        Assert(v(i).real() == T(i), "CSmallVector set");
    for (int i=0; i<N; ++i) 
        Assert(v(i).imag() == T(i+1234), "CSmallVector set");

    if (N % 2 == 0) {
        tmv::SmallVectorView<std::complex<T>,N/2,2> v1 = v.subVector(0,N,2);
        for (int i=0; i<N/2; ++i) 
            Assert(v1(i)==std::complex<T>(T(2*i),T(2*i+1234)),
                   "CSmallVector stride=2");

        for (int i=0; i<N/2; ++i) v1[i] = std::complex<T>(T(i),T(i+9876));
        for (int i=0; i<N/2; ++i) 
            Assert(v[2*i]==std::complex<T>(T(i),T(i+9876)),
                   "setting CSmallVector with stride = 2");

        for (int i=0; i<N; ++i) v(i) = std::complex<T>(T(i),T(i+1234));
    }

    v.swap(2,5);
    Assert(v[2] == std::complex<T>(5,5+1234),"Swap in CSmallVector");
    Assert(v[5] == std::complex<T>(2,2+1234),"Swap in CSmallVector");
    v.swap(2,5);

    tmv::SmallVector<std::complex<T>,N> v2 = v.conjugate();

    for (int i=0; i<N; ++i) 
        Assert(v2(i) == std::complex<T>(T(i),T(-i-1234)),
               "Conjugate CSmallVector");
    Assert(v2 == v.conjugate(),"Conjugate == CSmallVector");

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

    Assert(v2 == v.conjugateSelf(),"ConjugateSelf CSmallVector");

    tmv::SmallVector<T,N> a;
    for(int i=0;i<N;++i) a(i) = T(i+10);
    tmv::SmallVector<T,N> b;
    for(int i=0;i<N;++i) b(i) = T(-3*i+191);

    tmv::SmallVector<std::complex<T>,N> ca = a;
    tmv::SmallVector<std::complex<T>,N> temp;
    Assert(Norm(temp=ca-a) <= EPS*Norm(a),"Copy real V -> complex V");
    ca *= std::complex<T>(3,4);
    tmv::SmallVector<std::complex<T>,N> cb = b*std::complex<T>(3,4);

    std::complex<T> prod = 0;
    T normsum = 0, normdiff = 0;
    for(int i=0;i<N;++i) {
        prod += ca[i] * cb[i];
        normsum += std::norm(ca[i]+cb[i]);
        normdiff += std::norm(ca[i]-cb[i]);
    }
    normsum = tmv::TMV_SQRT(normsum);
    normdiff = tmv::TMV_SQRT(normdiff);
    Assert(std::abs(ca*cb - prod) <= EPS*Norm(ca)*Norm(cb),"CInner Product");
    Assert(std::abs(Norm(temp=ca+cb) - normsum) <= 
           EPS*std::abs(Norm(ca)+Norm(cb)),"CSmallVector Sum");
    Assert(std::abs(Norm(temp=ca-cb) - normdiff) <= 
           EPS*std::abs(Norm(ca)+Norm(cb)),"CSmallVector Diff");

    tmv::SmallVector<std::complex<T>,20> w;
    w << 3.3,1.2,5.4,-1.2,4.3,-9.4,0,-2,4,-11.5,
      -12,14,33,1,-9.3,-3.9,4.9,10,-31,1.e-33;
    tmv::SmallVector<T,20> iw;
    iw << 1.4,9.8,-0.2,-8.6,3.0,-4.4,3,9,-1.9,-11.4,
       11.1,-140,-23,11,5.2,-3.8,4.9,99,-71,-0.5;
    w.imagPart() = iw;

    tmv::SmallVector<std::complex<T>,20> origw = w;
    int perm[20];

    if (showacc)
        std::cout<<"unsorted w = "<<w<<std::endl;
    w.sort(perm);
    for(int i=1;i<20;++i) {
        Assert(w(i-1).real() <= w(i).real(),"Sort complex SmallVector");
    }
    if (showacc)
        std::cout<<"sorted w = "<<w<<std::endl;

    //w.sort(0);
    w.reversePermute(perm);
    Assert(w==origw,"Reverse permute sorted SmallVector = orig");
    w.sort(0);
    origw.permute(perm);
    Assert(w==origw,"Permute SmallVector = sorted SmallVector");
}

template <class T> static void TestSmallVectorIO()
{
    const int N = 100;

    tmv::SmallVector<T,N> v;
    tmv::SmallVector<std::complex<T>,N> cv;
    for (int i=0; i<N; ++i) {
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

    tmv::SmallVector<T,N> xv1;
    tmv::SmallVector<std::complex<T>,N> xcv1;
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

template <class T> void TestSmallVector()
{
#if 1
    TestSmallVectorReal<T>();
    TestSmallVectorComplex<T>();
    TestSmallVectorIO<T>();
    std::cout<<"SmallVector<"<<tmv::TMV_Text(T())<<"> passed all basic tests\n";
#endif

#if 1
    TestSmallVectorArith_1a<T>();
    TestSmallVectorArith_1b<T>();
    TestSmallVectorArith_2a<T>();
    TestSmallVectorArith_2b<T>();
    TestSmallVectorArith_2c<T>();
    TestSmallVectorArith_2d<T>();
    std::cout<<"SmallVector<"<<tmv::TMV_Text(T())<<"> passed all arithmetic tests\n";
#endif
}

#ifdef TEST_DOUBLE
template void TestSmallVector<double>();
#endif
#ifdef TEST_FLOAT
template void TestSmallVector<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSmallVector<long double>();
#endif
#ifdef TEST_INT
template void TestSmallVector<int>();
#endif
