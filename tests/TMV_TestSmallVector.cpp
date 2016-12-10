
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include "TMV.h"
#include <fstream>
#include <cstdio>

template <class T> 
static void TestSmallVectorReal()
{
    typedef typename tmv::Traits<T>::float_type FT;

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

    // Test assignments and constructors from arrays
    std::vector<T> qv(6);
    T qvar[] = { T(8), T(6), T(4), T(2), T(0), T(-2) };
    for(int i=0;i<6;++i) qv[i] = qvar[i];

    // Construct from C array
    tmv::SmallVector<T,6> q1;
    std::copy(qvar,qvar+6,q1.begin());
    // Construct from vector
    tmv::SmallVector<T,6> q2;
    std::copy(qv.begin(),qv.end(),q2.begin());
    // Assign SmallVectorView from vector
    tmv::SmallVector<T,50> q3x;
    typename tmv::SmallVector<T,50>::subvector_step_type q3 = 
        q3x.subVector(4,34,5);
    std::copy(qv.begin(),qv.end(),q3.begin());
    // Use op<< assignment
    tmv::SmallVector<T,6> q4;
    q4 << 8, 6, 4, 2, 0, -2;

    if (showacc) {
        std::cout<<"q1 = "<<q1<<std::endl;
        std::cout<<"q2 = "<<q2<<std::endl;
        std::cout<<"q3 = "<<q3<<std::endl;
        std::cout<<"q4 = "<<q4<<std::endl;
    }

    for(int i=0;i<6;++i) {
        Assert(q1(i) == T(8-2*i),"Create Vector from T*");
        Assert(q2(i) == T(8-2*i),"Create Vector from std::vector");
        Assert(q3(i) == T(8-2*i),"Create VectorView from std::vector");
        Assert(q4(i) == T(8-2*i),"Create Vector from << list");
    }

    // Test Basic Arithmetic
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

    T prod = 0;
    FT normsum = 0, normdiff = 0;
    for(int i=0;i<N;++i) {
        prod += a[i] * b[i];
        normsum += (a[i]+b[i])*(a[i]+b[i]);
        normdiff += (a[i]-b[i])*(a[i]-b[i]);
    }
    normsum = tmv::TMV_SQRT(normsum);
    normdiff = tmv::TMV_SQRT(normdiff);
    Assert(Equal2(a*b , prod, EPS*Norm(a)*Norm(b)),"Inner Product");
    Assert(Equal2(Norm(a+b) , normsum, EPS*(Norm1(a)+Norm1(b))),
           "SmallVector Sum");
    Assert(Equal2(Norm(a-b) , normdiff, EPS*(Norm1(a)+Norm1(b))),
           "SmallVector Diff");

    tmv::SmallVector<T,20> w;
    w << 3.3,1.2,5.4,-1.2,4.3,-9.4,0,-2,4,-11.5,
      -12,14,33,1,-9.3,-3.9,4.9,10,-31,1.e-33;

    tmv::SmallVector<T,20> origw = w;
    tmv::Permutation P(20);

    if (showacc)
        std::cout<<"unsorted w = "<<w<<std::endl;
    w.sort(P);
    for(int i=1;i<20;++i) {
        Assert(w(i-1) <= w(i),"Sort real SmallVector");
    }
    if (showacc)
        std::cout<<"sorted w = "<<w<<std::endl;

    w.sort();
    P.inverse().applyOnLeft(w);
    if (showacc)
        std::cout<<"Reverse permuted w = "<<w<<std::endl;
    Assert(w==origw,"Reverse permute sorted SmallVector = orig");
    w.sort();
    P.applyOnLeft(origw);
    if (showacc)
        std::cout<<"Sort permuted w = "<<origw<<std::endl;
    Assert(w==origw,"Permute SmallVector = sorted SmallVector");

}

template <class T> 
static void TestSmallVectorComplex()
{
    typedef std::complex<T> CT;
    typedef typename tmv::Traits<T>::float_type FT;

    const int N = 100;

    tmv::SmallVector<CT,N> v;
    for (int i=0; i<N; ++i) v(i) = CT(T(i),T(i+1234));

    for (int i=0; i<N; ++i) 
        Assert(v(i).real() == T(i), "CSmallVector set");
    for (int i=0; i<N; ++i) 
        Assert(v(i).imag() == T(i+1234), "CSmallVector set");

    if (N % 2 == 0) {
        tmv::SmallVectorView<CT,N/2,2> v1 = v.subVector(0,N,2);
        for (int i=0; i<N/2; ++i) 
            Assert(v1(i)==CT(T(2*i),T(2*i+1234)),
                   "CSmallVector stride=2");

        for (int i=0; i<N/2; ++i) v1[i] = CT(T(i),T(i+9876));
        for (int i=0; i<N/2; ++i) 
            Assert(v[2*i]==CT(T(i),T(i+9876)),
                   "setting CSmallVector with stride = 2");

        for (int i=0; i<N; ++i) v(i) = CT(T(i),T(i+1234));
    }

    v.swap(2,5);
    Assert(v[2] == CT(T(5),T(5+1234)),"Swap in CSmallVector");
    Assert(v[5] == CT(T(2),T(2+1234)),"Swap in CSmallVector");
    v.swap(2,5);

    tmv::SmallVector<CT,N> v2 = v.conjugate();

    for (int i=0; i<N; ++i) 
        Assert(v2(i) == CT(T(i),T(-i-1234)),
               "Conjugate CSmallVector");
    Assert(v2 == v.conjugate(),"Conjugate == CSmallVector");

    Assert(tmv::TMV_ABS((v*v2).imag()) <= EPS,"CSmallVector * CSmallVector");
    FT norm1 = tmv::TMV_SQRT((v*v2).real());
    FT norm2 = Norm(v);
    if (showacc) {
        std::cout<<"v = "<<v<<std::endl;
        std::cout<<"v2 = "<<v2<<std::endl;
        std::cout<<"v*v2 = "<<v*v2<<std::endl;
        std::cout<<"norm1 = "<<norm1<<std::endl;
        std::cout<<"norm2 = "<<norm2<<std::endl;
    }
    Assert(tmv::TMV_ABS(norm1 - norm2) <= EPS*norm1,"Norm CSmallVector");

    Assert(v2 == v.conjugateSelf(),"ConjugateSelf CSmallVector");

    tmv::SmallVector<T,N> a;
    for(int i=0;i<N;++i) a(i) = T(i+10);
    tmv::SmallVector<T,N> b;
    for(int i=0;i<N;++i) b(i) = T(-3*i+191);

    tmv::SmallVector<CT,N> ca = a;
    Assert(Equal(ca,a,EPS*Norm(a)),"Copy real V -> complex V");
    ca *= CT(T(3),T(4));
    tmv::SmallVector<CT,N> cb = b*CT(T(3),T(4));

    CT prod = 0;
    FT normsum = 0, normdiff = 0;
    for(int i=0;i<N;++i) {
        prod += ca[i] * cb[i];
        normsum += std::norm(ca[i]+cb[i]);
        normdiff += std::norm(ca[i]-cb[i]);
    }
    normsum = tmv::TMV_SQRT(normsum);
    normdiff = tmv::TMV_SQRT(normdiff);
    Assert(Equal2(ca*cb , prod, EPS*Norm(ca)*Norm(cb)),"CInner Product");
    Assert(Equal2(Norm(ca+cb) , normsum, EPS*(Norm(ca)+Norm(cb))),
           "CSmallVector Sum");
    Assert(Equal2(Norm(ca-cb) , normdiff, EPS*(Norm(ca)+Norm(cb))),
           "CSmallVector Diff");

    tmv::SmallVector<CT,20> w;
    w << 3.3,1.2,5.4,-1.2,4.3,-9.4,0,-2,4,-11.5,
      -12,14,33,1,-9.3,-3.9,4.9,10,-31,1.e-33;
    tmv::SmallVector<T,20> iw;
    iw << 1.4,9.8,-0.2,-8.6,3.0,-4.4,3,9,-1.9,-11.4,
       11.1,-140,-23,11,5.2,-3.8,4.9,99,-71,-0.5;
    w.imagPart() = iw;

    tmv::SmallVector<CT,20> origw = w;
    tmv::Permutation P(20);

    if (showacc)
        std::cout<<"unsorted w = "<<w<<std::endl;
    w.sort(P);
    for(int i=1;i<20;++i) {
        Assert(w(i-1).real() <= w(i).real(),"Sort complex SmallVector");
    }
    if (showacc)
        std::cout<<"sorted w = "<<w<<std::endl;

    P.inverse().applyOnLeft(w);
    Assert(w==origw,"Reverse permute sorted SmallVector = orig");
    w.sort();
    P.applyOnLeft(origw);
    Assert(w==origw,"Permute SmallVector = sorted SmallVector");
}

template <class T> 
static void TestSmallVectorIO()
{
    typedef std::complex<T> CT;

    const int NN = 20;

    if (showstartdone) {
        std::cout<<"Start Test SmallVector I/O"<<std::endl;
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"NN = "<<NN<<std::endl;
    }

    tmv::SmallVector<T,NN> v;
    tmv::SmallVector<CT,NN> cv;
    for (int i=0; i<NN; ++i) {
        v(i) = T(i+34);
        cv(i) = CT(T(i),T(100-i));
    }
    v(3) = T(1.e-30);
    cv(3) = CT(T(1.e-30),T(1.e-30));
    v(8) = T(9.e-3);
    cv(8) = CT(T(9.e-3),T(9.e-3));
    cv(9) = CT(T(9),T(9.e-3));
    v(12) = T(0.123456789);
    cv(12) = CT(T(3.123456789),T(600.987654321));

    // First check clipping function...
    tmv::SmallVector<T,NN> v2 = v;
    tmv::SmallVector<CT,NN> cv2 = cv;
    if (!std::numeric_limits<T>::is_integer) {
        v2.clip(T(1.e-2));
        cv2.clip(T(1.e-2));
    }
    tmv::SmallVector<T,NN> v3 = v;
    tmv::SmallVector<CT,NN> cv3 = cv;
    v3(3) = T(0);
    cv3(3) = T(0);
    v3(8) = T(0); // Others, esp. cv3(8), shouldn't get clipped.
    Assert(v2 == v3,"SmallVector clip");
    Assert(cv2 == cv3,"Complex SmallVector clip");

    // However, ThreshIO for complex works slightly differently than clip.
    // It clips _either_ the real or imag component, so now cv2(8) and 
    // cv2(9) need to be modified.
    cv2(8) = cv3(8) = T(0);
    cv2(9) = cv3(9) = T(9);

    // Write vectors with 4 different styles
    std::ofstream fout("tmvtest_smallvector_io.dat");
    Assert(fout,"Couldn't open tmvtest_smallvector_io.dat for output");
    fout << v << std::endl;
    fout << cv << std::endl;
    fout << tmv::CompactIO() << v << std::endl;
    fout << tmv::CompactIO() << cv << std::endl;
    fout << tmv::ThreshIO(1.e-2).setPrecision(12) << v << std::endl;
    fout << tmv::ThreshIO(1.e-2).setPrecision(12) << cv << std::endl;
    // Not a very pretty IO style, but it tests being able to read
    // a style that has no whitespace and has more than one 
    // character for some of the markup elements.
    tmv::IOStyle myStyle = 
        tmv::CompactIO().setThresh(1.e-2).setPrecision(4).
        markup("Start","[",",","]","---","Done");
    fout << myStyle << v << std::endl;
    fout << myStyle << cv << std::endl;
    fout.close();

    // When using (the default) prec(6), these will be the values read in.
    v(12) = T(0.123457);
    cv(12) = CT(T(3.12346),T(600.988));

    // When using prec(12), the full correct values will be read in. (v2,cv2)

    // When using prec(4), these will be the values read in.
    v3(12) = T(0.1235);
    if (std::numeric_limits<T>::is_integer) cv3(12) = CT(T(3),T(600));
    else cv3(12) = CT(T(3.123),T(601.0));

    // Read them back in
    tmv::SmallVector<T,NN> xv;
    tmv::SmallVector<CT,NN> xcv;
    std::ifstream fin("tmvtest_smallvector_io.dat");
    Assert(fin,"Couldn't open tmvtest_smallvector_io.dat for input");
    fin >> xv >> xcv;
    Assert(v == xv,"SmallVector I/O check normal");
    Assert(cv == xcv,"CSmallVector I/O check normal");
    fin >> tmv::CompactIO() >> xv >> tmv::CompactIO() >> xcv;
    Assert(v == xv,"SmallVector I/O check compact");
    Assert(cv == xcv,"CSmallVector I/O check compact");
    fin >> xv.view() >> xcv.view();
    Assert(v2 == xv,"SmallVector I/O check thresh");
    Assert(cv2 == xcv,"CSmallVector I/O check thresh");
    fin >> myStyle >> xv.view() >> myStyle >> xcv.view();
    Assert(v3 == xv,"SmallVector I/O check compact thresh & prec(4)");
    Assert(cv3 == xcv,"CSmallVector I/O check compact thresh & prec(4)");
    fin.close();

    // Read back into regular Vector
    // Also check switching the default IOStyle.
    tmv::CompactIO().makeDefault();
    tmv::Vector<T> zv1, zv2, zv3, zv4;
    tmv::Vector<CT > zcv1, zcv2, zcv3, zcv4;
    fin.open("tmvtest_smallvector_io.dat");
    Assert(fin,"Couldn't open tmvtest_smallvector_io.dat for input");
    fin >> tmv::NormalIO() >> zv1 >> tmv::NormalIO() >> zcv1;
    Assert(v == zv1,"SmallVector I/O check normal -> Vector");
    Assert(cv == zcv1,"CSmallVector I/O check normal -> Vector");
    fin >> zv2 >> zcv2;
    Assert(v == zv2,"SmallVector I/O check compact -> Vector");
    Assert(cv == zcv2,"CSmallVector I/O check compact -> Vector");
    fin >> tmv::NormalIO() >> zv3 >> tmv::NormalIO() >> zcv3;
    Assert(v2 == zv3,"SmallVector I/O check thresh -> Vector");
    Assert(cv2 == zcv3,"CSmallVector I/O check thresh -> Vector");
    fin >> myStyle >> zv4 >> myStyle >> zcv4;
    Assert(v3 == zv4,"SmallVector I/O check compact thresh -> Vector");
    Assert(cv3 == zcv4,"CSmallVector I/O check compact thresh -> Vector");
    fin.close();
    // Switch it back.
    tmv::IOStyle::revertDefault();

#if XTEST == 0
    std::remove("tmvtest_smallvector_io.dat");
#endif
}

template <class T> 
void TestAllSmallVector()
{
#if 1
    TestSmallVectorReal<T>();
    TestSmallVectorComplex<T>();
    TestSmallVectorIO<T>();
    std::cout<<"SmallVector<"<<Text(T())<<"> passed all basic tests\n";
#endif

#if 1
    TestSmallVectorArith_1a<T>();
    TestSmallVectorArith_1b<T>();
    TestSmallVectorArith_2a<T>();
    TestSmallVectorArith_2b<T>();
    TestSmallVectorArith_2c<T>();
    TestSmallVectorArith_2d<T>();
    std::cout<<"SmallVector<"<<Text(T())<<"> passed all arithmetic tests\n";
#endif
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
