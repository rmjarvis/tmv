
#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include "TMV_TestVectorArith.h"
#include <fstream>
#include <cstdio>
#include <vector>

#define N 100
#define NN 20

#define CT std::complex<T>

template <class T> 
static void TestSmallVectorReal()
{
    if (showstartdone) {
        std::cout<<"Start Test SmallVector Real"<<std::endl;
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
    }

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
    for (int i=1; i<=N; ++i) 
        Assert(af(i) == a(i-1),"FortranStyle SmallVector access");
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
    tmv::VectorView<T> q3 = q3x.subVector(4,34,5);
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
    tmv::Permutation p;

    if (showacc)
        std::cout<<"unsorted w = "<<w<<std::endl;
    w.sort(p);
    for(int i=1;i<NN;++i) {
        Assert(w(i-1) <= w(i),"Sort real SmallVector");
    }
    if (showacc)
        std::cout<<"sorted w = "<<w<<std::endl;

    w.sort();
    if (showacc)
        std::cout<<"after w.sort(): w = "<<w<<std::endl;
    w = p.inverse() * w;
    if (showacc)
        std::cout<<"Reverse permuted w = "<<w<<std::endl;
    Assert(w==origw,"Reverse permute sorted SmallVector = orig");
    w.sort();
    tmv::SmallVector<T,NN> w2 = p * origw;
    if (showacc)
        std::cout<<"Sort permuted w = "<<w2<<std::endl;
    Assert(w==w2,"Permute SmallVector = sorted SmallVector");
}

template <class T> 
static void TestSmallVectorComplex()
{
    if (showstartdone) {
        std::cout<<"Start Test SmallVector Complex"<<std::endl;
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
    }

    tmv::SmallVector<CT,N> v;
    for (int i=0; i<N; ++i) v(i) = CT(T(i),T(i+1234));

    for (int i=0; i<N; ++i) 
        Assert(real(v(i)) == T(i),"CSmallVector set");
    for (int i=0; i<N; ++i) 
        Assert(imag(v(i)) == T(i+1234),"CSmallVector set");

    tmv::VectorView<CT > v1 = v.subVector(0,N,2);
    for (int i=0; i<N/2; ++i) 
        Assert(v1(i)==CT(T(2*i),T(2*i+1234)),
               "CSmallVector stride=2");

    for (int i=0; i<N/2; ++i) v1[i] = CT(T(i),T(i+9876));
    for (int i=0; i<N/2; ++i) 
        Assert(v[2*i]==CT(T(i),T(i+9876)),
               "setting CSmallVector with stride = 2");

    for (int i=0; i<N; ++i) v(i) = CT(T(i),T(i+1234));

    v.swap(2,5);
    Assert(v[2] == CT(5,5+1234),"Swap in CSmallVector");
    Assert(v[5] == CT(2,2+1234),"Swap in CSmallVector");
    v.swap(2,5);

    tmv::SmallVector<CT,N> v2 = v.conjugate();

    for (int i=0; i<N; ++i) 
        Assert(v2(i) == CT(T(i),T(-i-1234)),
               "Conjugate CSmallVector");
    Assert(v2 == v.conjugate(),"Conjugate == CSmallVector");

    tmv::SmallVector<CT,N> v3;
    for (int i=0; i<N; ++i) v3(i) = CT(i+10,2*i);
    v3(23) = CT(40*N,9*N);
    v3(42) = CT(0,1);
    v3(15) = CT(-32*N,24*N);
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

    CT prod_act(0);
    for (int i=0; i<N; ++i) prod_act += v[i] * v2[i];
    CT prod = v*v2;
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


    CT sum_act(0);
    for (int i=0; i<N; ++i) sum_act += v[i];
    CT sumel = v.sumElements();
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

    tmv::SmallVector<CT,N> ca = a;
    Assert(Equal(ca,a,EPS),"Copy real V -> complex V");

    ca *= CT(3,4);
    tmv::SmallVector<CT,N> cb = b*CT(3,4);
    prod = T(29)*T(25)*CT(-28,96);
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

    tmv::SmallVector<CT,NN> w;
    w <<
        33,12,54,-12,43,-94,0,-20,40,-115,
        -120,140,330,10,-93,-39,49,100,-310,1;
    tmv::SmallVector<T,NN> iw;
    iw <<
        14,98,-2,-86,30,-44,30,90,-19,-114,
        111,-1400,-230,110,52,-38,49,990,-710,-5;
    w.imagPart() = iw;

    tmv::SmallVector<CT,NN> origw = w;
    tmv::Permutation p;

    if (showacc)
        std::cout<<"unsorted w = "<<w<<std::endl;
    w.sort(p);
    for(int i=1;i<NN;++i) {
        Assert(real(w(i-1)) <= real(w(i)),"Sort complex SmallVector");
    }
    if (showacc)
        std::cout<<"sorted w = "<<w<<std::endl;

    w = p.inverse() * w;
    Assert(w==origw,"Reverse permute sorted SmallVector = orig");
    w.sort();
    tmv::SmallVector<CT,NN> w2 = p * origw;
    Assert(w==w2,"Permute SmallVector = sorted SmallVector");
}

template <class T> 
static void TestSmallVectorArith()
{
    if (showstartdone) {
        std::cout<<"Start Test SmallVector Arith"<<std::endl;
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"N = "<<N<<std::endl;
    }

    tmv::SmallVector<T,NN> a;
    for(int i=0;i<NN;++i) a(i) = T(i+10);
    tmv::SmallVector<T,NN> b;
    for(int i=0;i<NN;++i) b(i) = T(-3*i+2);

    tmv::SmallVector<CT,NN> ca = a*CT(2,-1);;
    tmv::SmallVector<CT,NN> cb = b*CT(-5,1);

    TestVectorArith1(a,ca,"SmallVector C");
    TestVectorArith2(a,ca,b,cb,"SmallVector CC");

    tmv::SmallVector<T,NN,tmv::FortranStyle> af = a;
    tmv::SmallVector<T,NN,tmv::FortranStyle> bf = b;
    tmv::SmallVector<CT,NN,tmv::FortranStyle> caf = ca;
    tmv::SmallVector<CT,NN,tmv::FortranStyle> cbf = cb;

    TestVectorArith1(af,caf,"SmallVector F");
    TestVectorArith2(af,caf,bf,cbf,"SmallVector FF");
    TestVectorArith2(a,ca,bf,cbf,"SmallVector CF");
    TestVectorArith2(af,caf,b,cb,"SmallVector FC");

    tmv::VectorView<T> av = a.view();
    tmv::VectorView<CT > cav = ca.view();
    tmv::VectorView<T> bv = b.view();
    tmv::VectorView<CT > cbv = cb.view();

    TestVectorArith2(av,cav,b,cb,"SmallVector/Vector");
    TestVectorArith2(a,ca,bv,cbv,"Vector/SmallVector");
}

template <class T> 
static void TestSmallVectorIO()
{
    if (showstartdone) {
        std::cout<<"Start Test SmallVector I/O"<<std::endl;
        std::cout<<"T = "<<tmv::TMV_Text(T())<<std::endl;
        std::cout<<"NN = "<<NN<<std::endl;
    }

    tmv::SmallVector<T,NN> v;
    tmv::SmallVector<CT,NN> cv;
    for (int i=0; i<NN; ++i) {
        v(i) = T(i+34);
        cv(i) = CT(T(i),T(N-i));
    }
    v(3) = T(1.e-30);
    cv(3) = CT(T(1.e-30),T(1.e-30));
    v(8) = T(9.e-3);
    cv(8) = CT(T(9.e-3),T(9.e-3));
    cv(9) = CT(T(9),T(9.e-3));
    v(12) = T(0.123456789);
    cv(12) = CT(T(3.123456789),T(6.987654321));

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
    Assert(bool(fout),"Couldn't open tmvtest_smallvector_io.dat for output");
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
    cv(12) = CT(T(3.12346),T(6.98765));

    // When using prec(12), the full correct values will be read in. (v2,cv2)

    // When using prec(4), these will be the values read in.
    v3(12) = T(0.1235);
    if (std::numeric_limits<T>::is_integer) cv3(12) = CT(3,6);
    else cv3(12) = CT(T(3.123),T(6.988));

    // Read them back in
    tmv::SmallVector<T,NN> xv;
    tmv::SmallVector<CT,NN> xcv;
    std::ifstream fin("tmvtest_smallvector_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_smallvector_io.dat for input");
    fin >> xv >> xcv;
    Assert(EqualIO(v,xv,EPS),"SmallVector I/O check normal");
    Assert(EqualIO(cv,xcv,EPS),"CSmallVector I/O check normal");
    fin >> tmv::CompactIO() >> xv >> tmv::CompactIO() >> xcv;
    Assert(EqualIO(v,xv,EPS),"SmallVector I/O check compact");
    Assert(EqualIO(cv,xcv,EPS),"CSmallVector I/O check compact");
    fin >> xv.view() >> xcv.view();
    Assert(EqualIO(v2,xv,EPS),"SmallVector I/O check thresh");
    Assert(EqualIO(cv2,xcv,EPS),"CSmallVector I/O check thresh");
    fin >> myStyle >> xv.view() >> myStyle >> xcv.view();
    Assert(EqualIO(v3,xv,EPS),"SmallVector I/O check compact thresh & prec(4)");
    Assert(EqualIO(cv3,xcv,EPS),"CSmallVector I/O check compact thresh & prec(4)");
    fin.close();

    // Read back into regular Vector
    // Also check switching the default IOStyle.
    tmv::CompactIO().makeDefault();
    tmv::Vector<T> zv1, zv2, zv3, zv4;
    tmv::Vector<CT > zcv1, zcv2, zcv3, zcv4;
    fin.open("tmvtest_smallvector_io.dat");
    Assert(bool(fin),"Couldn't open tmvtest_smallvector_io.dat for input");
    fin >> tmv::NormalIO() >> zv1 >> tmv::NormalIO() >> zcv1;
    Assert(EqualIO(v,zv1,EPS),"SmallVector I/O check normal -> Vector");
    Assert(EqualIO(cv,zcv1,EPS),"CSmallVector I/O check normal -> Vector");
    fin >> zv2 >> zcv2;
    Assert(EqualIO(v,zv2,EPS),"SmallVector I/O check compact -> Vector");
    Assert(EqualIO(cv,zcv2,EPS),"CSmallVector I/O check compact -> Vector");
    fin >> tmv::NormalIO() >> zv3 >> tmv::NormalIO() >> zcv3;
    Assert(EqualIO(v2,zv3,EPS),"SmallVector I/O check thresh -> Vector");
    Assert(EqualIO(cv2,zcv3,EPS),"CSmallVector I/O check thresh -> Vector");
    fin >> myStyle >> zv4 >> myStyle >> zcv4;
    Assert(EqualIO(v3,zv4,EPS),"SmallVector I/O check compact thresh -> Vector");
    Assert(EqualIO(cv3,zcv4,EPS),"CSmallVector I/O check compact thresh -> Vector");
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
