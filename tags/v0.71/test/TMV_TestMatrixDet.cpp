
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"

template <class T> 
void TestMatrixDet()
{
    tmv::Matrix<T> m0(0,0);

    tmv::Matrix<T> m1(1,1);
    m1 << 
        6;

    tmv::Matrix<T> m2(2,2);
    m2 << 
        6, 1,
        3, 2;

    tmv::Matrix<T> m3(3,3);
    m3 << 
        6, 1, 2,
        3, 2, 1,
        5, 6, 7;

    tmv::Matrix<T> m4(4,4);
    m4 << 
        6, 1, 2, 4,
        3, 2, 1, 8,
        5, 6, 7, 1,
        2, 9, 8, 3;

    tmv::Matrix<T> m5(5,5);
    m5 << 
        6, 1, 2, 4, 1,
        3, 2, 1, 8, 4,
        5, 6, 7, 1, 2,
        2, 9, 8, 3, 7,
        1, 4, 5, 9, 3;
    
    tmv::Matrix<T> m6(5,5);
    m6 << 
        6, 5, 2, 4, 2,
        3, 0, 9, 3, 6,
        5, 5, 7, 1, 8,
        2, 10, 8, 3, 4,
        1, 15, 5, 9, 4;

    tmv::Matrix<T> m7(5,5);
    m7 << 
        6, 1, 0, 4, 1,
        3, 0, 0, 8, 0,
        5, 6, 7, 1, 2,
        2, 9, 0, 3, 0,
        0, 0, 0, 9, 0;

    tmv::Matrix<T> m8(5,5);
    m8 << 
        6, 1, 2, 4, 1,
        3, 2, 8, 8, 4,
        5, 6, 4, 1, 2,
        2, 9, 8, 3, 4,
        1, 4, 6, 9, 3;

    tmv::Matrix<T> m9 = m5.row(0) ^ m5.row(0);

    tmv::Matrix<T> m10 = m5.rowRange(0,4).transpose() * m5.rowRange(0,4);

    tmv::Matrix<T> m11(10,10);
    m11 << 
        0, 2, 3, 0, 0, 2, 6, 0, 1, 0,
        0, 0, 2, 6, 0, 1, 0, 2, 1, 0,
        1, 4, 0, 0, 3, 0, 3, 0, 0, 3,
        7, 0, 0, 0, 8, 0, 4, 0, 0, 5,
        0, 0, 0, 3, 0, 1, 0, 3, 1, 0,
        0, 2, 1, 0, 4, 0, 9, 0, 0, 6,
        0, 4, 7, 6, 0, 1, 0, 0, 4, 0,
        0, 0, 0, 1, 0, 5, 0, 1, 0, 0, 
        4, 0, 0, 0, 6, 0, 0, 0, 0, 4,
        0, 1, 8, 0, 0, 0, 3, 0, 9, 2;

    tmv::Matrix<T> m12(10,10);
    m12 << 
        9, 2, 3, 5, 3, 2, 6, 1, 1, 8,
        5, 6, 2, 6, 3, 1, 2, 2, 1, 8,
        1, 4, 4, 7, 3, 3, 3, 1, 9, 3,
        7, 8, 4, 5, 8, 8, 4, 6, 3, 5,
        2, 1, 2, 3, 8, 1, 4, 3, 1, 6,
        4, 2, 1, 3, 4, 3, 9, 1, 4, 6,
        7, 4, 7, 6, 8, 1, 3, 8, 4, 4,
        8, 8, 3, 1, 1, 5, 8, 1, 3, 4,
        4, 1, 4, 2, 6, 1, 5, 8, 3, 4,
        3, 1, 8, 8, 3, 2, 3, 8, 9, 2;

    tmv::Matrix<std::complex<T> > cm1(5,5);
    cm1.realPart() = m5;
    cm1.imagPart() = m7.transpose();

    tmv::Matrix<std::complex<T> > cm2(10,10);
    for(int i=0;i<10;++i) for(int j=0;j<10;++j) 
        if (m11(i,j) == T(0)) 
            cm2(i,j) = T(0);
        else
            cm2(i,j) = std::complex<T>(m11(i,j)-5, 3-m11(i,j));

    tmv::Matrix<std::complex<T> > cm3(10,10);
    cm3.realPart() = m12;  cm3.realPart().addToAll(-4);
    cm3.imagPart() = T(2)*m12.transpose();
    for(int i=0;i<10;++i) for(int j=0;j<10;++j) {
        // The intent here is to get the determinant as close to MAX_INT
        // for a 32 bit integer as possible (in either component) without
        // going over.  i.e. close to 2e9.
        // This is to prove that the algorithm works for any matrix that
        // has its answer expressible as an int (or complex<int> in this case).
        if (cm3.realPart()(i,j) >= T(5)) cm3.realPart()(i,j) = T(0);
        if (cm3.imagPart()(i,j) > T(5)) cm3.imagPart()(i,j) -= T(10);
    }

    if (showacc) {
        std::cout<<"m0: "<<m0.det()<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS2(m0.det())<<"  "<<EPS<<std::endl;
        std::cout<<"m1: "<<m1.det()<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS2(m1.det()-6)<<"  "<<
            6*NormSq(m1)*EPS<<std::endl;
        std::cout<<"m2: "<<m2.det()<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS2(m2.det()-9)<<"  "<<
            9*NormSq(m2)*EPS<<std::endl;
        std::cout<<"m3: "<<m3.det()<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS2(m3.det()-48)<<"  "<<
            48*NormSq(m3)*EPS<<std::endl;
        std::cout<<"m4: "<<m4.det()<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS2(m4.det()+66)<<"  "<<
            66*NormSq(m4)*EPS<<std::endl;
        std::cout<<"m5: "<<m5.det()<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS2(m5.det()+1118)<<"  "<<
            1000*NormSq(m5)*EPS<<std::endl;
        std::cout<<"m6: "<<m6.det()<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS2(m6.det()+15360)<<"  "<<
            15000*NormSq(m6)*EPS<<std::endl;
        std::cout<<"m7: "<<m7.det()<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS2(m7.det()+1701)<<"  "<<
            1700*NormSq(m7)*EPS<<std::endl;
        std::cout<<"m8: "<<m8.det()<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS2(m8.det())<<"  "<<
            1000*NormSq(m8)*EPS<<std::endl;
        std::cout<<"m9: "<<m9.det()<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS2(m9.det())<<"  "<<
            1000*NormSq(m9)*EPS<<std::endl;
        std::cout<<"m10: "<<m10.det()<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS2(m10.det())<<"  "<<
            1000*NormSq(m10)*EPS<<std::endl;
        std::cout<<"m11: "<<m11.det()<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS2(m11.det()-251656)<<"  "<<
            250000*NormSq(m11)*EPS<<std::endl;
        std::cout<<"m12: "<<m12.det()<<std::endl;
        std::cout<<"diff = "<<tmv::TMV_ABS2(m12.det()-20762581)<<"  "<<
            20000000*NormSq(m12)*EPS<<std::endl;
        std::cout<<"cm1: "<<cm1.det()<<std::endl;
        std::cout<<"diff = "<<
            tmv::TMV_ABS2(cm1.det()-std::complex<T>(-19796,-2165))<<"  "<<
            20000*NormSq(cm1)*EPS<<std::endl;
        std::cout<<"cm2: "<<cm2.det()<<std::endl;
        std::cout<<"diff = "<<
            tmv::TMV_ABS2(cm2.det()-std::complex<T>(2386176,-1084032))<<"  "<<
            3000000*NormSq(cm2)*EPS<<std::endl;
        std::cout<<"cm3: "<<cm3.det()<<std::endl;
        std::cout<<"diff = "<<
            tmv::TMV_ABS2(cm3.det()-std::complex<T>(-1923395548,811584412))<<"  "<<
            2000000000*NormSq(cm3)*EPS<<std::endl;
    }
    Assert(Equal2(m0.det(),1,EPS),"0x0 determinant");
    Assert(Equal2(m1.det(),6,6*NormSq(m1)*EPS),"1x1 determinant");
    Assert(Equal2(m2.det(),9,9*NormSq(m2)*EPS),"2x2 determinant");
    Assert(Equal2(m3.det(),48,48*NormSq(m3)*EPS),"3x3 determinant");
    Assert(Equal2(m4.det(),-66,66*NormSq(m4)*EPS),"4x4 determinant");
    Assert(Equal2(m5.det(),-1118,1000*NormSq(m5)*EPS),"5x5 determinant");
    Assert(Equal2(m6.det(),-15360,15000*NormSq(m6)*EPS),
           "5x5 determinant with GCD factors");
    Assert(Equal2(m7.det(),-1701,1700*NormSq(m7)*EPS),
           "5x5 determinant permuted triangle");
    Assert(Equal2(m8.det(),0,1000*NormSq(m8)*EPS),"5x5 determinant singular");
    Assert(Equal2(m9.det(),0,1000*NormSq(m9)*EPS),
           "5x5 determinant rank 1 outer product");
    Assert(Equal2(m10.det(),0,1000*NormSq(m10)*EPS),
           "5x5 determinant rank 4 outer product");
    Assert(Equal2(m11.det(),251656,250000*NormSq(m11)*EPS),
           "10x10 determinant permuted band");
    Assert(Equal2(m12.det(),20762581,20000000*NormSq(m11)*EPS),
           "10x10 determinant full");
    Assert(Equal2(cm1.det(),std::complex<T>(-19796,-2165),
                  20000*NormSq(cm1)*EPS),
           "5x5 determinant complex");
    Assert(Equal2(cm2.det(),std::complex<T>(2386176,-1084032),
                  4000000*NormSq(cm2)*EPS),
           "10x10 determinant complex permuted band");
    Assert(Equal2(cm3.det(),std::complex<T>(-1923395548,811584412),
                  2000000000*NormSq(cm3)*EPS),
           "10x10 determinant complex full");

    std::cout<<"Matrix<"<<tmv::TMV_Text(T())<<"> passed all determinant tests\n";
}

#ifdef TEST_DOUBLE
template void TestMatrixDet<double>();
#endif
#ifdef TEST_FLOAT
template void TestMatrixDet<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestMatrixDet<long double>();
#endif
#ifdef TEST_INT
template void TestMatrixDet<int>();
#endif
