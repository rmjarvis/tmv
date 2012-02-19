#define START 0

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymArith.h"

#define NOELEMMULT

#include "TMV_TestMatrixArith.h"

template <class T> 
void TestSymMatrixArith_B2()
{
#if (XTEST & 2)
    std::vector<tmv::SymMatrixView<T> > s;
    std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
    MakeSymList(s,cs,InDef);

    const int N = s[0].size();

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);
    tmv::Matrix<std::complex<T> > ca1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) = 
        std::complex<T>(3+i-5*j,2-3*i);

    tmv::Vector<T> v = a1.col(0);
    tmv::Vector<std::complex<T> > cv = ca1.col(0);

    tmv::MatrixView<T> a1v = a1.view();
    tmv::MatrixView<std::complex<T> > ca1v = ca1.view();

    tmv::Matrix<T> a2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = T(1-3*i+6*j);
    tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) ca2(i,j) = 
        std::complex<T>(1-3*i+6*j,-4+2*j);

    tmv::Vector<T> v5(5*N);
    tmv::VectorView<T> vs = v5.subVector(0,5*N,5);
    vs = v;
    tmv::Vector<std::complex<T> > cv5(5*N);
    tmv::VectorView<std::complex<T> > cvs = cv5.subVector(0,5*N,5);
    cvs = cv;

    tmv::Matrix<T,tmv::RowMajor> a3 = a2.rowRange(0,N);
    tmv::Matrix<std::complex<T> > ca3 = a3 * std::complex<T>(-3,4);
    tmv::Matrix<T,tmv::RowMajor> a4 = a1.colRange(0,0);
    tmv::Matrix<std::complex<T> > ca4 = a4;

    tmv::MatrixView<T> a3v = a3.view();
    tmv::MatrixView<T> a4v = a4.view();
    tmv::MatrixView<std::complex<T> > ca3v = ca3.view();
    tmv::MatrixView<std::complex<T> > ca4v = ca4.view();

    for(size_t i=START;i<s.size();i++) {
        if (showstartdone) {
            std::cout<<"Start loop i = "<<i<<std::endl;
            std::cout<<"si = "<<s[i]<<std::endl;
        }

        tmv::SymMatrixView<T> si = s[i];
        tmv::SymMatrixView<std::complex<T> > csi = cs[i];

        TestMatrixArith4(a1v,ca1v,si,csi,"SquareM/Sym");
        TestMatrixArith5(a1v,ca1v,si,csi,"SquareM/Sym");
        TestMatrixArith6x(a1v,ca1v,si,csi,"SquareM/Sym");
        TestMatrixArith4(a3v,ca3v,si,csi,"NonSquareM/Sym");
        TestMatrixArith5(a3v,ca3v,si,csi,"NonSquareM/Sym");
        TestMatrixArith6x(a3v,ca3v,si,csi,"NonSquareM/Sym");
        TestMatrixArith4(a4v,ca4v,si,csi,"DegenerateM/Sym");
        TestMatrixArith5(a4v,ca4v,si,csi,"DegenerateM/Sym");
        TestMatrixArith6x(a4v,ca4v,si,csi,"DegenerateM/Sym");
    }
#endif
}

#ifdef TEST_DOUBLE
template void TestSymMatrixArith_B2<double>();
#endif
#ifdef TEST_FLOAT
template void TestSymMatrixArith_B2<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymMatrixArith_B2<long double>();
#endif
#ifdef TEST_INT
template void TestSymMatrixArith_B2<int>();
#endif
