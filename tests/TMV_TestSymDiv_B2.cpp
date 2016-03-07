
#define START 0

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymArith.h"

#define NOLDIVEQ
#define NORDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> 
void TestSymDiv_B2(tmv::DivType dt, PosDefCode pdc)
{
    const int N = 10;

    std::vector<tmv::SymMatrixView<T> > s;
    std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
    MakeSymList(s,cs,pdc);

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(1-3*i+j);
    a1.diag().addToAll(T(10)*N);
    a1 /= T(10);
    tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);
    a1.diag().addToAll(T(10)*N);
    ca1.diag().addToAll(T(10)*N);

    tmv::MatrixView<T> a1v = a1.view();
    tmv::MatrixView<std::complex<T> > ca1v = ca1.view();

#if (XTEST & 2)
    tmv::Matrix<T> a3 = a1.colRange(0,N/2);
    tmv::Matrix<std::complex<T> > ca3 = ca1.colRange(0,N/2);
    tmv::Matrix<T> a4 = a1.rowRange(0,N/2);
    tmv::Matrix<std::complex<T> > ca4 = ca1.rowRange(0,N/2);
    tmv::Matrix<T> a5(2*N,N);
    a5.rowRange(0,N) = a1;
    a5.rowRange(N,2*N) = a1;
    tmv::Matrix<std::complex<T> > ca5(2*N,N);
    ca5.rowRange(0,N) = ca1;
    ca5.rowRange(N,2*N) = ca1;
    tmv::Matrix<T> a6 = a5.transpose();
    tmv::Matrix<std::complex<T> > ca6 = ca5.transpose();

    tmv::MatrixView<T> a3v = a3.view();
    tmv::MatrixView<T> a4v = a4.view();
    tmv::MatrixView<T> a5v = a5.view();
    tmv::MatrixView<T> a6v = a6.view();
    tmv::MatrixView<std::complex<T> > ca3v = ca3.view();
    tmv::MatrixView<std::complex<T> > ca4v = ca4.view();
    tmv::MatrixView<std::complex<T> > ca5v = ca5.view();
    tmv::MatrixView<std::complex<T> > ca6v = ca6.view();
#endif

    for(size_t i=START;i<s.size();i++) {
        if (showstartdone)
            std::cout<<"Start loop: i = "<<i<<", si = "<<
                tmv::TMV_Text(s[i])<<"  "<<s[i]<<std::endl;
        tmv::SymMatrixView<T> si = s[i];
        tmv::SymMatrixView<std::complex<T> > csi = cs[i];

        si.saveDiv();
        csi.saveDiv();

        TestMatrixDivArith1(dt,a1v,si,ca1v,csi,"Sym/SquareMatrix");
        if (dt == tmv::LU) continue;
#if (XTEST & 2)
        TestMatrixDivArith1(dt,a3v,si,ca3v,csi,"Sym/NonSquareMatrix");
        TestMatrixDivArith1(dt,a4v,si,ca4v,csi,"Sym/NonSquareMatrix");
        TestMatrixDivArith1(dt,a5v,si,ca5v,csi,"Sym/NonSquareMatrix");
        TestMatrixDivArith1(dt,a6v,si,ca6v,csi,"Sym/NonSquareMatrix");
#endif
    }
}

#ifdef TEST_DOUBLE
template void TestSymDiv_B2<double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef TEST_FLOAT
template void TestSymDiv_B2<float>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef TEST_LONGDOUBLE
template void TestSymDiv_B2<long double>(tmv::DivType dt, PosDefCode pdc);
#endif
