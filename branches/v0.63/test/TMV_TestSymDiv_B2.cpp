
#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_Sym.h"
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
    std::vector<tmv::BaseMatrix<T>*> B;
    std::vector<tmv::BaseMatrix<std::complex<T> >*> CB;
    MakeSymList(s,cs,B,CB,pdc);

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(1-3*i+j);
    a1.diag().AddToAll(T(10)*N);
    a1 /= T(10);
    tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);
    a1.diag().AddToAll(T(10)*N);
    ca1.diag().AddToAll(T(10)*N);

    tmv::MatrixView<T> a1v = a1.View();
    tmv::MatrixView<std::complex<T> > ca1v = ca1.View();

#ifdef XTEST
    tmv::Matrix<T> a3 = a1.Cols(0,N/2);
    tmv::Matrix<std::complex<T> > ca3 = ca1.Cols(0,N/2);
    tmv::Matrix<T> a4 = a1.Rows(0,N/2);
    tmv::Matrix<std::complex<T> > ca4 = ca1.Rows(0,N/2);
    tmv::Matrix<T> a5(2*N,N);
    a5.Rows(0,N) = a1;
    a5.Rows(N,2*N) = a1;
    tmv::Matrix<std::complex<T> > ca5(2*N,N);
    ca5.Rows(0,N) = ca1;
    ca5.Rows(N,2*N) = ca1;
    tmv::Matrix<T> a6 = a5.Transpose();
    tmv::Matrix<std::complex<T> > ca6 = ca5.Transpose();

    tmv::MatrixView<T> a3v = a3.View();
    tmv::MatrixView<T> a4v = a4.View();
    tmv::MatrixView<T> a5v = a5.View();
    tmv::MatrixView<T> a6v = a6.View();
    tmv::MatrixView<std::complex<T> > ca3v = ca3.View();
    tmv::MatrixView<std::complex<T> > ca4v = ca4.View();
    tmv::MatrixView<std::complex<T> > ca5v = ca5.View();
    tmv::MatrixView<std::complex<T> > ca6v = ca6.View();
#endif

    for(size_t i=START;i<s.size();i++) {
        if (showstartdone)
            std::cout<<"Start loop: i = "<<i<<", si = "<<
                tmv::TMV_Text(s[i])<<"  "<<s[i]<<std::endl;
        const tmv::SymMatrixView<T>& si = s[i];
        const tmv::SymMatrixView<std::complex<T> >& csi = cs[i];

        si.SaveDiv();
        csi.SaveDiv();

        if (csi.issym()) {
            tmv::SymMatrix<T> sx = si;
            tmv::SymMatrix<std::complex<T> > csx = csi;

            TestMatrixDivArith1<T>(dt,sx,csx,a1v,si,ca1v,csi,
                                   "Sym/SquareMatrix");
            if (dt == tmv::LU) continue;
#ifdef XTEST
            TestMatrixDivArith1<T>(dt,sx,csx,a3v,si,ca3v,csi,
                                   "Sym/NonSquareMatrix");
            TestMatrixDivArith1<T>(dt,sx,csx,a4v,si,ca4v,csi,
                                   "Sym/NonSquareMatrix");
            TestMatrixDivArith1<T>(dt,sx,csx,a5v,si,ca5v,csi,
                                   "Sym/NonSquareMatrix");
            TestMatrixDivArith1<T>(dt,sx,csx,a6v,si,ca6v,csi,
                                   "Sym/NonSquareMatrix");
#endif
        } else {
            tmv::HermMatrix<T> hx = si;
            tmv::HermMatrix<std::complex<T> > chx = csi;

            TestMatrixDivArith1<T>(dt,hx,chx,a1v,si,ca1v,csi,
                                   "Herm/SquareMatrix");
            if (dt == tmv::LU) continue;
#ifdef XTEST
            TestMatrixDivArith1<T>(dt,hx,chx,a3v,si,ca3v,csi,
                                   "Herm/NonSquareMatrix");
            TestMatrixDivArith1<T>(dt,hx,chx,a4v,si,ca4v,csi,
                                   "Herm/NonSquareMatrix");
            TestMatrixDivArith1<T>(dt,hx,chx,a5v,si,ca5v,csi,
                                   "Herm/NonSquareMatrix");
            TestMatrixDivArith1<T>(dt,hx,chx,a6v,si,ca6v,csi,
                                   "Herm/NonSquareMatrix");
#endif
        }
    }
    for(size_t i=0;i<B.size();++i) delete B[i];
    for(size_t i=0;i<CB.size();++i) delete CB[i];
}

#ifdef INST_DOUBLE
template void TestSymDiv_B2<double>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_FLOAT
template void TestSymDiv_B2<float>(tmv::DivType dt, PosDefCode pdc);
#endif
#ifdef INST_LONGDOUBLE
template void TestSymDiv_B2<long double>(tmv::DivType dt, PosDefCode pdc);
#endif
