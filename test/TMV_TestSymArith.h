
template <class T>
inline void CopyBackM(
    const tmv::Matrix<std::complex<T> >& m0,
    tmv::SymMatrixView<std::complex<T> >& m)
{ 
    if (m.issym()) m = SymMatrixViewOf(m0,m.uplo()); 
    else m = HermMatrixViewOf(m0,m.uplo());
}

#define MakeSymList(s,cs,pdc) \
    std::vector<tmv::SymMatrix<T,tmv::Upper|tmv::RowMajor> > S1; \
    std::vector<tmv::SymMatrix<T,tmv::Upper|tmv::ColMajor> > S2; \
    std::vector<tmv::SymMatrix<T,tmv::Lower|tmv::RowMajor> > S3; \
    std::vector<tmv::SymMatrix<T,tmv::Lower|tmv::ColMajor> > S4; \
    std::vector<tmv::HermMatrix<T,tmv::Upper|tmv::RowMajor> > S5; \
    std::vector<tmv::HermMatrix<T,tmv::Upper|tmv::ColMajor> > S6; \
    std::vector<tmv::HermMatrix<T,tmv::Lower|tmv::RowMajor> > S7; \
    std::vector<tmv::HermMatrix<T,tmv::Lower|tmv::ColMajor> > S8; \
    std::vector<tmv::Matrix<T> > S9; \
    std::vector<tmv::SymMatrix<std::complex<T>,tmv::Upper|tmv::RowMajor> > CS1; \
    std::vector<tmv::SymMatrix<std::complex<T>,tmv::Upper|tmv::ColMajor> > CS2; \
    std::vector<tmv::SymMatrix<std::complex<T>,tmv::Lower|tmv::RowMajor> > CS3; \
    std::vector<tmv::SymMatrix<std::complex<T>,tmv::Lower|tmv::ColMajor> > CS4; \
    std::vector<tmv::HermMatrix<std::complex<T>,tmv::Upper|tmv::RowMajor> > CS5; \
    std::vector<tmv::HermMatrix<std::complex<T>,tmv::Upper|tmv::ColMajor> > CS6; \
    std::vector<tmv::HermMatrix<std::complex<T>,tmv::Lower|tmv::RowMajor> > CS7; \
    std::vector<tmv::HermMatrix<std::complex<T>,tmv::Lower|tmv::ColMajor> > CS8; \
    std::vector<tmv::Matrix<std::complex<T> > > CS9; \
    DoMakeSymList(s,cs,pdc,S1,S2,S3,S4,S5,S6,S7,S8,S9, \
                  CS1,CS2,CS3,CS4,CS5,CS6,CS7,CS8,CS9);

template <class T> 
inline void DoMakeSymList(
    std::vector<tmv::SymMatrixView<T> >& s,
    std::vector<tmv::SymMatrixView<std::complex<T> > >& cs,
    PosDefCode pdc,
    std::vector<tmv::SymMatrix<T,tmv::Upper|tmv::RowMajor> >& SUR,
    std::vector<tmv::SymMatrix<T,tmv::Upper|tmv::ColMajor> >& SUC,
    std::vector<tmv::SymMatrix<T,tmv::Lower|tmv::RowMajor> >& SLR,
    std::vector<tmv::SymMatrix<T,tmv::Lower|tmv::ColMajor> >& SLC,
    std::vector<tmv::HermMatrix<T,tmv::Upper|tmv::RowMajor> >& HUR,
    std::vector<tmv::HermMatrix<T,tmv::Upper|tmv::ColMajor> >& HUC,
    std::vector<tmv::HermMatrix<T,tmv::Lower|tmv::RowMajor> >& HLR,
    std::vector<tmv::HermMatrix<T,tmv::Lower|tmv::ColMajor> >& HLC,
    std::vector<tmv::Matrix<T> >& M,
    std::vector<tmv::SymMatrix<std::complex<T>,tmv::Upper|tmv::RowMajor> >& CSUR,
    std::vector<tmv::SymMatrix<std::complex<T>,tmv::Upper|tmv::ColMajor> >& CSUC,
    std::vector<tmv::SymMatrix<std::complex<T>,tmv::Lower|tmv::RowMajor> >& CSLR,
    std::vector<tmv::SymMatrix<std::complex<T>,tmv::Lower|tmv::ColMajor> >& CSLC,
    std::vector<tmv::HermMatrix<std::complex<T>,tmv::Upper|tmv::RowMajor> >& CHUR,
    std::vector<tmv::HermMatrix<std::complex<T>,tmv::Upper|tmv::ColMajor> >& CHUC,
    std::vector<tmv::HermMatrix<std::complex<T>,tmv::Lower|tmv::RowMajor> >& CHLR,
    std::vector<tmv::HermMatrix<std::complex<T>,tmv::Lower|tmv::ColMajor> >& CHLC,
    std::vector<tmv::Matrix<std::complex<T> > >& CM)
{
    const int N = std::numeric_limits<T>::is_integer ? 6 : 10;
    const int RESERVE = 20;
    SUR.reserve(RESERVE);
    SUC.reserve(RESERVE);
    SLR.reserve(RESERVE);
    SLC.reserve(RESERVE);
    HUR.reserve(RESERVE);
    HUC.reserve(RESERVE);
    HLR.reserve(RESERVE);
    HLC.reserve(RESERVE);
    M.reserve(RESERVE);
    CSUR.reserve(RESERVE);
    CSUC.reserve(RESERVE);
    CSLR.reserve(RESERVE);
    CSLC.reserve(RESERVE);
    CHUR.reserve(RESERVE);
    CHUC.reserve(RESERVE);
    CHLR.reserve(RESERVE);
    CHLC.reserve(RESERVE);
    CM.reserve(RESERVE);
 
    tmv::Matrix<T> a1(N,N);
    tmv::Matrix<std::complex<T> > ca1(N,N);
    tmv::Matrix<T> a2(2*N,2*N);
    tmv::Matrix<std::complex<T> > ca2(2*N,2*N);

    if (std::numeric_limits<T>::is_integer) {
        for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
            a1(i,j) = T(-1+i/2);
            ca1(i,j) = std::complex<T>(-1+i/2,-2+j/3);
        }
        for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) {
            a2(i,j) = T(-3+i);
            ca2(i,j) = std::complex<T>(-3+i,-7+j);
        }
        a1.diag().addToAll(3);
        a2.diag().addToAll(5);
        ca1.diag().addToAll(3);
        ca2.diag().addToAll(5);
    } else {
        for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
            a1(i,j) = T(3+i-5*j);
            ca1(i,j) = std::complex<T>(3+i-5*j,2-3*i);
        }
        for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) {
            a2(i,j) = T(1-3*i+6*j);
            ca2(i,j) = std::complex<T>(1-3*i+6*j,-4+2*j);
        }

        if (pdc == PosDef) {
            a1 /= T(N*N);
            a2 /= T(N*N);
            a1.diag().addToAll(T(1));
            a2.diag().addToAll(T(1));
            for(int i=0;i<N;++i) a1.diag()(i) += T(i);
            for(int i=0;i<2*N;++i) a2.diag()(i) += T(i);
            ca1 /= T(N*N);
            ca2 /= T(N*N);
            ca1.diag().addToAll(T(1));
            ca2.diag().addToAll(T(1));
            for(int i=0;i<N;++i) ca1.diag()(i) += T(i);
            for(int i=0;i<2*N;++i) ca2.diag()(i) += T(i);
        } else if (pdc == InDef) {
            a1 /= T(N);
            a2 /= T(N);
            a1.diag(0,2,4).setZero();
            a2.diag(0,2,4).setZero();
            a1(4,0) = a1(0,4) = T(50);
            a1.diag(-1,6,9).addToAll(T(10));
            a2.diag(-1,6,9).addToAll(T(10));
            a1.diag(1,6,9).addToAll(T(10));
            a2.diag(1,6,9).addToAll(T(10));
            a1.row(9,0,5).addToAll(T(10));
            a2.row(9,0,5).addToAll(T(10));
            a1(3,3) = T(0);
            a2(3,3) = T(0);
            if (N > 10) {
                a1.diag(0,10,N) *= T(0.0001);
                a2.diag(0,10,N) *= T(0.0001);
            }
            ca1 /= T(N);
            ca2 /= T(N);
            ca1.diag(0,2,4).setZero();
            ca2.diag(0,2,4).setZero();
            ca1(4,0) = ca1(0,4) = T(50);
            ca1.diag(-1,6,9).addToAll(T(10));
            ca2.diag(-1,6,9).addToAll(T(10));
            ca1.diag(1,6,9).addToAll(T(10));
            ca2.diag(1,6,9).addToAll(T(10));
            ca1.row(9,0,5).addToAll(T(10));
            ca2.row(9,0,5).addToAll(T(10));
            ca1(3,3) = T(0);
            ca2(3,3) = T(0);
            if (N > 10) {
                ca1.diag(0,10,N) *= T(0.0001);
                ca2.diag(0,10,N) *= T(0.0001);
            }
        } else {
            a1.row(2).setZero();
            a1.col(2).setZero();
            a1.row(4) = a1.row(5);
            a1.col(4) = a1.col(5);
            a1(4,5) = a1(5,4) = a1(5,5) = a1(4,4);
            a2.row(2).setZero();
            a2.col(2).setZero();
            a2.row(4) = a2.row(5);
            a2.col(4) = a2.col(5);
            a2(4,5) = a2(5,4) = a2(5,5) = a2(4,4);
            ca1.row(2).setZero();
            ca1.col(2).setZero();
            ca1.row(4) = ca1.row(5);
            ca1.col(4) = ca1.col(5);
            ca1(4,5) = ca1(5,4) = ca1(5,5) = ca1(4,4);
            ca2.row(2).setZero();
            ca2.col(2).setZero();
            ca2.row(4) = ca2.row(5);
            ca2.col(4) = ca2.col(5);
            ca2(4,5) = ca2(5,4) = ca2(5,5) = ca2(4,4);
        }
    }

    SUR.push_back(tmv::SymMatrix<T,tmv::Upper|tmv::RowMajor>(a1));
    HUR.push_back(tmv::HermMatrix<T,tmv::Upper|tmv::RowMajor>(a1));
    CSUR.push_back(tmv::SymMatrix<std::complex<T>,tmv::Upper|tmv::RowMajor>(ca1));
    CHUR.push_back(tmv::HermMatrix<std::complex<T>,tmv::Upper|tmv::RowMajor>(ca1));
    s.push_back(SUR.back().view());
    s.push_back(HUR.back().view());
    cs.push_back(CSUR.back().view());
    cs.push_back(CHUR.back().view());

    SUC.push_back(tmv::SymMatrix<T,tmv::Upper|tmv::ColMajor>(a1));
    HUC.push_back(tmv::HermMatrix<T,tmv::Upper|tmv::ColMajor>(a1));
    CSUC.push_back(tmv::SymMatrix<std::complex<T>,tmv::Upper|tmv::ColMajor>(ca1));
    CHUC.push_back(tmv::HermMatrix<std::complex<T>,tmv::Upper|tmv::ColMajor>(ca1));
    s.push_back(SUC.back().view());
    s.push_back(HUC.back().view());
    cs.push_back(CSUC.back().view());
    cs.push_back(CHUC.back().view());

#if (XTEST & 2)
    SLR.push_back(tmv::SymMatrix<T,tmv::Lower|tmv::RowMajor>(a1));
    HLR.push_back(tmv::HermMatrix<T,tmv::Lower|tmv::RowMajor>(a1));
    CSLR.push_back(tmv::SymMatrix<std::complex<T>,tmv::Lower|tmv::RowMajor>(ca1));
    CHLR.push_back(tmv::HermMatrix<std::complex<T>,tmv::Lower|tmv::RowMajor>(ca1));
    s.push_back(SLR.back().view());
    s.push_back(HLR.back().view());
    cs.push_back(CSLR.back().view());
    cs.push_back(CHLR.back().view());

    SLC.push_back(tmv::SymMatrix<T,tmv::Lower|tmv::ColMajor>(a1));
    HLC.push_back(tmv::HermMatrix<T,tmv::Lower|tmv::ColMajor>(a1));
    CSLC.push_back(tmv::SymMatrix<std::complex<T>,tmv::Lower|tmv::ColMajor>(ca1));
    CHLC.push_back(tmv::HermMatrix<std::complex<T>,tmv::Lower|tmv::ColMajor>(ca1));
    s.push_back(SLC.back().view());
    s.push_back(HLC.back().view());
    cs.push_back(CSLC.back().view());
    cs.push_back(CHLC.back().view());

    M.push_back(tmv::Matrix<T>(a2));
    CM.push_back(tmv::Matrix<std::complex<T> >(ca2));
    s.push_back(SymMatrixViewOf(M.back(),tmv::Upper).subSymMatrix(0,2*N,2));
    cs.push_back(SymMatrixViewOf(CM.back(),tmv::Upper).subSymMatrix(0,2*N,2));

    M.push_back(tmv::Matrix<T>(a2));
    CM.push_back(tmv::Matrix<std::complex<T> >(ca2));
    CM.back().diag().imagPart().setZero();
    s.push_back(HermMatrixViewOf(M.back(),tmv::Upper).subSymMatrix(0,2*N,2));
    cs.push_back(HermMatrixViewOf(CM.back(),tmv::Upper).subSymMatrix(0,2*N,2));
#endif

    TMVAssert(int(SUR.size()) <= RESERVE);
    TMVAssert(int(SUC.size()) <= RESERVE);
    TMVAssert(int(SLR.size()) <= RESERVE);
    TMVAssert(int(SLC.size()) <= RESERVE);
    TMVAssert(int(HUR.size()) <= RESERVE);
    TMVAssert(int(HUC.size()) <= RESERVE);
    TMVAssert(int(HLR.size()) <= RESERVE);
    TMVAssert(int(HLC.size()) <= RESERVE);
    TMVAssert(int(M.size()) <= RESERVE);
    TMVAssert(int(CSUR.size()) <= RESERVE);
    TMVAssert(int(CSUC.size()) <= RESERVE);
    TMVAssert(int(CSLR.size()) <= RESERVE);
    TMVAssert(int(CSLC.size()) <= RESERVE);
    TMVAssert(int(CHUR.size()) <= RESERVE);
    TMVAssert(int(CHUC.size()) <= RESERVE);
    TMVAssert(int(CHLR.size()) <= RESERVE);
    TMVAssert(int(CHLC.size()) <= RESERVE);
    TMVAssert(int(CM.size()) <= RESERVE);
}


