
template <class T>
inline void CopyBackM(
    const tmv::BandMatrix<std::complex<T> >& m0,
    tmv::SymBandMatrixView<std::complex<T> >& m)
{ 
    if (m.issym()) m = SymBandMatrixViewOf(m0,m.uplo()); 
    else m = HermBandMatrixViewOf(m0,m.uplo());
}

#define MakeSymBandList(s,cs,pdc) \
    std::vector<tmv::SymBandMatrix<T,tmv::Upper|tmv::RowMajor> > sB1; \
    std::vector<tmv::SymBandMatrix<T,tmv::Upper|tmv::ColMajor> > sB2; \
    std::vector<tmv::SymBandMatrix<T,tmv::Upper|tmv::DiagMajor> > sB3; \
    std::vector<tmv::SymBandMatrix<T,tmv::Lower|tmv::RowMajor> > sB4; \
    std::vector<tmv::SymBandMatrix<T,tmv::Lower|tmv::ColMajor> > sB5; \
    std::vector<tmv::SymBandMatrix<T,tmv::Lower|tmv::DiagMajor> > sB6; \
    std::vector<tmv::HermBandMatrix<T,tmv::Upper|tmv::RowMajor> > sB7; \
    std::vector<tmv::HermBandMatrix<T,tmv::Upper|tmv::ColMajor> > sB8; \
    std::vector<tmv::HermBandMatrix<T,tmv::Upper|tmv::DiagMajor> > sB9; \
    std::vector<tmv::HermBandMatrix<T,tmv::Lower|tmv::RowMajor> > sB10; \
    std::vector<tmv::HermBandMatrix<T,tmv::Lower|tmv::ColMajor> > sB11; \
    std::vector<tmv::HermBandMatrix<T,tmv::Lower|tmv::DiagMajor> > sB12; \
    std::vector<tmv::Matrix<T> > sB13; \
    std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Upper|tmv::RowMajor> > CsB1; \
    std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Upper|tmv::ColMajor> > CsB2; \
    std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Upper|tmv::DiagMajor> > CsB3; \
    std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Lower|tmv::RowMajor> > CsB4; \
    std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Lower|tmv::ColMajor> > CsB5; \
    std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Lower|tmv::DiagMajor> > CsB6; \
    std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Upper|tmv::RowMajor> > CsB7; \
    std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Upper|tmv::ColMajor> > CsB8; \
    std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Upper|tmv::DiagMajor> > CsB9; \
    std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Lower|tmv::RowMajor> > CsB10; \
    std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Lower|tmv::ColMajor> > CsB11; \
    std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Lower|tmv::DiagMajor> > CsB12; \
    std::vector<tmv::Matrix<std::complex<T> > > CsB13; \
    DoMakeSymBandList( \
        s,cs,pdc,sB1,sB2,sB3,sB4,sB5,sB6,sB7,sB8,sB9,sB10,sB11,sB12,sB13, \
        CsB1,CsB2,CsB3,CsB4,CsB5,CsB6,CsB7,CsB8,CsB9,CsB10,CsB11,CsB12,CsB13);

template <class T> inline void DoMakeSymBandList(
    std::vector<tmv::SymBandMatrixView<T> >& s,
    std::vector<tmv::SymBandMatrixView<std::complex<T> > >& cs,
    PosDefCode pdc,
    std::vector<tmv::SymBandMatrix<T,tmv::Upper|tmv::RowMajor> >& SUR,
    std::vector<tmv::SymBandMatrix<T,tmv::Upper|tmv::ColMajor> >& SUC,
    std::vector<tmv::SymBandMatrix<T,tmv::Upper|tmv::DiagMajor> >& SUD,
    std::vector<tmv::SymBandMatrix<T,tmv::Lower|tmv::RowMajor> >& SLR,
    std::vector<tmv::SymBandMatrix<T,tmv::Lower|tmv::ColMajor> >& SLC,
    std::vector<tmv::SymBandMatrix<T,tmv::Lower|tmv::DiagMajor> >& SLD,
    std::vector<tmv::HermBandMatrix<T,tmv::Upper|tmv::RowMajor> >& HUR,
    std::vector<tmv::HermBandMatrix<T,tmv::Upper|tmv::ColMajor> >& HUC,
    std::vector<tmv::HermBandMatrix<T,tmv::Upper|tmv::DiagMajor> >& HUD,
    std::vector<tmv::HermBandMatrix<T,tmv::Lower|tmv::RowMajor> >& HLR,
    std::vector<tmv::HermBandMatrix<T,tmv::Lower|tmv::ColMajor> >& HLC,
    std::vector<tmv::HermBandMatrix<T,tmv::Lower|tmv::DiagMajor> >& HLD,
    std::vector<tmv::Matrix<T> >& M,
    std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Upper|tmv::RowMajor> >& CSUR,
    std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Upper|tmv::ColMajor> >& CSUC,
    std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Upper|tmv::DiagMajor> >& CSUD,
    std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Lower|tmv::RowMajor> >& CSLR,
    std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Lower|tmv::ColMajor> >& CSLC,
    std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Lower|tmv::DiagMajor> >& CSLD,
    std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Upper|tmv::RowMajor> >& CHUR,
    std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Upper|tmv::ColMajor> >& CHUC,
    std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Upper|tmv::DiagMajor> >& CHUD,
    std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Lower|tmv::RowMajor> >& CHLR,
    std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Lower|tmv::ColMajor> >& CHLC,
    std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Lower|tmv::DiagMajor> >& CHLD,
    std::vector<tmv::Matrix<std::complex<T> > >& CM)
{
    // The integer det function has trouble with the 10x10 determinant,
    // since it overflows.  6x6 is ok.
    const int N=std::numeric_limits<T>::is_integer ? 6 : 10;
    const int RESERVE = 20;
    SUR.reserve(RESERVE);
    SUC.reserve(RESERVE);
    SUD.reserve(RESERVE);
    SLR.reserve(RESERVE);
    SLC.reserve(RESERVE);
    SLD.reserve(RESERVE);
    HUR.reserve(RESERVE);
    HUC.reserve(RESERVE);
    HUD.reserve(RESERVE);
    HLR.reserve(RESERVE);
    HLC.reserve(RESERVE);
    HLD.reserve(RESERVE);
    M.reserve(RESERVE);
    CSUR.reserve(RESERVE);
    CSUC.reserve(RESERVE);
    CSUD.reserve(RESERVE);
    CSLR.reserve(RESERVE);
    CSLC.reserve(RESERVE);
    CSLD.reserve(RESERVE);
    CHUR.reserve(RESERVE);
    CHUC.reserve(RESERVE);
    CHUD.reserve(RESERVE);
    CHLR.reserve(RESERVE);
    CHLC.reserve(RESERVE);
    CHLD.reserve(RESERVE);
    CM.reserve(RESERVE);
 
 
    tmv::Matrix<T> a1(N,N);
    tmv::Matrix<T> a2(2*N,2*N);
    tmv::Matrix<std::complex<T> > ca1(N,N);
    tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
    tmv::Vector<T> v1(N);
    tmv::Vector<T> v2(N-1);
    tmv::Vector<std::complex<T> > cv1(N);
    tmv::Vector<std::complex<T> > cv2(N-1);
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
        for (int i=0; i<N; ++i) v1(i) = T(1+i);
        for (int i=0; i<N-1; ++i) v2(i) = T(-3+2*i);
        for (int i=0; i<N; ++i) cv1(i) = std::complex<T>(1+i,4+i);
        for (int i=0; i<N-1; ++i) cv2(i) = std::complex<T>(3-2*i,6-i);
    } else if (pdc == PosDef) {
        for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
            a1(i,j) = T(3+i+5*j);
            ca1(i,j) = std::complex<T>(3+i+5*j,2+3*i);
        }
        for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) {
            a2(i,j) = T(1+3*i+6*j);
            ca2(i,j) = std::complex<T>(1+3*i+6*j,4+2*j);
        }
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
        for (int i=0; i<N; ++i) v1(i) = T(16+3*i);
        for (int i=0; i<N-1; ++i) v2(i) = T(6+i);
        for (int i=0; i<N; ++i) cv1(i) = std::complex<T>(16+3*i,i+4);
        for (int i=0; i<N-1; ++i) cv2(i) = std::complex<T>(2*i+3,+6+i);
    } else if (pdc == InDef) {
        for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
            a1(i,j) = T(3+i+5*j);
            ca1(i,j) = std::complex<T>(3+i+5*j,2+3*i);
        }
        for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) {
            a2(i,j) = T(1+3*i+6*j);
            ca2(i,j) = std::complex<T>(1+3*i+6*j,4+2*j);
        }
        a1 /= T(N*N);
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
        ca1(4,0) = a1(0,4) = T(50);
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
        for (int i=0; i<N; ++i) v1(i) = T(16-3*i);
        for (int i=0; i<N-1; ++i) v2(i) = T(-6+i);
        for (int i=0; i<N; ++i) cv1(i) = std::complex<T>(16-3*i,i+4);
        for (int i=0; i<N-1; ++i) cv2(i) = std::complex<T>(2*i-3,-6+i);
    } else {
        for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
            a1(i,j) = T(8+i-5*j);
            ca1(i,j) = std::complex<T>(8+i-5*j,9-3*i);
        }
        for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) {
            a2(i,j) = T(-6-3*i+6*j);
            ca2(i,j) = std::complex<T>(-6-3*i+6*j,-4+2*j);
        }
        a1 /= T(N*N);
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
        for (int i=0; i<N; ++i) v1(i) = T(18-3*i);
        for (int i=0; i<N-1; ++i) v2(i) = T(-6+i);
        for (int i=0; i<N; ++i) cv1(i) = std::complex<T>(18-3*i,i-6);
        for (int i=0; i<N-1; ++i) cv2(i) = std::complex<T>(2*i-3,-6+i);
        v1.subVector(N/2,3*N/4).setZero();
        v2.subVector(N/2,3*N/4).setZero();
        cv1.subVector(N/2,3*N/4).setZero();
        cv2.subVector(N/2,3*N/4).setZero();
    }
    tmv::Vector<std::complex<T> > cv1r = cv1;
    cv1r.imagPart().setZero();

    SUR.push_back(tmv::SymBandMatrix<T,tmv::Upper|tmv::RowMajor>(a1,3));
    HUR.push_back(tmv::HermBandMatrix<T,tmv::Upper|tmv::RowMajor>(a1,3));
    CSUR.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper|tmv::RowMajor>(ca1,3));
    CHUR.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper|tmv::RowMajor>(ca1,3));
    s.push_back(SUR.back().view());
    s.push_back(HUR.back().view());
    cs.push_back(CSUR.back().view());
    cs.push_back(CHUR.back().view());

    SUC.push_back(tmv::SymBandMatrix<T,tmv::Upper|tmv::ColMajor>(a1,3));
    HUC.push_back(tmv::HermBandMatrix<T,tmv::Upper|tmv::ColMajor>(a1,3));
    CSUC.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper|tmv::ColMajor>(ca1,3));
    CHUC.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper|tmv::ColMajor>(ca1,3));
    s.push_back(SUC.back().view());
    s.push_back(HUC.back().view());
    cs.push_back(CSUC.back().view());
    cs.push_back(CHUC.back().view());

    SUD.push_back(tmv::SymBandMatrix<T,tmv::Upper|tmv::DiagMajor>(a1,3));
    HUD.push_back(tmv::HermBandMatrix<T,tmv::Upper|tmv::DiagMajor>(a1,3));
    CSUD.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper|tmv::DiagMajor>(ca1,3));
    CHUD.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper|tmv::DiagMajor>(ca1,3));
    s.push_back(SUD.back().view());
    s.push_back(HUD.back().view());
    cs.push_back(CSUD.back().view());
    cs.push_back(CHUD.back().view());

    SUD.push_back(tmv::SymTriDiagMatrix(v1,v2));
    HUD.push_back(tmv::HermTriDiagMatrix(v1,v2,tmv::Upper));
    CSUD.push_back(tmv::SymTriDiagMatrix(cv1,cv2));
    CHUD.push_back(tmv::HermTriDiagMatrix(v1,cv2,tmv::Upper));
    s.push_back(SUD.back().view());
    s.push_back(HUD.back().view());
    cs.push_back(CSUD.back().view());
    cs.push_back(CHUD.back().view());

    SUD.push_back(tmv::SymBandMatrix<T,tmv::Upper|tmv::DiagMajor>(
            DiagMatrixViewOf(v1)));
    HUD.push_back(tmv::HermBandMatrix<T,tmv::Upper|tmv::DiagMajor>(
            DiagMatrixViewOf(v1)));
    CSUD.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper|tmv::DiagMajor>(
            DiagMatrixViewOf(cv1)));
    CHUD.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper|tmv::DiagMajor>(
            DiagMatrixViewOf(cv1)));
    s.push_back(SUD.back().view());
    s.push_back(HUD.back().view());
    cs.push_back(CSUD.back().view());
    cs.push_back(CHUD.back().view());

#if (XTEST & 2)
    SLC.push_back(tmv::SymBandMatrix<T,tmv::Lower|tmv::ColMajor>(a1,N-2));
    HLC.push_back(tmv::HermBandMatrix<T,tmv::Lower|tmv::ColMajor>(a1,N-2));
    CSLC.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower|tmv::ColMajor>(ca1,N-2));
    CHLC.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower|tmv::ColMajor>(ca1,N-2));
    s.push_back(SLC.back().view());
    s.push_back(HLC.back().view());
    cs.push_back(CSLC.back().view());
    cs.push_back(CHLC.back().view());

    SLR.push_back(tmv::SymBandMatrix<T,tmv::Lower|tmv::RowMajor>(a1,N-2));
    HLR.push_back(tmv::HermBandMatrix<T,tmv::Lower|tmv::RowMajor>(a1,N-2));
    CSLR.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower|tmv::RowMajor>(ca1,N-2));
    CHLR.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower|tmv::RowMajor>(ca1,N-2));
    s.push_back(SLR.back().view());
    s.push_back(HLR.back().view());
    cs.push_back(CSLR.back().view());
    cs.push_back(CHLR.back().view());

    SLD.push_back(tmv::SymBandMatrix<T,tmv::Lower|tmv::DiagMajor>(a1,N-2));
    HLD.push_back(tmv::HermBandMatrix<T,tmv::Lower|tmv::DiagMajor>(a1,N-2));
    CSLD.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower|tmv::DiagMajor>(ca1,N-2));
    CHLD.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower|tmv::DiagMajor>(ca1,N-2));
    s.push_back(SLD.back().view());
    s.push_back(HLD.back().view());
    cs.push_back(CSLD.back().view());
    cs.push_back(CHLD.back().view());

    SLD.push_back(tmv::SymBandMatrix<T,tmv::Lower|tmv::DiagMajor>(
            tmv::SymTriDiagMatrix(v1,v2)));
    HLD.push_back(tmv::SymBandMatrix<T,tmv::Lower|tmv::DiagMajor>(
            tmv::HermTriDiagMatrix(v1,v2,tmv::Lower)));
    CSLD.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower|tmv::DiagMajor>(
            tmv::SymTriDiagMatrix(cv1,cv2)));
    CHLD.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower|tmv::DiagMajor>(
            tmv::HermTriDiagMatrix(cv1r,cv2,tmv::Lower)));
    s.push_back(SLD.back().view());
    s.push_back(HLD.back().view());
    cs.push_back(CSLD.back().view());
    cs.push_back(CHLD.back().view());

    SLC.push_back(tmv::SymBandMatrix<T,tmv::Lower|tmv::ColMajor>(
            DiagMatrixViewOf(v1)));
    HLC.push_back(tmv::HermBandMatrix<T,tmv::Lower|tmv::ColMajor>(
            DiagMatrixViewOf(v1)));
    CSLC.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower|tmv::ColMajor>(
            DiagMatrixViewOf(cv1)));
    CHLC.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower|tmv::ColMajor>(
            DiagMatrixViewOf(cv1)));
    s.push_back(SLC.back().view());
    s.push_back(HLC.back().view());
    cs.push_back(CSLC.back().view());
    cs.push_back(CHLC.back().view());

    M.push_back(a2);
    CM.push_back(ca2);
    s.push_back(SymBandMatrixViewOf(M.back(),tmv::Upper,6).subSymBandMatrix(0,2*N,3,2));
    cs.push_back(SymBandMatrixViewOf(CM.back(),tmv::Upper,6).subSymBandMatrix(0,2*N,3,2));

    M.push_back(a2);
    CM.push_back(ca2);
    CM.back().diag().imagPart().setZero();
    s.push_back(HermBandMatrixViewOf(M.back(),tmv::Upper,6).subSymBandMatrix(0,2*N,3,2));
    cs.push_back(HermBandMatrixViewOf(CM.back(),tmv::Upper,6).subSymBandMatrix(0,2*N,3,2));
#endif
    TMVAssert(int(SUR.size()) <= RESERVE);
    TMVAssert(int(SUC.size()) <= RESERVE);
    TMVAssert(int(SUD.size()) <= RESERVE);
    TMVAssert(int(SLR.size()) <= RESERVE);
    TMVAssert(int(SLC.size()) <= RESERVE);
    TMVAssert(int(SLD.size()) <= RESERVE);
    TMVAssert(int(HUR.size()) <= RESERVE);
    TMVAssert(int(HUC.size()) <= RESERVE);
    TMVAssert(int(HUD.size()) <= RESERVE);
    TMVAssert(int(HLR.size()) <= RESERVE);
    TMVAssert(int(HLC.size()) <= RESERVE);
    TMVAssert(int(HLD.size()) <= RESERVE);
    TMVAssert(int(M.size()) <= RESERVE);
    TMVAssert(int(CSUR.size()) <= RESERVE);
    TMVAssert(int(CSUC.size()) <= RESERVE);
    TMVAssert(int(CSUD.size()) <= RESERVE);
    TMVAssert(int(CSLR.size()) <= RESERVE);
    TMVAssert(int(CSLC.size()) <= RESERVE);
    TMVAssert(int(CSLD.size()) <= RESERVE);
    TMVAssert(int(CHUR.size()) <= RESERVE);
    TMVAssert(int(CHUC.size()) <= RESERVE);
    TMVAssert(int(CHUD.size()) <= RESERVE);
    TMVAssert(int(CHLR.size()) <= RESERVE);
    TMVAssert(int(CHLC.size()) <= RESERVE);
    TMVAssert(int(CHLD.size()) <= RESERVE);
    TMVAssert(int(CM.size()) <= RESERVE);
}

