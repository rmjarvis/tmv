
template <class T>
inline void CopyBackM(
    const tmv::BandMatrix<std::complex<T> >& m0,
    tmv::SymBandMatrixView<std::complex<T> >& m)
{ 
    if (m.issym()) m = SymBandMatrixViewOf(m0,m.uplo()); 
    else m = HermBandMatrixViewOf(m0,m.uplo());
}



template <class T> inline void MakeSymBandList(
    std::vector<tmv::SymBandMatrixView<T> >& s,
    std::vector<tmv::SymBandMatrixView<std::complex<T> > >& cs,
    std::vector<tmv::BaseMatrix<T>*>& B,
    std::vector<tmv::BaseMatrix<std::complex<T> >*>& CB,
    PosDefCode pdc)
{
    // The integer det function has trouble with the 10x10 determinant,
    // since it overflows.  6x6 is ok.
    const int N=std::numeric_limits<T>::is_integer ? 6 : 10;

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

    tmv::SymBandMatrix<T,tmv::Upper,tmv::RowMajor>* SUR = new
        tmv::SymBandMatrix<T,tmv::Upper,tmv::RowMajor>(a1,3);
    tmv::HermBandMatrix<T,tmv::Upper,tmv::RowMajor>* HUR = new
        tmv::HermBandMatrix<T,tmv::Upper,tmv::RowMajor>(a1,3);
    tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor>* CSUR = new
        tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor>(ca1,3);
    tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor>* CHUR = new
        tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor>(ca1,3);
    s.push_back(SUR->view());
    s.push_back(HUR->view());
    cs.push_back(CSUR->view());
    cs.push_back(CHUR->view());
    B.push_back(SUR);
    B.push_back(HUR);
    CB.push_back(CSUR);
    CB.push_back(CHUR);

    tmv::SymBandMatrix<T,tmv::Upper,tmv::ColMajor>* SUC = new
        tmv::SymBandMatrix<T,tmv::Upper,tmv::ColMajor>(a1,3);
    tmv::HermBandMatrix<T,tmv::Upper,tmv::ColMajor>* HUC = new
        tmv::HermBandMatrix<T,tmv::Upper,tmv::ColMajor>(a1,3);
    tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor>* CSUC = new
        tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor>(ca1,3);
    tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor>* CHUC = new
        tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor>(ca1,3);
    s.push_back(SUC->view());
    s.push_back(HUC->view());
    cs.push_back(CSUC->view());
    cs.push_back(CHUC->view());
    B.push_back(SUC);
    B.push_back(HUC);
    CB.push_back(CSUC);
    CB.push_back(CHUC);

    tmv::SymBandMatrix<T,tmv::Upper,tmv::DiagMajor>* SUD1 = new
        tmv::SymBandMatrix<T,tmv::Upper,tmv::DiagMajor>(a1,3);
    tmv::HermBandMatrix<T,tmv::Upper,tmv::DiagMajor>* HUD1 = new
        tmv::HermBandMatrix<T,tmv::Upper,tmv::DiagMajor>(a1,3);
    tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>* CSUD1 = new
        tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(ca1,3);
    tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>* CHUD1 = new
        tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(ca1,3);
    s.push_back(SUD1->view());
    s.push_back(HUD1->view());
    cs.push_back(CSUD1->view());
    cs.push_back(CHUD1->view());
    B.push_back(SUD1);
    B.push_back(HUD1);
    CB.push_back(CSUD1);
    CB.push_back(CHUD1);

    tmv::SymBandMatrix<T,tmv::Upper,tmv::DiagMajor>* SUD2 = new
        tmv::SymBandMatrix<T,tmv::Upper,tmv::DiagMajor>(
            tmv::SymTriDiagMatrix(v1,v2));
    tmv::HermBandMatrix<T,tmv::Upper,tmv::DiagMajor>* HUD2 = new
        tmv::HermBandMatrix<T,tmv::Upper,tmv::DiagMajor>(
            tmv::HermTriDiagMatrix(v1,v2));
    tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>* CSUD2 = new
        tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(
            tmv::SymTriDiagMatrix(cv1,cv2));
    tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>* CHUD2 = new
        tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(
            tmv::HermTriDiagMatrix(v1,cv2,tmv::Upper));
    s.push_back(SUD2->view());
    s.push_back(HUD2->view());
    cs.push_back(CSUD2->view());
    cs.push_back(CHUD2->view());
    B.push_back(SUD2);
    B.push_back(HUD2);
    CB.push_back(CSUD2);
    CB.push_back(CHUD2);

    tmv::SymBandMatrix<T,tmv::Upper,tmv::DiagMajor>* SUD3 = new
        tmv::SymBandMatrix<T,tmv::Upper,tmv::DiagMajor>(
            DiagMatrixViewOf(v1));
    tmv::HermBandMatrix<T,tmv::Upper,tmv::DiagMajor>* HUD3 = new
        tmv::HermBandMatrix<T,tmv::Upper,tmv::DiagMajor>(
            DiagMatrixViewOf(v1));
    tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>* CSUD3 = new
        tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(
            DiagMatrixViewOf(cv1));
    tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>* CHUD3 = new
        tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(
            DiagMatrixViewOf(cv1));
    s.push_back(SUD3->view());
    s.push_back(HUD3->view());
    cs.push_back(CSUD3->view());
    cs.push_back(CHUD3->view());
    B.push_back(SUD3);
    B.push_back(HUD3);
    CB.push_back(CSUD3);
    CB.push_back(CHUD3);

#if (XTEST & 2)
    tmv::SymBandMatrix<T,tmv::Lower,tmv::ColMajor>* SLC = new
        tmv::SymBandMatrix<T,tmv::Lower,tmv::ColMajor>(a1,N-2);
    tmv::HermBandMatrix<T,tmv::Lower,tmv::ColMajor>* HLC = new
        tmv::HermBandMatrix<T,tmv::Lower,tmv::ColMajor>(a1,N-2);
    tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>* CSLC = new
        tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(ca1,N-2);
    tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>* CHLC = new
        tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(ca1,N-2);
    s.push_back(SLC->view());
    s.push_back(HLC->view());
    cs.push_back(CSLC->view());
    cs.push_back(CHLC->view());
    B.push_back(SLC);
    B.push_back(HLC);
    CB.push_back(CSLC);
    CB.push_back(CHLC);

    tmv::SymBandMatrix<T,tmv::Lower,tmv::RowMajor>* SLR = new
        tmv::SymBandMatrix<T,tmv::Lower,tmv::RowMajor>(a1,N-2);
    tmv::HermBandMatrix<T,tmv::Lower,tmv::RowMajor>* HLR = new
        tmv::HermBandMatrix<T,tmv::Lower,tmv::RowMajor>(a1,N-2);
    tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor>* CSLR = new
        tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor>(ca1,N-2);
    tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor>* CHLR = new
        tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor>(ca1,N-2);
    s.push_back(SLR->view());
    s.push_back(HLR->view());
    cs.push_back(CSLR->view());
    cs.push_back(CHLR->view());
    B.push_back(SLR);
    B.push_back(HLR);
    CB.push_back(CSLR);
    CB.push_back(CHLR);

    tmv::SymBandMatrix<T,tmv::Lower,tmv::DiagMajor>* SLD1 = new
        tmv::SymBandMatrix<T,tmv::Lower,tmv::DiagMajor>(a1,N-2);
    tmv::HermBandMatrix<T,tmv::Lower,tmv::DiagMajor>* HLD1 = new
        tmv::HermBandMatrix<T,tmv::Lower,tmv::DiagMajor>(a1,N-2);
    tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor>* CSLD1 = new
        tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor>(ca1,N-2);
    tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor>* CHLD1 = new
        tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor>(ca1,N-2);
    s.push_back(SLD1->view());
    s.push_back(HLD1->view());
    cs.push_back(CSLD1->view());
    cs.push_back(CHLD1->view());
    B.push_back(SLD1);
    B.push_back(HLD1);
    CB.push_back(CSLD1);
    CB.push_back(CHLD1);

    tmv::SymBandMatrix<T,tmv::Lower,tmv::DiagMajor>* SLD2 = new
        tmv::SymBandMatrix<T,tmv::Lower,tmv::DiagMajor>(
            tmv::SymTriDiagMatrix(v1,v2));
    tmv::HermBandMatrix<T,tmv::Lower,tmv::DiagMajor>* HLD2 = new
        tmv::HermBandMatrix<T,tmv::Lower,tmv::DiagMajor>(
            tmv::HermTriDiagMatrix(v1,v2));
    tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor>* CSLD2 = new
        tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor>(
            tmv::SymTriDiagMatrix(cv1,cv2));
    tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor>* CHLD2 = new
        tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor>(
            tmv::HermTriDiagMatrix(cv1r,cv2,tmv::Lower));
    s.push_back(SLD2->view());
    s.push_back(HLD2->view());
    cs.push_back(CSLD2->view());
    cs.push_back(CHLD2->view());
    B.push_back(SLD2);
    B.push_back(HLD2);
    CB.push_back(CSLD2);
    CB.push_back(CHLD2);

    tmv::SymBandMatrix<T,tmv::Lower,tmv::ColMajor>* SLC2 = new
        tmv::SymBandMatrix<T,tmv::Lower,tmv::ColMajor>(
            DiagMatrixViewOf(v1));
    tmv::HermBandMatrix<T,tmv::Lower,tmv::ColMajor>* HLC2 = new
        tmv::HermBandMatrix<T,tmv::Lower,tmv::ColMajor>(
            DiagMatrixViewOf(v1));
    tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>* CSLC2 = new
        tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(
            DiagMatrixViewOf(cv1));
    tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>* CHLC2 = new
        tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(
            DiagMatrixViewOf(cv1));
    s.push_back(SLC2->view());
    s.push_back(SLC2->view());
    cs.push_back(CSLC2->view());
    cs.push_back(CSLC2->view());
    B.push_back(SLC2);
    B.push_back(HLC2);
    CB.push_back(CSLC2);
    CB.push_back(CHLC2);

    tmv::Matrix<T>* a2a = new tmv::Matrix<T>(a2);
    tmv::Matrix<T>* a2b = new tmv::Matrix<T>(a2);
    tmv::Matrix<std::complex<T> >* ca2a = new tmv::Matrix<std::complex<T> >(ca2);
    tmv::Matrix<std::complex<T> >* ca2b = new tmv::Matrix<std::complex<T> >(ca2);
    ca2b->diag().imagPart().setZero();
    s.push_back(SymBandMatrixViewOf(*a2a,tmv::Upper,6).subSymBandMatrix(0,2*N,3,2));
    s.push_back(HermBandMatrixViewOf(*a2b,tmv::Upper,6).subSymBandMatrix(0,2*N,3,2));
    cs.push_back(SymBandMatrixViewOf(*ca2a,tmv::Upper,6).subSymBandMatrix(0,2*N,3,2));
    cs.push_back(HermBandMatrixViewOf(*ca2b,tmv::Upper,6).subSymBandMatrix(0,2*N,3,2));
    B.push_back(a2a);
    B.push_back(a2b);
    CB.push_back(ca2a);
    CB.push_back(ca2b);
#endif
}

