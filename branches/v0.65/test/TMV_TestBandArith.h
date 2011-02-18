#define MakeBandList(b,cb) \
    std::vector<tmv::BandMatrix<T,tmv::RowMajor> > B1; \
    std::vector<tmv::BandMatrix<T,tmv::ColMajor> > B2; \
    std::vector<tmv::BandMatrix<T,tmv::DiagMajor> > B3; \
    std::vector<tmv::Matrix<T> > B4; \
    std::vector<tmv::BandMatrix<std::complex<T>,tmv::RowMajor> > CB1; \
    std::vector<tmv::BandMatrix<std::complex<T>,tmv::ColMajor> > CB2; \
    std::vector<tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> > CB3; \
    std::vector<tmv::Matrix<std::complex<T> > > CB4; \
    DoMakeBandList(b,cb,B1,B2,B3,B4,CB1,CB2,CB3,CB4); 

template <class T> 
inline void DoMakeBandList(
    std::vector<tmv::BandMatrixView<T> >& b,
    std::vector<tmv::BandMatrixView<std::complex<T> > >& cb,
    std::vector<tmv::BandMatrix<T,tmv::RowMajor> >& BR,
    std::vector<tmv::BandMatrix<T,tmv::ColMajor> >& BC,
    std::vector<tmv::BandMatrix<T,tmv::DiagMajor> >& BD,
    std::vector<tmv::Matrix<T> >& M,
    std::vector<tmv::BandMatrix<std::complex<T>,tmv::RowMajor> >& CBR,
    std::vector<tmv::BandMatrix<std::complex<T>,tmv::ColMajor> >& CBC,
    std::vector<tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> >& CBD,
    std::vector<tmv::Matrix<std::complex<T> > >& CM)
{
    const int N = std::numeric_limits<T>::is_integer ? 6 : 10;
    const int RESERVE = 20;
    // Need to do this so the push_back's don't invalidate the views that
    // are saved in b and cb.
    BR.reserve(RESERVE);
    BC.reserve(RESERVE);
    BD.reserve(RESERVE);
    M.reserve(RESERVE);
    CBR.reserve(RESERVE);
    CBC.reserve(RESERVE);
    CBD.reserve(RESERVE);
    CM.reserve(RESERVE);

    tmv::Matrix<T> a1(N,N);
    tmv::Matrix<std::complex<T> > ca1(N,N);
    tmv::Matrix<T> a2(2*N,2*N);
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
    } else {
        for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
            a1(i,j) = T(3+i-5*j);
            ca1(i,j) = std::complex<T>(3+i-5*j,4-8*i-j);
        }
        for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) {
            a2(i,j) = T(1-3*i+6*j);
            ca2(i,j) = std::complex<T>(1-3*i+6*j,8+2*i-6*j);
        }
        a1.diag().addToAll(T(3*N));
        a2.diag().addToAll(T(3*N));
        ca1.diag().addToAll(T(3*N));
        ca2.diag().addToAll(T(3*N));
        for (int i=0; i<N; ++i) v1(i) = T(16-3*i); 
        for (int i=0; i<N-1; ++i) v2(i) = T(-6+i); 
        for (int i=0; i<N; ++i) cv1(i) = std::complex<T>(16-3*i,i+4); 
        for (int i=0; i<N-1; ++i) cv2(i) = std::complex<T>(2*i-3,-6+i); 
    }

    BR.push_back(tmv::BandMatrix<T,tmv::RowMajor>(a1,3,1)); 
    CBR.push_back(tmv::BandMatrix<std::complex<T>,tmv::RowMajor>(ca1,3,1));
    b.push_back(BR.back().view());
    cb.push_back(CBR.back().view());
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(a1,3,1));
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(ca1,3,1));
    b.push_back(BD.back().view());
    cb.push_back(CBD.back().view());
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(BR.back(),1,1));
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(CBR.back(),1,1));
    b.push_back(BD.back().view());
    cb.push_back(CBD.back().view());
    BC.push_back(tmv::BandMatrix<T,tmv::ColMajor>(BR.back(),3,0));
    CBC.push_back(tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(CBR.back(),3,0));
    b.push_back(BC.back().view());
    cb.push_back(CBC.back().view());

#if (XTEST & 2)
    BC.push_back(tmv::BandMatrix<T,tmv::ColMajor>(a2,6,6)); 
    CBC.push_back(tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(ca2,6,6));
    b.push_back(BC.back().subBandMatrix(0,2*N,0,N,3,3));
    cb.push_back(CBC.back().subBandMatrix(0,2*N,0,N,3,3));
    BC.push_back(tmv::BandMatrix<T,tmv::ColMajor>(BC.back()));
    CBC.push_back(tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(CBC.back()));
    b.push_back(BC.back().subBandMatrix(0,N+2,0,N,4,4));
    cb.push_back(CBC.back().subBandMatrix(0,N+2,0,N,4,4));
    BC.push_back(tmv::BandMatrix<T,tmv::ColMajor>(BC.back()));
    CBC.push_back(tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(CBC.back()));
    b.push_back(BC.back().subBandMatrix(0,2*N,0,2*N,3,3,2,2));
    cb.push_back(CBC.back().subBandMatrix(0,2*N,0,2*N,3,3,2,2));
    BC.push_back(tmv::BandMatrix<T,tmv::ColMajor>(BC.back()));
    CBC.push_back(tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(CBC.back()));
    b.push_back(BC.back().subBandMatrix(0,N,0,2*N,3,3));
    cb.push_back(CBC.back().subBandMatrix(0,N,0,2*N,3,3));
    BC.push_back(tmv::BandMatrix<T,tmv::ColMajor>(BC.back()));
    CBC.push_back(tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(CBC.back())); 
    b.push_back(BC.back().subBandMatrix(0,N,0,N+2,4,4));
    cb.push_back(CBC.back().subBandMatrix(0,N,0,N+2,4,4));
    BC.push_back(tmv::BandMatrix<T,tmv::ColMajor>(a1,3,1));
    CBC.push_back(tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(ca1,3,1));
    b.push_back(BC.back().view());
    cb.push_back(CBC.back().view());
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(a1,1,3));
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(ca1,1,3));
    b.push_back(BD.back().view());
    cb.push_back(CBD.back().view());
    BC.push_back(tmv::BandMatrix<T,tmv::ColMajor>(a1,0,3));
    CBC.push_back(tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(ca1,0,3));
    b.push_back(BC.back().view());
    cb.push_back(CBC.back().view());
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(tmv::UpperBiDiagMatrix(v1,v2))); 
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(tmv::UpperBiDiagMatrix(cv1,cv2)));
    b.push_back(BD.back().view());
    cb.push_back(CBD.back().view());
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(tmv::LowerBiDiagMatrix(v2,v1)));
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(tmv::LowerBiDiagMatrix(cv2,cv1)));
    b.push_back(BD.back().view());
    cb.push_back(CBD.back().view());
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(tmv::TriDiagMatrix(v2,v1,v2))); 
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(tmv::TriDiagMatrix(cv2,cv1,cv2)));
    b.push_back(BD.back().view());
    cb.push_back(CBD.back().view());
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(tmv::UpperBiDiagMatrix(v1,v1)));
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(tmv::UpperBiDiagMatrix(cv1,cv1)));
    b.push_back(BD.back().view());
    cb.push_back(CBD.back().view());
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(tmv::LowerBiDiagMatrix(v1,v1)));
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(tmv::LowerBiDiagMatrix(cv1,cv1)));
    b.push_back(BD.back().view());
    cb.push_back(CBD.back().view());
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(tmv::TriDiagMatrix(v1,v1,v2)));
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(tmv::TriDiagMatrix(cv1,cv1,cv2)));
    b.push_back(BD.back().view());
    cb.push_back(CBD.back().view());
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(tmv::TriDiagMatrix(v2,v1,v1)));
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(tmv::TriDiagMatrix(cv2,cv1,cv1)));
    b.push_back(BD.back().view());
    cb.push_back(CBD.back().view());
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(a1,1,N-1));
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(ca1,1,N-1));
    b.push_back(BD[10].view());
    cb.push_back(CBD[10].view());
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(a1,3,N-2));
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(ca1,3,N-2));
    b.push_back(BD[11].view());
    cb.push_back(CBD[11].view());
    BR.push_back(tmv::BandMatrix<T,tmv::RowMajor>(a1,0,0));
    CBR.push_back(tmv::BandMatrix<std::complex<T>,tmv::RowMajor>(ca1,0,0));
    b.push_back(BR.back().view());
    cb.push_back(CBR.back().view());
    M.push_back(tmv::Matrix<T>(a1));
    CM.push_back(tmv::Matrix<std::complex<T> >(ca1));
    b.push_back(BandMatrixViewOf(M.back(),3,N-1));
    cb.push_back(BandMatrixViewOf(CM.back(),3,N-1));
#endif

    // Make sure we didn't exceed the reserve amount.  Otherwise the 
    // views in b and cb are invalid.
    TMVAssert(int(BR.size()) <= RESERVE);
    TMVAssert(int(BC.size()) <= RESERVE);
    TMVAssert(int(BD.size()) <= RESERVE);
    TMVAssert(int(M.size()) <= RESERVE);
    TMVAssert(int(CBR.size()) <= RESERVE);
    TMVAssert(int(CBC.size()) <= RESERVE);
    TMVAssert(int(CBD.size()) <= RESERVE);
    TMVAssert(int(CM.size()) <= RESERVE);
}
