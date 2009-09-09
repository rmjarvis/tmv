template <class T> inline void MakeBandList(
    std::vector<tmv::BandMatrixView<T> >& b,
    std::vector<tmv::BandMatrixView<std::complex<T> > >& cb)
{
  static std::vector<tmv::BandMatrix<T,tmv::RowMajor> > BR;
  static std::vector<tmv::BandMatrix<T,tmv::ColMajor> > BC;
  static std::vector<tmv::BandMatrix<T,tmv::DiagMajor> > BD;
  static std::vector<tmv::BandMatrix<std::complex<T>,tmv::RowMajor> > CBR;
  static std::vector<tmv::BandMatrix<std::complex<T>,tmv::ColMajor> > CBC;
  static std::vector<tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> > CBD;
  const int N = 10;
  static tmv::Matrix<T> a1(N,N);
  static tmv::Matrix<std::complex<T> > ca1(N,N);

  if (BR.size() == 0) {

    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-5*j;
    a1.diag().AddToAll(T(3*N));
    tmv::Matrix<T> a2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-3*i+6*j;
    a2.diag().AddToAll(T(3*N));
    tmv::Vector<T> v1(N);
    tmv::Vector<T> v2(N-1);
    for (int i=0; i<N; ++i) v1(i) = 16.-3*i; 
    for (int i=0; i<N-1; ++i) v2(i) = -6.+i; 

    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
      ca1(i,j) = std::complex<T>(3.+i-5*j,4.-8*i-j);
    ca1.diag().AddToAll(T(3*N));
    tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) 
      ca2(i,j) = std::complex<T>(1.-3*i+6*j,8.+2*i-6*j);
    ca2.diag().AddToAll(T(3*N));
    tmv::Vector<std::complex<T> > cv1(N);
    tmv::Vector<std::complex<T> > cv2(N-1);
    for (int i=0; i<N; ++i) cv1(i) = std::complex<T>(16.-3*i,i+4.); 
    for (int i=0; i<N-1; ++i) cv2(i) = std::complex<T>(2*i-3.,-6.+i); 

    BR.push_back(tmv::BandMatrix<T,tmv::RowMajor>(a1,3,1)); // 0
    CBR.push_back(tmv::BandMatrix<std::complex<T>,tmv::RowMajor>(ca1,3,1));
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(a1,3,1)); // 0
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(ca1,3,1));
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(BR[0],1,1)); // 1
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(CBR[0],1,1));
    BC.push_back(tmv::BandMatrix<T,tmv::ColMajor>(BR[0],3,0)); // 0
    CBC.push_back(tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(CBR[0],3,0));
#ifdef XTEST
    BC.push_back(tmv::BandMatrix<T,tmv::ColMajor>(a2,6,6)); // 1
    CBC.push_back(tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(ca2,6,6));
    BC.push_back(BC[1]); // 2
    CBC.push_back(CBC[1]);
    BC.push_back(BC[1]); // 3
    CBC.push_back(CBC[1]);
    BC.push_back(BC[1]); // 4
    CBC.push_back(CBC[1]);
    BC.push_back(BC[1]); // 5
    CBC.push_back(CBC[1]); 
    BC.push_back(tmv::BandMatrix<T,tmv::ColMajor>(a1,3,1)); // 6
    CBC.push_back(tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(ca1,3,1));
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(a1,1,3)); // 2
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(ca1,1,3));
    BC.push_back(tmv::BandMatrix<T,tmv::ColMajor>(a1,0,3)); // 7
    CBC.push_back(tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(ca1,0,3));
    BD.push_back(tmv::UpperBiDiagMatrix(v1,v2)); // 3
    CBD.push_back(tmv::UpperBiDiagMatrix(cv1,cv2));
    BD.push_back(tmv::LowerBiDiagMatrix(v2,v1)); // 4
    CBD.push_back(tmv::LowerBiDiagMatrix(cv2,cv1));
    BD.push_back(tmv::TriDiagMatrix(v2,v1,v2)); // 5
    CBD.push_back(tmv::TriDiagMatrix(cv2,cv1,cv2));
    BD.push_back(tmv::UpperBiDiagMatrix(v1,v1)); // 6
    CBD.push_back(tmv::UpperBiDiagMatrix(cv1,cv1));
    BD.push_back(tmv::LowerBiDiagMatrix(v1,v1)); // 7
    CBD.push_back(tmv::LowerBiDiagMatrix(cv1,cv1));
    BD.push_back(tmv::TriDiagMatrix(v1,v1,v2)); // 8
    CBD.push_back(tmv::TriDiagMatrix(cv1,cv1,cv2));
    BD.push_back(tmv::TriDiagMatrix(v2,v1,v1)); // 9
    CBD.push_back(tmv::TriDiagMatrix(cv2,cv1,cv1));
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(a1,1,N-1)); // 10
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(ca1,1,N-1));
    BD.push_back(tmv::BandMatrix<T,tmv::DiagMajor>(a1,3,N-2)); // 11
    CBD.push_back(tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(ca1,3,N-2));
    BR.push_back(tmv::BandMatrix<T>(a1,0,0)); // 1
    CBR.push_back(tmv::BandMatrix<std::complex<T> >(ca1,0,0));
#endif
  }

  b.push_back(BR[0].View());
  cb.push_back(CBR[0].View());
  b.push_back(BD[0].View());
  cb.push_back(CBD[0].View());
  b.push_back(BD[1].View());
  cb.push_back(CBD[1].View());
  b.push_back(BC[0].View());
  cb.push_back(CBC[0].View());
#ifdef XTEST
  b.push_back(BC[1].SubBandMatrix(0,2*N,0,N,3,3));
  cb.push_back(CBC[1].SubBandMatrix(0,2*N,0,N,3,3));
  b.push_back(BC[2].SubBandMatrix(0,N+2,0,N,4,4));
  cb.push_back(CBC[2].SubBandMatrix(0,N+2,0,N,4,4));
  b.push_back(BC[3].SubBandMatrix(0,2*N,0,2*N,3,3,2,2));
  cb.push_back(CBC[3].SubBandMatrix(0,2*N,0,2*N,3,3,2,2));
  b.push_back(BC[4].SubBandMatrix(0,N,0,2*N,3,3));
  cb.push_back(CBC[4].SubBandMatrix(0,N,0,2*N,3,3));
  b.push_back(BC[5].SubBandMatrix(0,N,0,N+2,4,4));
  cb.push_back(CBC[5].SubBandMatrix(0,N,0,N+2,4,4));
  b.push_back(BC[6].View());
  cb.push_back(CBC[6].View());
  b.push_back(BD[2].View());
  cb.push_back(CBD[2].View());
  b.push_back(BC[7].View());
  cb.push_back(CBC[7].View());
  b.push_back(BD[3].View());
  cb.push_back(CBD[3].View());
  b.push_back(BD[4].View());
  cb.push_back(CBD[4].View());
  b.push_back(BD[5].View());
  cb.push_back(CBD[5].View());
  b.push_back(BD[6].View());
  cb.push_back(CBD[6].View());
  b.push_back(BD[7].View());
  cb.push_back(CBD[7].View());
  b.push_back(BD[8].View());
  cb.push_back(CBD[8].View());
  b.push_back(BD[9].View());
  cb.push_back(CBD[9].View());
  b.push_back(BD[10].View());
  cb.push_back(CBD[10].View());
  b.push_back(BD[11].View());
  cb.push_back(CBD[11].View());
  b.push_back(BR[1].View());
  cb.push_back(CBR[1].View());
  b.push_back(BandMatrixViewOf(a1,3,N-1));
  cb.push_back(BandMatrixViewOf(ca1,3,N-1));
#endif
}
