template <class T> inline void MakeBandList(
    std::vector<tmv::BandMatrixView<T> >& b,
    std::vector<tmv::BandMatrixView<std::complex<T> > >& cb,
    std::vector<tmv::BaseMatrix<T>*>& B,
    std::vector<tmv::BaseMatrix<std::complex<T> >*>& CB)
{
  std::vector<tmv::BandMatrix<T,tmv::RowMajor>*> BR;
  std::vector<tmv::BandMatrix<T,tmv::ColMajor>*> BC;
  std::vector<tmv::BandMatrix<T,tmv::DiagMajor>*> BD;
  std::vector<tmv::Matrix<T>*> M;
  std::vector<tmv::BandMatrix<std::complex<T>,tmv::RowMajor>*> CBR;
  std::vector<tmv::BandMatrix<std::complex<T>,tmv::ColMajor>*> CBC;
  std::vector<tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>*> CBD;
  std::vector<tmv::Matrix<std::complex<T> >*> CM;
  const int N = 10;
  tmv::Matrix<T> a1(N,N);
  tmv::Matrix<std::complex<T> > ca1(N,N);

  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);
  a1.diag().addToAll(T(3*N));
  tmv::Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = T(1-3*i+6*j);
  a2.diag().addToAll(T(3*N));
  tmv::Vector<T> v1(N);
  tmv::Vector<T> v2(N-1);
  for (int i=0; i<N; ++i) v1(i) = T(16-3*i); 
  for (int i=0; i<N-1; ++i) v2(i) = T(-6+i); 

  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
    ca1(i,j) = std::complex<T>(3+i-5*j,4-8*i-j);
  ca1.diag().addToAll(T(3*N));
  tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) 
    ca2(i,j) = std::complex<T>(1-3*i+6*j,8+2*i-6*j);
  ca2.diag().addToAll(T(3*N));
  tmv::Vector<std::complex<T> > cv1(N);
  tmv::Vector<std::complex<T> > cv2(N-1);
  for (int i=0; i<N; ++i) cv1(i) = std::complex<T>(16-3*i,i+4); 
  for (int i=0; i<N-1; ++i) cv2(i) = std::complex<T>(2*i-3,-6+i); 

  BR.push_back(new tmv::BandMatrix<T,tmv::RowMajor>(a1,3,1)); 
  CBR.push_back(new tmv::BandMatrix<std::complex<T>,tmv::RowMajor>(ca1,3,1));
  b.push_back(BR.back()->view());
  cb.push_back(CBR.back()->view());
  BD.push_back(new tmv::BandMatrix<T,tmv::DiagMajor>(a1,3,1));
  CBD.push_back(new tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(ca1,3,1));
  b.push_back(BD.back()->view());
  cb.push_back(CBD.back()->view());
  BD.push_back(new tmv::BandMatrix<T,tmv::DiagMajor>(*BR.back(),1,1));
  CBD.push_back(new tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(*CBR.back(),1,1));
  b.push_back(BD.back()->view());
  cb.push_back(CBD.back()->view());
  BC.push_back(new tmv::BandMatrix<T,tmv::ColMajor>(*BR.back(),3,0));
  CBC.push_back(new tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(*CBR.back(),3,0));
  b.push_back(BC.back()->view());
  cb.push_back(CBC.back()->view());

#ifdef XTEST
  BC.push_back(new tmv::BandMatrix<T,tmv::ColMajor>(a2,6,6)); 
  CBC.push_back(new tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(ca2,6,6));
  b.push_back(BC.back()->SubBandMatrix(0,2*N,0,N,3,3));
  cb.push_back(CBC.back()->SubBandMatrix(0,2*N,0,N,3,3));
  BC.push_back(new tmv::BandMatrix<T,tmv::ColMajor>(*BC.back()));
  CBC.push_back(new tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(*CBC.back()));
  b.push_back(BC.back()->SubBandMatrix(0,N+2,0,N,4,4));
  cb.push_back(CBC.back()->SubBandMatrix(0,N+2,0,N,4,4));
  BC.push_back(new tmv::BandMatrix<T,tmv::ColMajor>(*BC.back()));
  CBC.push_back(new tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(*CBC.back()));
  b.push_back(BC.back()->SubBandMatrix(0,2*N,0,2*N,3,3,2,2));
  cb.push_back(CBC.back()->SubBandMatrix(0,2*N,0,2*N,3,3,2,2));
  BC.push_back(new tmv::BandMatrix<T,tmv::ColMajor>(*BC.back()));
  CBC.push_back(new tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(*CBC.back()));
  b.push_back(BC.back()->SubBandMatrix(0,N,0,2*N,3,3));
  cb.push_back(CBC.back()->SubBandMatrix(0,N,0,2*N,3,3));
  BC.push_back(new tmv::BandMatrix<T,tmv::ColMajor>(*BC.back()));
  CBC.push_back(new tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(*CBC.back())); 
  b.push_back(BC.back()->SubBandMatrix(0,N,0,N+2,4,4));
  cb.push_back(CBC.back()->SubBandMatrix(0,N,0,N+2,4,4));
  BC.push_back(new tmv::BandMatrix<T,tmv::ColMajor>(a1,3,1));
  CBC.push_back(new tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(ca1,3,1));
  b.push_back(BC.back()->view());
  cb.push_back(CBC.back()->view());
  BD.push_back(new tmv::BandMatrix<T,tmv::DiagMajor>(a1,1,3));
  CBD.push_back(new tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(ca1,1,3));
  b.push_back(BD.back()->view());
  cb.push_back(CBD.back()->view());
  BC.push_back(new tmv::BandMatrix<T,tmv::ColMajor>(a1,0,3));
  CBC.push_back(new tmv::BandMatrix<std::complex<T>,tmv::ColMajor>(ca1,0,3));
  b.push_back(BC.back()->view());
  cb.push_back(CBC.back()->view());
  BD.push_back(new tmv::BandMatrix<T,tmv::DiagMajor>(tmv::UpperBiDiagMatrix(v1,v2))); 
  CBD.push_back(new tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(tmv::UpperBiDiagMatrix(cv1,cv2)));
  b.push_back(BD.back()->view());
  cb.push_back(CBD.back()->view());
  BD.push_back(new tmv::BandMatrix<T,tmv::DiagMajor>(tmv::LowerBiDiagMatrix(v2,v1)));
  CBD.push_back(new tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(tmv::LowerBiDiagMatrix(cv2,cv1)));
  b.push_back(BD.back()->view());
  cb.push_back(CBD.back()->view());
  BD.push_back(new tmv::BandMatrix<T,tmv::DiagMajor>(tmv::TriDiagMatrix(v2,v1,v2))); 
  CBD.push_back(new tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(tmv::TriDiagMatrix(cv2,cv1,cv2)));
  b.push_back(BD.back()->view());
  cb.push_back(CBD.back()->view());
  BD.push_back(new tmv::BandMatrix<T,tmv::DiagMajor>(tmv::UpperBiDiagMatrix(v1,v1)));
  CBD.push_back(new tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(tmv::UpperBiDiagMatrix(cv1,cv1)));
  b.push_back(BD.back()->view());
  cb.push_back(CBD.back()->view());
  BD.push_back(new tmv::BandMatrix<T,tmv::DiagMajor>(tmv::LowerBiDiagMatrix(v1,v1)));
  CBD.push_back(new tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(tmv::LowerBiDiagMatrix(cv1,cv1)));
  b.push_back(BD.back()->view());
  cb.push_back(CBD.back()->view());
  BD.push_back(new tmv::BandMatrix<T,tmv::DiagMajor>(tmv::TriDiagMatrix(v1,v1,v2)));
  CBD.push_back(new tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(tmv::TriDiagMatrix(cv1,cv1,cv2)));
  b.push_back(BD.back()->view());
  cb.push_back(CBD.back()->view());
  BD.push_back(new tmv::BandMatrix<T,tmv::DiagMajor>(tmv::TriDiagMatrix(v2,v1,v1)));
  CBD.push_back(new tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(tmv::TriDiagMatrix(cv2,cv1,cv1)));
  b.push_back(BD.back()->view());
  cb.push_back(CBD.back()->view());
  BD.push_back(new tmv::BandMatrix<T,tmv::DiagMajor>(a1,1,N-1));
  CBD.push_back(new tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(ca1,1,N-1));
  b.push_back(BD[10]->view());
  cb.push_back(CBD[10]->view());
  BD.push_back(new tmv::BandMatrix<T,tmv::DiagMajor>(a1,3,N-2));
  CBD.push_back(new tmv::BandMatrix<std::complex<T>,tmv::DiagMajor>(ca1,3,N-2));
  b.push_back(BD[11]->view());
  cb.push_back(CBD[11]->view());
  BR.push_back(new tmv::BandMatrix<T,tmv::RowMajor>(a1,0,0));
  CBR.push_back(new tmv::BandMatrix<std::complex<T>,tmv::RowMajor>(ca1,0,0));
  b.push_back(BR.back()->view());
  cb.push_back(CBR.back()->view());
  M.push_back(new tmv::Matrix<T>(a1));
  CM.push_back(new tmv::Matrix<std::complex<T> >(ca1));
  b.push_back(BandMatrixViewOf(*M.back(),3,N-1));
  cb.push_back(BandMatrixViewOf(*CM.back(),3,N-1));
#endif
  B.insert(B.end(),BR.begin(),BR.end());
  B.insert(B.end(),BC.begin(),BC.end());
  B.insert(B.end(),BD.begin(),BD.end());
  B.insert(B.end(),M.begin(),M.end());
  CB.insert(CB.end(),CBR.begin(),CBR.end());
  CB.insert(CB.end(),CBC.begin(),CBC.end());
  CB.insert(CB.end(),CBD.begin(),CBD.end());
  CB.insert(CB.end(),CM.begin(),CM.end());
}
