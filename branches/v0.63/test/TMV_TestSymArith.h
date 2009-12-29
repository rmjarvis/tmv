// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
template <class T> inline void MakeSymList(
    std::vector<tmv::SymMatrixView<T> >& s,
    std::vector<tmv::SymMatrixView<std::complex<T> > >& cs,
    std::vector<tmv::BaseMatrix<T>*>& B,
    std::vector<tmv::BaseMatrix<std::complex<T> >*>& CB,
    PosDefCode pdc)
{
  const int N=10;

  tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-5*j);
  tmv::Matrix<std::complex<T> > ca1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) = 
    std::complex<T>(3+i-5*j,2-3*i);
  tmv::Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = T(1-3*i+6*j);
  tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) ca2(i,j) = 
    std::complex<T>(1-3*i+6*j,-4+2*j);

  if (pdc == PosDef) {
    a1 /= T(N*N);
    a2 /= T(N*N);
    a1.diag().AddToAll(T(1));
    a2.diag().AddToAll(T(1));
    for(int i=0;i<N;++i) a1.diag()(i) += T(i);
    for(int i=0;i<2*N;++i) a2.diag()(i) += T(i);
    ca1 /= T(N*N);
    ca2 /= T(N*N);
    ca1.diag().AddToAll(T(1));
    ca2.diag().AddToAll(T(1));
    for(int i=0;i<N;++i) ca1.diag()(i) += T(i);
    for(int i=0;i<2*N;++i) ca2.diag()(i) += T(i);
  } else if (pdc == InDef) {
    a1 /= T(N);
    a2 /= T(N);
    a1.diag(0,2,4).Zero();
    a2.diag(0,2,4).Zero();
    a1(4,0) = a1(0,4) = T(50);
    a1.diag(-1,6,9).AddToAll(T(10));
    a2.diag(-1,6,9).AddToAll(T(10));
    a1.diag(1,6,9).AddToAll(T(10));
    a2.diag(1,6,9).AddToAll(T(10));
    a1.row(9,0,5).AddToAll(T(10));
    a2.row(9,0,5).AddToAll(T(10));
    a1(3,3) = T(0);
    a2(3,3) = T(0);
    if (N > 10) {
      a1.diag(0,10,N) *= T(0.0001);
      a2.diag(0,10,N) *= T(0.0001);
    }
    ca1 /= T(N);
    ca2 /= T(N);
    ca1.diag(0,2,4).Zero();
    ca2.diag(0,2,4).Zero();
    ca1(4,0) = ca1(0,4) = T(50);
    ca1.diag(-1,6,9).AddToAll(T(10));
    ca2.diag(-1,6,9).AddToAll(T(10));
    ca1.diag(1,6,9).AddToAll(T(10));
    ca2.diag(1,6,9).AddToAll(T(10));
    ca1.row(9,0,5).AddToAll(T(10));
    ca2.row(9,0,5).AddToAll(T(10));
    ca1(3,3) = T(0);
    ca2(3,3) = T(0);
    if (N > 10) {
      ca1.diag(0,10,N) *= T(0.0001);
      ca2.diag(0,10,N) *= T(0.0001);
    }
  } else {
    a1.row(2).Zero();
    a1.col(2).Zero();
    a1.row(4) = a1.row(5);
    a1.col(4) = a1.col(5);
    a1(4,5) = a1(5,4) = a1(5,5) = a1(4,4);
    a2.row(2).Zero();
    a2.col(2).Zero();
    a2.row(4) = a2.row(5);
    a2.col(4) = a2.col(5);
    a2(4,5) = a2(5,4) = a2(5,5) = a2(4,4);
    ca1.row(2).Zero();
    ca1.col(2).Zero();
    ca1.row(4) = ca1.row(5);
    ca1.col(4) = ca1.col(5);
    ca1(4,5) = ca1(5,4) = ca1(5,5) = ca1(4,4);
    ca2.row(2).Zero();
    ca2.col(2).Zero();
    ca2.row(4) = ca2.row(5);
    ca2.col(4) = ca2.col(5);
    ca2(4,5) = ca2(5,4) = ca2(5,5) = ca2(4,4);
  }

  tmv::SymMatrix<T,tmv::Upper,tmv::RowMajor>* SUR = new
  tmv::SymMatrix<T,tmv::Upper,tmv::RowMajor>(a1);
  tmv::HermMatrix<T,tmv::Upper,tmv::RowMajor>* HUR = new
  tmv::HermMatrix<T,tmv::Upper,tmv::RowMajor>(a1);
  tmv::SymMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor>* CSUR = new
  tmv::SymMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor>(ca1);
  tmv::HermMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor>* CHUR = new
  tmv::HermMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor>(ca1);
  s.push_back(SUR->View());
  s.push_back(HUR->View());
  cs.push_back(CSUR->View());
  cs.push_back(CHUR->View());
  B.push_back(SUR);
  B.push_back(HUR);
  CB.push_back(CSUR);
  CB.push_back(CHUR);

  tmv::SymMatrix<T,tmv::Upper,tmv::ColMajor>* SUC = new
  tmv::SymMatrix<T,tmv::Upper,tmv::ColMajor>(a1);
  tmv::HermMatrix<T,tmv::Upper,tmv::ColMajor>* HUC = new
  tmv::HermMatrix<T,tmv::Upper,tmv::ColMajor>(a1);
  tmv::SymMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor>* CSUC = new
  tmv::SymMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor>(ca1);
  tmv::HermMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor>* CHUC = new
  tmv::HermMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor>(ca1);
  s.push_back(SUC->View());
  s.push_back(HUC->View());
  cs.push_back(CSUC->View());
  cs.push_back(CHUC->View());
  B.push_back(SUC);
  B.push_back(HUC);
  CB.push_back(CSUC);
  CB.push_back(CHUC);

#ifdef XTEST
  tmv::SymMatrix<T,tmv::Lower,tmv::RowMajor>* SLR = new
  tmv::SymMatrix<T,tmv::Lower,tmv::RowMajor>(a1);
  tmv::HermMatrix<T,tmv::Lower,tmv::RowMajor>* HLR = new
  tmv::HermMatrix<T,tmv::Lower,tmv::RowMajor>(a1);
  tmv::SymMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor>* CSLR = new
  tmv::SymMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor>(ca1);
  tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor>* CHLR = new
  tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor>(ca1);
  s.push_back(SLR->View());
  s.push_back(HLR->View());
  cs.push_back(CSLR->View());
  cs.push_back(CHLR->View());
  B.push_back(SLR);
  B.push_back(HLR);
  CB.push_back(CSLR);
  CB.push_back(CHLR);
  tmv::SymMatrix<T,tmv::Lower,tmv::ColMajor>* SLC = new
  tmv::SymMatrix<T,tmv::Lower,tmv::ColMajor>(a1);
  tmv::HermMatrix<T,tmv::Lower,tmv::ColMajor>* HLC = new
  tmv::HermMatrix<T,tmv::Lower,tmv::ColMajor>(a1);
  tmv::SymMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>* CSLC = new
  tmv::SymMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(ca1);
  tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>* CHLC = new
  tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(ca1);
  s.push_back(SLC->View());
  s.push_back(HLC->View());
  cs.push_back(CSLC->View());
  cs.push_back(CHLC->View());
  B.push_back(SLC);
  B.push_back(HLC);
  CB.push_back(CSLC);
  CB.push_back(CHLC);
  tmv::Matrix<T>* a2a = new tmv::Matrix<T>(a2);
  tmv::Matrix<T>* a2b = new tmv::Matrix<T>(a2);
  tmv::Matrix<std::complex<T> >* ca2a = new tmv::Matrix<std::complex<T> >(ca2);
  tmv::Matrix<std::complex<T> >* ca2b = new tmv::Matrix<std::complex<T> >(ca2);
  ca2b->diag().Imag().Zero();
  s.push_back(SymMatrixViewOf(*a2a,tmv::Upper).SubSymMatrix(0,2*N,2));
  s.push_back(HermMatrixViewOf(*a2b,tmv::Upper).SubSymMatrix(0,2*N,2));
  cs.push_back(SymMatrixViewOf(*ca2a,tmv::Upper).SubSymMatrix(0,2*N,2));
  cs.push_back(HermMatrixViewOf(*ca2b,tmv::Upper).SubSymMatrix(0,2*N,2));
  B.push_back(a2a);
  B.push_back(a2b);
  CB.push_back(ca2a);
  CB.push_back(ca2b);
#endif
}


