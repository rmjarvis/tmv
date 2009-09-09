template <class T> inline void MakeSymBandList_PosDef(
  std::vector<tmv::SymBandMatrixView<T> >& s,
  std::vector<tmv::SymBandMatrixView<std::complex<T> > >& cs)
{
  static std::vector<tmv::SymBandMatrix<T,tmv::Upper,tmv::RowMajor> > SUR;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> > CSUR;
  static std::vector<tmv::SymBandMatrix<T,tmv::Upper,tmv::ColMajor> > SUC;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> > CSUC;
  static std::vector<tmv::SymBandMatrix<T,tmv::Upper,tmv::DiagMajor> > SUD;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor> > CSUD;
  static std::vector<tmv::HermBandMatrix<T,tmv::Upper,tmv::RowMajor> > HUR;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> > CHUR;
  static std::vector<tmv::HermBandMatrix<T,tmv::Upper,tmv::ColMajor> > HUC;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> > CHUC;
  static std::vector<tmv::HermBandMatrix<T,tmv::Upper,tmv::DiagMajor> > HUD;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor> > CHUD;
  static std::vector<tmv::SymBandMatrix<T,tmv::Lower,tmv::RowMajor> > SLR;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> > CSLR;
  static std::vector<tmv::SymBandMatrix<T,tmv::Lower,tmv::ColMajor> > SLC;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> > CSLC;
  static std::vector<tmv::SymBandMatrix<T,tmv::Lower,tmv::DiagMajor> > SLD;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor> > CSLD;
  static std::vector<tmv::HermBandMatrix<T,tmv::Lower,tmv::RowMajor> > HLR;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> > CHLR;
  static std::vector<tmv::HermBandMatrix<T,tmv::Lower,tmv::ColMajor> > HLC;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> > CHLC;
  static std::vector<tmv::HermBandMatrix<T,tmv::Lower,tmv::DiagMajor> > HLD;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor> > CHLD;
  const int N=10;
  static tmv::Matrix<T> a2a(2*N,2*N);
  static tmv::Matrix<T> a2b(2*N,2*N);
  static tmv::Matrix<std::complex<T> > ca2a(2*N,2*N);
  static tmv::Matrix<std::complex<T> > ca2b(2*N,2*N);

  if (SUR.size() == 0) {

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i+5*j;
    tmv::Matrix<T> a2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.+3*i+6*j;
    tmv::Matrix<std::complex<T> > ca1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) =
      std::complex<T>(3.+i+5.*j,2.+3.*i);
    tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) ca2(i,j) =
      std::complex<T>(1.+3*i+6.*j,4.+2.*j);

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

    tmv::Vector<T> v1(N);
    tmv::Vector<T> v2(N-1);
    for (int i=0; i<N; ++i) v1(i) = 16.+3*i;
    for (int i=0; i<N-1; ++i) v2(i) = +6.+i;

    tmv::Vector<std::complex<T> > cv1(N);
    tmv::Vector<std::complex<T> > cv2(N-1);
    for (int i=0; i<N; ++i) cv1(i) = std::complex<T>(16.+3*i,i+4.);
    for (int i=0; i<N-1; ++i) cv2(i) = std::complex<T>(2*i+3.,+6.+i);
    tmv::Vector<std::complex<T> > cv1r = cv1;
    cv1r.Imag().Zero();

    SUR.push_back(tmv::SymBandMatrix<T,tmv::Upper,tmv::RowMajor>(a1,3));
    HUR.push_back(tmv::HermBandMatrix<T,tmv::Upper,tmv::RowMajor>(a1,3));
    CSUR.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor>(ca1,3));
    CHUR.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor>(ca1,3));

    SUC.push_back(tmv::SymBandMatrix<T,tmv::Upper,tmv::ColMajor>(a1,3));
    HUC.push_back(tmv::HermBandMatrix<T,tmv::Upper,tmv::ColMajor>(a1,3));
    CSUC.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor>(ca1,3));
    CHUC.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor>(ca1,3));

    SUD.push_back(tmv::SymBandMatrix<T,tmv::Upper,tmv::DiagMajor>(a1,3));
    HUD.push_back(tmv::HermBandMatrix<T,tmv::Upper,tmv::DiagMajor>(a1,3));
    CSUD.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(ca1,3));
    CHUD.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(ca1,3));

    SUD.push_back(tmv::SymTriDiagMatrix<tmv::Upper>(v1,v2));
    HUD.push_back(tmv::HermTriDiagMatrix<tmv::Upper>(v1,v2));
    CSUD.push_back(tmv::SymTriDiagMatrix<tmv::Upper>(cv1,cv2));
    CHUD.push_back(tmv::HermTriDiagMatrix<tmv::Upper>(cv1r,cv2));

    SUD.push_back(tmv::SymBandMatrix<T,tmv::Upper,tmv::DiagMajor>(a1,0));
    HUD.push_back(tmv::HermBandMatrix<T,tmv::Upper,tmv::DiagMajor>(a1,0));
    CSUD.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(ca1,0));
    CHUD.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(ca1,0));

#ifdef XTEST
    SLC.push_back(tmv::SymBandMatrix<T,tmv::Lower,tmv::ColMajor>(a1,7));
    HLC.push_back(tmv::HermBandMatrix<T,tmv::Lower,tmv::ColMajor>(a1,7));
    CSLC.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(ca1,7));
    CHLC.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(ca1,7));

    SLR.push_back(tmv::SymBandMatrix<T,tmv::Lower,tmv::RowMajor>(a1,7));
    HLR.push_back(tmv::HermBandMatrix<T,tmv::Lower,tmv::RowMajor>(a1,7));
    CSLR.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor>(ca1,7));
    CHLR.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor>(ca1,7));

    SLD.push_back(tmv::SymBandMatrix<T,tmv::Lower,tmv::DiagMajor>(a1,7));
    HLD.push_back(tmv::HermBandMatrix<T,tmv::Lower,tmv::DiagMajor>(a1,7));
    CSLD.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor>(ca1,7));
    CHLD.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor>(ca1,7));

    SLD.push_back(tmv::SymTriDiagMatrix<tmv::Lower>(v1,v2));
    HLD.push_back(tmv::HermTriDiagMatrix<tmv::Lower>(v1,v2));
    CSLD.push_back(tmv::SymTriDiagMatrix<tmv::Lower>(cv1,cv2));
    CHLD.push_back(tmv::HermTriDiagMatrix<tmv::Lower>(cv1r,cv2));

    SLC.push_back(tmv::SymBandMatrix<T,tmv::Lower,tmv::ColMajor>(a1,0));
    HLC.push_back(tmv::HermBandMatrix<T,tmv::Lower,tmv::ColMajor>(a1,0));
    CSLC.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(ca1,0));
    CHLC.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(ca1,0));

    a2a = a2;
    a2b = a2;
    ca2a = ca2;
    ca2b = ca2;
    ca2b.diag().Imag().Zero();
#endif
  }

  s.push_back(SUR[0].View());
  s.push_back(HUR[0].View());
  cs.push_back(CSUR[0].View());
  cs.push_back(CHUR[0].View());
  s.push_back(SUC[0].View());
  s.push_back(HUC[0].View());
  cs.push_back(CSUC[0].View());
  cs.push_back(CHUC[0].View());
  s.push_back(SUD[0].View());
  s.push_back(HUD[0].View());
  cs.push_back(CSUD[0].View());
  cs.push_back(CHUD[0].View());
  s.push_back(SUD[1].View());
  s.push_back(HUD[1].View());
  cs.push_back(CSUD[1].View());
  cs.push_back(CHUD[1].View());
  s.push_back(SUD[2].View());
  s.push_back(HUD[2].View());
  cs.push_back(CSUD[2].View());
  cs.push_back(CHUD[2].View());
#ifdef XTEST
  s.push_back(SLC[0].View());
  s.push_back(HLC[0].View());
  cs.push_back(CSLC[0].View());
  cs.push_back(CHLC[0].View());
  s.push_back(SLR[0].View());
  s.push_back(HLR[0].View());
  cs.push_back(CSLR[0].View());
  cs.push_back(CHLR[0].View());
  s.push_back(SLD[0].View());
  s.push_back(HLD[0].View());
  cs.push_back(CSLD[0].View());
  cs.push_back(CHLD[0].View());
  s.push_back(SLD[1].View());
  s.push_back(HLD[1].View());
  cs.push_back(CSLD[1].View());
  cs.push_back(CHLD[1].View());
  s.push_back(SLC[1].View());
  s.push_back(HLC[1].View());
  cs.push_back(CSLC[1].View());
  cs.push_back(CHLC[1].View());
  s.push_back(SymBandMatrixViewOf(a2a,tmv::Upper,6).SubSymBandMatrix(0,2*N,3,2));
  s.push_back(HermBandMatrixViewOf(a2b,tmv::Upper,6).SubSymBandMatrix(0,2*N,3,2));
  cs.push_back(SymBandMatrixViewOf(ca2a,tmv::Upper,6).SubSymBandMatrix(0,2*N,3,2));
  cs.push_back(HermBandMatrixViewOf(ca2b,tmv::Upper,6).SubSymBandMatrix(0,2*N,3,2));
#endif
}

template <class T> inline void MakeSymBandList_InDef(
  std::vector<tmv::SymBandMatrixView<T> >& s,
  std::vector<tmv::SymBandMatrixView<std::complex<T> > >& cs)
{
  static std::vector<tmv::SymBandMatrix<T,tmv::Upper,tmv::RowMajor> > SUR;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> > CSUR;
  static std::vector<tmv::SymBandMatrix<T,tmv::Upper,tmv::ColMajor> > SUC;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> > CSUC;
  static std::vector<tmv::SymBandMatrix<T,tmv::Upper,tmv::DiagMajor> > SUD;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor> > CSUD;
  static std::vector<tmv::HermBandMatrix<T,tmv::Upper,tmv::RowMajor> > HUR;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> > CHUR;
  static std::vector<tmv::HermBandMatrix<T,tmv::Upper,tmv::ColMajor> > HUC;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> > CHUC;
  static std::vector<tmv::HermBandMatrix<T,tmv::Upper,tmv::DiagMajor> > HUD;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor> > CHUD;
  static std::vector<tmv::SymBandMatrix<T,tmv::Lower,tmv::RowMajor> > SLR;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> > CSLR;
  static std::vector<tmv::SymBandMatrix<T,tmv::Lower,tmv::ColMajor> > SLC;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> > CSLC;
  static std::vector<tmv::SymBandMatrix<T,tmv::Lower,tmv::DiagMajor> > SLD;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor> > CSLD;
  static std::vector<tmv::HermBandMatrix<T,tmv::Lower,tmv::RowMajor> > HLR;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> > CHLR;
  static std::vector<tmv::HermBandMatrix<T,tmv::Lower,tmv::ColMajor> > HLC;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> > CHLC;
  static std::vector<tmv::HermBandMatrix<T,tmv::Lower,tmv::DiagMajor> > HLD;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor> > CHLD;
  const int N=10;
  static tmv::Matrix<T> a2a(2*N,2*N);
  static tmv::Matrix<T> a2b(2*N,2*N);
  static tmv::Matrix<std::complex<T> > ca2a(2*N,2*N);
  static tmv::Matrix<std::complex<T> > ca2b(2*N,2*N);

  if (SUR.size() == 0) {

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-5*j;
    tmv::Matrix<T> a2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-3*i+6*j;
    tmv::Matrix<std::complex<T> > ca1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) =
      std::complex<T>(3.+i-5.*j,2.-3.*i);
    tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) ca2(i,j) =
      std::complex<T>(1.-3*i+6.*j,-4.+2.*j);

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
    ca1(4,0) = a1(0,4) = T(50);
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

    tmv::Vector<T> v1(N);
    tmv::Vector<T> v2(N-1);
    for (int i=0; i<N; ++i) v1(i) = 16.-3*i;
    for (int i=0; i<N-1; ++i) v2(i) = -6.+i;

    tmv::Vector<std::complex<T> > cv1(N);
    tmv::Vector<std::complex<T> > cv2(N-1);
    for (int i=0; i<N; ++i) cv1(i) = std::complex<T>(16.-3*i,i+4.);
    for (int i=0; i<N-1; ++i) cv2(i) = std::complex<T>(2*i-3.,-6.+i);
    tmv::Vector<std::complex<T> > cv1r = cv1;
    cv1r.Imag().Zero();

    SUR.push_back(tmv::SymBandMatrix<T,tmv::Upper,tmv::RowMajor>(a1,3));
    HUR.push_back(tmv::HermBandMatrix<T,tmv::Upper,tmv::RowMajor>(a1,3));
    CSUR.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor>(ca1,3));
    CHUR.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor>(ca1,3));

    SUC.push_back(tmv::SymBandMatrix<T,tmv::Upper,tmv::ColMajor>(a1,3));
    HUC.push_back(tmv::HermBandMatrix<T,tmv::Upper,tmv::ColMajor>(a1,3));
    CSUC.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor>(ca1,3));
    CHUC.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor>(ca1,3));

    SUD.push_back(tmv::SymBandMatrix<T,tmv::Upper,tmv::DiagMajor>(a1,3));
    HUD.push_back(tmv::HermBandMatrix<T,tmv::Upper,tmv::DiagMajor>(a1,3));
    CSUD.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(ca1,3));
    CHUD.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(ca1,3));

    SUD.push_back(tmv::SymTriDiagMatrix<tmv::Upper>(v1,v2));
    HUD.push_back(tmv::HermTriDiagMatrix<tmv::Upper>(v1,v2));
    CSUD.push_back(tmv::SymTriDiagMatrix<tmv::Upper>(cv1,cv2));
    CHUD.push_back(tmv::HermTriDiagMatrix<tmv::Upper>(cv1r,cv2));

    SUD.push_back(tmv::SymBandMatrix<T,tmv::Upper,tmv::DiagMajor>(SUD[1],0));
    HUD.push_back(tmv::HermBandMatrix<T,tmv::Upper,tmv::DiagMajor>(HUD[1],0));
    CSUD.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(CSUD[1],0));
    CHUD.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(CHUD[1],0));

#ifdef XTEST
    SLC.push_back(tmv::SymBandMatrix<T,tmv::Lower,tmv::ColMajor>(a1,7));
    HLC.push_back(tmv::HermBandMatrix<T,tmv::Lower,tmv::ColMajor>(a1,7));
    CSLC.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(ca1,7));
    CHLC.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(ca1,7));

    SLR.push_back(tmv::SymBandMatrix<T,tmv::Lower,tmv::RowMajor>(a1,7));
    HLR.push_back(tmv::HermBandMatrix<T,tmv::Lower,tmv::RowMajor>(a1,7));
    CSLR.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor>(ca1,7));
    CHLR.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor>(ca1,7));

    SLD.push_back(tmv::SymBandMatrix<T,tmv::Lower,tmv::DiagMajor>(a1,7));
    HLD.push_back(tmv::HermBandMatrix<T,tmv::Lower,tmv::DiagMajor>(a1,7));
    CSLD.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor>(ca1,7));
    CHLD.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor>(ca1,7));

    SLD.push_back(tmv::SymTriDiagMatrix<tmv::Lower>(v1,v2));
    HLD.push_back(tmv::HermTriDiagMatrix<tmv::Lower>(v1,v2));
    CSLD.push_back(tmv::SymTriDiagMatrix<tmv::Lower>(cv1,cv2));
    CHLD.push_back(tmv::HermTriDiagMatrix<tmv::Lower>(cv1r,cv2));

    SLC.push_back(tmv::SymBandMatrix<T,tmv::Lower,tmv::ColMajor>(SUD[1],0));
    HLC.push_back(tmv::HermBandMatrix<T,tmv::Lower,tmv::ColMajor>(HUD[1],0));
    CSLC.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(CSUD[1],0));
    CHLC.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(CHUD[1],0));

    a2a = a2;
    a2b = a2;
    ca2a = ca2;
    ca2b = ca2;
    ca2b.diag().Imag().Zero();
#endif
  }

  s.push_back(SUR[0].View());
  s.push_back(HUR[0].View());
  cs.push_back(CSUR[0].View());
  cs.push_back(CHUR[0].View());
  s.push_back(SUC[0].View());
  s.push_back(HUC[0].View());
  cs.push_back(CSUC[0].View());
  cs.push_back(CHUC[0].View());
  s.push_back(SUD[0].View());
  s.push_back(HUD[0].View());
  cs.push_back(CSUD[0].View());
  cs.push_back(CHUD[0].View());
  s.push_back(SUD[1].View());
  s.push_back(HUD[1].View());
  cs.push_back(CSUD[1].View());
  cs.push_back(CHUD[1].View());
  s.push_back(SUD[2].View());
  s.push_back(HUD[2].View());
  cs.push_back(CSUD[2].View());
  cs.push_back(CHUD[2].View());
#ifdef XTEST
  s.push_back(SLC[0].View());
  s.push_back(HLC[0].View());
  cs.push_back(CSLC[0].View());
  cs.push_back(CHLC[0].View());
  s.push_back(SLR[0].View());
  s.push_back(HLR[0].View());
  cs.push_back(CSLR[0].View());
  cs.push_back(CHLR[0].View());
  s.push_back(SLD[0].View());
  s.push_back(HLD[0].View());
  cs.push_back(CSLD[0].View());
  cs.push_back(CHLD[0].View());
  s.push_back(SLD[1].View());
  s.push_back(HLD[1].View());
  cs.push_back(CSLD[1].View());
  cs.push_back(CHLD[1].View());
  s.push_back(SLC[1].View());
  s.push_back(HLC[1].View());
  cs.push_back(CSLC[1].View());
  cs.push_back(CHLC[1].View());
  s.push_back(SymBandMatrixViewOf(a2a,tmv::Upper,6).SubSymBandMatrix(0,2*N,3,2));
  s.push_back(HermBandMatrixViewOf(a2b,tmv::Upper,6).SubSymBandMatrix(0,2*N,3,2));
  cs.push_back(SymBandMatrixViewOf(ca2a,tmv::Upper,6).SubSymBandMatrix(0,2*N,3,2));
  cs.push_back(HermBandMatrixViewOf(ca2b,tmv::Upper,6).SubSymBandMatrix(0,2*N,3,2));
#endif
}

template <class T> inline void MakeSymBandList_Sing(
  std::vector<tmv::SymBandMatrixView<T> >& s,
  std::vector<tmv::SymBandMatrixView<std::complex<T> > >& cs)
{
  static std::vector<tmv::SymBandMatrix<T,tmv::Upper,tmv::RowMajor> > SUR;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> > CSUR;
  static std::vector<tmv::SymBandMatrix<T,tmv::Upper,tmv::ColMajor> > SUC;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> > CSUC;
  static std::vector<tmv::SymBandMatrix<T,tmv::Upper,tmv::DiagMajor> > SUD;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor> > CSUD;
  static std::vector<tmv::HermBandMatrix<T,tmv::Upper,tmv::RowMajor> > HUR;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> > CHUR;
  static std::vector<tmv::HermBandMatrix<T,tmv::Upper,tmv::ColMajor> > HUC;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> > CHUC;
  static std::vector<tmv::HermBandMatrix<T,tmv::Upper,tmv::DiagMajor> > HUD;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor> > CHUD;
  static std::vector<tmv::SymBandMatrix<T,tmv::Lower,tmv::RowMajor> > SLR;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> > CSLR;
  static std::vector<tmv::SymBandMatrix<T,tmv::Lower,tmv::ColMajor> > SLC;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> > CSLC;
  static std::vector<tmv::SymBandMatrix<T,tmv::Lower,tmv::DiagMajor> > SLD;
  static std::vector<tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor> > CSLD;
  static std::vector<tmv::HermBandMatrix<T,tmv::Lower,tmv::RowMajor> > HLR;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> > CHLR;
  static std::vector<tmv::HermBandMatrix<T,tmv::Lower,tmv::ColMajor> > HLC;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> > CHLC;
  static std::vector<tmv::HermBandMatrix<T,tmv::Lower,tmv::DiagMajor> > HLD;
  static std::vector<tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor> > CHLD;
  const int N=10;
  static tmv::Matrix<T> a2a(2*N,2*N);
  static tmv::Matrix<T> a2b(2*N,2*N);
  static tmv::Matrix<std::complex<T> > ca2a(2*N,2*N);
  static tmv::Matrix<std::complex<T> > ca2b(2*N,2*N);

  if (SUR.size() == 0) {

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-5*j;
    tmv::Matrix<T> a2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-3*i+6*j;
    tmv::Matrix<std::complex<T> > ca1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) =
      std::complex<T>(3.+i-5.*j,2.-3.*i);
    tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
    for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) ca2(i,j) =
      std::complex<T>(1.-3*i+6.*j,-4.+2.*j);

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

    tmv::Vector<T> v1(N);
    tmv::Vector<T> v2(N-1);
    for (int i=0; i<N; ++i) v1(i) = 16.-3*i;
    for (int i=0; i<N-1; ++i) v2(i) = -6.+i;
    v1.SubVector(N/2,3*N/4).Zero();
    v2.SubVector(N/2,3*N/4).Zero();

    tmv::Vector<std::complex<T> > cv1(N);
    tmv::Vector<std::complex<T> > cv2(N-1);
    for (int i=0; i<N; ++i) cv1(i) = std::complex<T>(16.-3*i,i+4.);
    for (int i=0; i<N-1; ++i) cv2(i) = std::complex<T>(2*i-3.,-6.+i);
    cv1.SubVector(N/2,3*N/4).Zero();
    cv2.SubVector(N/2,3*N/4).Zero();
    tmv::Vector<std::complex<T> > cv1r = cv1;
    cv1r.Imag().Zero();

    SUR.push_back(tmv::SymBandMatrix<T,tmv::Upper,tmv::RowMajor>(a1,3));
    HUR.push_back(tmv::HermBandMatrix<T,tmv::Upper,tmv::RowMajor>(a1,3));
    CSUR.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor>(ca1,3));
    CHUR.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor>(ca1,3));

    SUC.push_back(tmv::SymBandMatrix<T,tmv::Upper,tmv::ColMajor>(a1,3));
    HUC.push_back(tmv::HermBandMatrix<T,tmv::Upper,tmv::ColMajor>(a1,3));
    CSUC.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor>(ca1,3));
    CHUC.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor>(ca1,3));

    SUD.push_back(tmv::SymBandMatrix<T,tmv::Upper,tmv::DiagMajor>(a1,3));
    HUD.push_back(tmv::HermBandMatrix<T,tmv::Upper,tmv::DiagMajor>(a1,3));
    CSUD.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(ca1,3));
    CHUD.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(ca1,3));

    SUD.push_back(tmv::SymTriDiagMatrix<tmv::Upper>(v1,v2));
    HUD.push_back(tmv::HermTriDiagMatrix<tmv::Upper>(v1,v2));
    CSUD.push_back(tmv::SymTriDiagMatrix<tmv::Upper>(cv1,cv2));
    CHUD.push_back(tmv::HermTriDiagMatrix<tmv::Upper>(cv1r,cv2));

    SUD.push_back(tmv::SymBandMatrix<T,tmv::Upper,tmv::DiagMajor>(a1,0));
    HUD.push_back(tmv::HermBandMatrix<T,tmv::Upper,tmv::DiagMajor>(a1,0));
    CSUD.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(ca1,0));
    CHUD.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Upper,tmv::DiagMajor>(ca1,0));

#ifdef XTEST
    SLC.push_back(tmv::SymBandMatrix<T,tmv::Lower,tmv::ColMajor>(a1,7));
    HLC.push_back(tmv::HermBandMatrix<T,tmv::Lower,tmv::ColMajor>(a1,7));
    CSLC.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(ca1,7));
    CHLC.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(ca1,7));

    SLR.push_back(tmv::SymBandMatrix<T,tmv::Lower,tmv::RowMajor>(a1,7));
    HLR.push_back(tmv::HermBandMatrix<T,tmv::Lower,tmv::RowMajor>(a1,7));
    CSLR.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor>(ca1,7));
    CHLR.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor>(ca1,7));

    SLD.push_back(tmv::SymBandMatrix<T,tmv::Lower,tmv::DiagMajor>(a1,7));
    HLD.push_back(tmv::HermBandMatrix<T,tmv::Lower,tmv::DiagMajor>(a1,7));
    CSLD.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor>(ca1,7));
    CHLD.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::DiagMajor>(ca1,7));

    SLD.push_back(tmv::SymTriDiagMatrix<tmv::Lower>(v1,v2));
    HLD.push_back(tmv::HermTriDiagMatrix<tmv::Lower>(v1,v2));
    CSLD.push_back(tmv::SymTriDiagMatrix<tmv::Lower>(cv1,cv2));
    CHLD.push_back(tmv::HermTriDiagMatrix<tmv::Lower>(cv1r,cv2));

    SLC.push_back(tmv::SymBandMatrix<T,tmv::Lower,tmv::ColMajor>(a1,0));
    HLC.push_back(tmv::HermBandMatrix<T,tmv::Lower,tmv::ColMajor>(a1,0));
    CSLC.push_back(tmv::SymBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(ca1,0));
    CHLC.push_back(tmv::HermBandMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor>(ca1,0));

    a2a = a2;
    a2b = a2;
    ca2a = ca2;
    ca2b = ca2;
    ca2b.diag().Imag().Zero();
#endif
  }

  s.push_back(SUR[0].View());
  s.push_back(HUR[0].View());
  cs.push_back(CSUR[0].View());
  cs.push_back(CHUR[0].View());
  s.push_back(SUC[0].View());
  s.push_back(HUC[0].View());
  cs.push_back(CSUC[0].View());
  cs.push_back(CHUC[0].View());
  s.push_back(SUD[0].View());
  s.push_back(HUD[0].View());
  cs.push_back(CSUD[0].View());
  cs.push_back(CHUD[0].View());
  s.push_back(SUD[1].View());
  s.push_back(HUD[1].View());
  cs.push_back(CSUD[1].View());
  cs.push_back(CHUD[1].View());
  s.push_back(SUD[2].View());
  s.push_back(HUD[2].View());
  cs.push_back(CSUD[2].View());
  cs.push_back(CHUD[2].View());
#ifdef XTEST
  s.push_back(SLC[0].View());
  s.push_back(HLC[0].View());
  cs.push_back(CSLC[0].View());
  cs.push_back(CHLC[0].View());
  s.push_back(SLR[0].View());
  s.push_back(HLR[0].View());
  cs.push_back(CSLR[0].View());
  cs.push_back(CHLR[0].View());
  s.push_back(SLD[0].View());
  s.push_back(HLD[0].View());
  cs.push_back(CSLD[0].View());
  cs.push_back(CHLD[0].View());
  s.push_back(SLD[1].View());
  s.push_back(HLD[1].View());
  cs.push_back(CSLD[1].View());
  cs.push_back(CHLD[1].View());
  s.push_back(SLC[1].View());
  s.push_back(HLC[1].View());
  cs.push_back(CSLC[1].View());
  cs.push_back(CHLC[1].View());
  s.push_back(SymBandMatrixViewOf(a2a,tmv::Upper,6).SubSymBandMatrix(0,2*N,3,2));
  s.push_back(HermBandMatrixViewOf(a2b,tmv::Upper,6).SubSymBandMatrix(0,2*N,3,2));
  cs.push_back(SymBandMatrixViewOf(ca2a,tmv::Upper,6).SubSymBandMatrix(0,2*N,3,2));
  cs.push_back(HermBandMatrixViewOf(ca2b,tmv::Upper,6).SubSymBandMatrix(0,2*N,3,2));
#endif
}

template <class T> inline void MakeSymBandList(
  std::vector<tmv::SymBandMatrixView<T> >& s,
  std::vector<tmv::SymBandMatrixView<std::complex<T> > >& cs, PosDefCode pdc)
{
  if (pdc == PosDef) {
    MakeSymBandList_PosDef(s,cs);
  } else if (pdc == InDef) {
    MakeSymBandList_InDef(s,cs);
  } else {
    MakeSymBandList_Sing(s,cs);
  }
}
