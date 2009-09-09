///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


template <class T> inline void MakeSymList_PosDef(
  std::vector<tmv::SymMatrixView<T> >& s,
  std::vector<tmv::SymMatrixView<std::complex<T> > >& cs)
{
  const int N=10;

  static tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-5*j;
  static tmv::Matrix<std::complex<T> > ca1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) = 
    std::complex<T>(3.+i-5*j,2.-3.*i);
  static tmv::Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-3*i+6*j;
  static tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) ca2(i,j) = 
    std::complex<T>(1.-3.*i+6*j,-4.+2.*j);

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

  static tmv::SymMatrix<T,tmv::Upper,tmv::RowMajor> S1(a1);
  static tmv::HermMatrix<T,tmv::Upper,tmv::RowMajor> H1(a1);
  static tmv::SymMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> CS1(ca1);
  static tmv::HermMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> CH1(ca1);
  s.push_back(S1.View());
  s.push_back(H1.View());
  cs.push_back(CS1.View());
  cs.push_back(CH1.View());

  static tmv::SymMatrix<T,tmv::Upper,tmv::ColMajor> S2(a1);
  static tmv::HermMatrix<T,tmv::Upper,tmv::ColMajor> H2(a1);
  static tmv::SymMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> CS2(ca1);
  static tmv::HermMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> CH2(ca1);
  s.push_back(S2.View());
  s.push_back(H2.View());
  cs.push_back(CS2.View());
  cs.push_back(CH2.View());

#ifdef XTEST
  static tmv::SymMatrix<T,tmv::Lower,tmv::RowMajor> S3(a1);
  static tmv::HermMatrix<T,tmv::Lower,tmv::RowMajor> H3(a1);
  static tmv::SymMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> CS3(ca1);
  static tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> CH3(ca1);
  static tmv::SymMatrix<T,tmv::Lower,tmv::ColMajor> S4(a1);
  static tmv::HermMatrix<T,tmv::Lower,tmv::ColMajor> H4(a1);
  static tmv::SymMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> CS4(ca1);
  static tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> CH4(ca1);
  static tmv::Matrix<T> a2a = a2;
  static tmv::Matrix<T> a2b = a2;
  static tmv::Matrix<std::complex<T> > ca2a = ca2;
  static tmv::Matrix<std::complex<T> > ca2b = ca2;
  ca2b.diag().Imag().Zero();
  s.push_back(S3.View());
  s.push_back(H3.View());
  cs.push_back(CS3.View());
  cs.push_back(CH3.View());
  s.push_back(S4.View());
  s.push_back(H4.View());
  cs.push_back(CS4.View());
  cs.push_back(CH4.View());
  s.push_back(SymMatrixViewOf(a2a,tmv::Upper).SubSymMatrix(0,2*N,2));
  s.push_back(HermMatrixViewOf(a2b,tmv::Upper).SubSymMatrix(0,2*N,2));
  cs.push_back(SymMatrixViewOf(ca2a,tmv::Upper).SubSymMatrix(0,2*N,2));
  cs.push_back(HermMatrixViewOf(ca2b,tmv::Upper).SubSymMatrix(0,2*N,2));
#endif
}

template <class T> inline void MakeSymList_InDef(
  std::vector<tmv::SymMatrixView<T> >& s,
  std::vector<tmv::SymMatrixView<std::complex<T> > >& cs)
{
  const int N=10;

  static tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-5*j;
  static tmv::Matrix<std::complex<T> > ca1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) = 
    std::complex<T>(3.+i-5*j,2.-3.*i);
  static tmv::Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-3*i+6*j;
  static tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) ca2(i,j) = 
    std::complex<T>(1.-3.*i+6*j,-4.+2.*j);

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

  static tmv::SymMatrix<T,tmv::Upper,tmv::RowMajor> S1(a1);
  static tmv::HermMatrix<T,tmv::Upper,tmv::RowMajor> H1(a1);
  static tmv::SymMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> CS1(ca1);
  static tmv::HermMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> CH1(ca1);
  s.push_back(S1.View());
  s.push_back(H1.View());
  cs.push_back(CS1.View());
  cs.push_back(CH1.View());

  static tmv::SymMatrix<T,tmv::Upper,tmv::ColMajor> S2(a1);
  static tmv::HermMatrix<T,tmv::Upper,tmv::ColMajor> H2(a1);
  static tmv::SymMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> CS2(ca1);
  static tmv::HermMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> CH2(ca1);
  s.push_back(S2.View());
  s.push_back(H2.View());
  cs.push_back(CS2.View());
  cs.push_back(CH2.View());

#ifdef XTEST
  static tmv::SymMatrix<T,tmv::Lower,tmv::RowMajor> S3(a1);
  static tmv::HermMatrix<T,tmv::Lower,tmv::RowMajor> H3(a1);
  static tmv::SymMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> CS3(ca1);
  static tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> CH3(ca1);
  static tmv::SymMatrix<T,tmv::Lower,tmv::ColMajor> S4(a1);
  static tmv::HermMatrix<T,tmv::Lower,tmv::ColMajor> H4(a1);
  static tmv::SymMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> CS4(ca1);
  static tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> CH4(ca1);
  static tmv::Matrix<T> a2a = a2;
  static tmv::Matrix<T> a2b = a2;
  static tmv::Matrix<std::complex<T> > ca2a = ca2;
  static tmv::Matrix<std::complex<T> > ca2b = ca2;
  ca2b.diag().Imag().Zero();
  s.push_back(S3.View());
  s.push_back(H3.View());
  cs.push_back(CS3.View());
  cs.push_back(CH3.View());
  s.push_back(S4.View());
  s.push_back(H4.View());
  cs.push_back(CS4.View());
  cs.push_back(CH4.View());
  s.push_back(SymMatrixViewOf(a2a,tmv::Upper).SubSymMatrix(0,2*N,2));
  s.push_back(HermMatrixViewOf(a2b,tmv::Upper).SubSymMatrix(0,2*N,2));
  cs.push_back(SymMatrixViewOf(ca2a,tmv::Upper).SubSymMatrix(0,2*N,2));
  cs.push_back(HermMatrixViewOf(ca2b,tmv::Upper).SubSymMatrix(0,2*N,2));
#endif
}

template <class T> inline void MakeSymList_Sing(
  std::vector<tmv::SymMatrixView<T> >& s,
  std::vector<tmv::SymMatrixView<std::complex<T> > >& cs)
{
  const int N=10;

  static tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-5*j;
  static tmv::Matrix<std::complex<T> > ca1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) ca1(i,j) = 
    std::complex<T>(3.+i-5*j,2.-3.*i);
  static tmv::Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-3*i+6*j;
  static tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) ca2(i,j) = 
    std::complex<T>(1.-3.*i+6*j,-4.+2.*j);

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

  static tmv::SymMatrix<T,tmv::Upper,tmv::RowMajor> S1(a1);
  static tmv::HermMatrix<T,tmv::Upper,tmv::RowMajor> H1(a1);
  static tmv::SymMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> CS1(ca1);
  static tmv::HermMatrix<std::complex<T>,tmv::Upper,tmv::RowMajor> CH1(ca1);
  s.push_back(S1.View());
  s.push_back(H1.View());
  cs.push_back(CS1.View());
  cs.push_back(CH1.View());

  static tmv::SymMatrix<T,tmv::Upper,tmv::ColMajor> S2(a1);
  static tmv::HermMatrix<T,tmv::Upper,tmv::ColMajor> H2(a1);
  static tmv::SymMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> CS2(ca1);
  static tmv::HermMatrix<std::complex<T>,tmv::Upper,tmv::ColMajor> CH2(ca1);
  s.push_back(S2.View());
  s.push_back(H2.View());
  cs.push_back(CS2.View());
  cs.push_back(CH2.View());

#ifdef XTEST
  static tmv::SymMatrix<T,tmv::Lower,tmv::RowMajor> S3(a1);
  static tmv::HermMatrix<T,tmv::Lower,tmv::RowMajor> H3(a1);
  static tmv::SymMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> CS3(ca1);
  static tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::RowMajor> CH3(ca1);
  static tmv::SymMatrix<T,tmv::Lower,tmv::ColMajor> S4(a1);
  static tmv::HermMatrix<T,tmv::Lower,tmv::ColMajor> H4(a1);
  static tmv::SymMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> CS4(ca1);
  static tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> CH4(ca1);
  static tmv::Matrix<T> a2a = a2;
  static tmv::Matrix<T> a2b = a2;
  static tmv::Matrix<std::complex<T> > ca2a = ca2;
  static tmv::Matrix<std::complex<T> > ca2b = ca2;
  ca2b.diag().Imag().Zero();
  s.push_back(S3.View());
  s.push_back(H3.View());
  cs.push_back(CS3.View());
  cs.push_back(CH3.View());
  s.push_back(S4.View());
  s.push_back(H4.View());
  cs.push_back(CS4.View());
  cs.push_back(CH4.View());
  s.push_back(SymMatrixViewOf(a2a,tmv::Upper).SubSymMatrix(0,2*N,2));
  s.push_back(HermMatrixViewOf(a2b,tmv::Upper).SubSymMatrix(0,2*N,2));
  cs.push_back(SymMatrixViewOf(ca2a,tmv::Upper).SubSymMatrix(0,2*N,2));
  cs.push_back(HermMatrixViewOf(ca2b,tmv::Upper).SubSymMatrix(0,2*N,2));
#endif
}

template <class T> inline void MakeSymList(
  std::vector<tmv::SymMatrixView<T> >& s,
  std::vector<tmv::SymMatrixView<std::complex<T> > >& cs, PosDefCode pdc)
{
  if (pdc == PosDef) {
    MakeSymList_PosDef(s,cs);
  } else if (pdc == InDef) {
    MakeSymList_InDef(s,cs);
  } else {
    MakeSymList_Sing(s,cs);
  }
}
