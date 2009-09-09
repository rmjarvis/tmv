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


template <class T> inline void MakeBandList(
    std::vector<tmv::BandMatrixView<T> >& b,
    std::vector<tmv::BandMatrixView<std::complex<T> > >& cb)
{
  const int N = 10;

  static tmv::Matrix<T> a1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = 3.+i-5*j;
  a1.diag().AddToAll(T(3*N));
  static tmv::Matrix<T> a2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) a2(i,j) = 1.-3*i+6*j;
  a2.diag().AddToAll(T(3*N));
  static tmv::Vector<T> v1(N);
  static tmv::Vector<T> v2(N-1);
  for (int i=0; i<N; ++i) v1(i) = 16.-3*i; 
  for (int i=0; i<N-1; ++i) v2(i) = -6.+i; 

  static tmv::Matrix<std::complex<T> > ca1(N,N);
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
    ca1(i,j) = std::complex<T>(3.+i-5*j,4.-8*i-j);
  ca1.diag().AddToAll(T(3*N));
  static tmv::Matrix<std::complex<T> > ca2(2*N,2*N);
  for (int i=0; i<2*N; ++i) for (int j=0; j<2*N; ++j) 
    ca2(i,j) = std::complex<T>(1.-3*i+6*j,8.+2*i-6*j);
  ca2.diag().AddToAll(T(3*N));
  static tmv::Vector<std::complex<T> > cv1(N);
  static tmv::Vector<std::complex<T> > cv2(N-1);
  for (int i=0; i<N; ++i) cv1(i) = std::complex<T>(16.-3*i,i+4.); 
  for (int i=0; i<N-1; ++i) cv2(i) = std::complex<T>(2*i-3.,-6.+i); 

  static tmv::BandMatrix<T,tmv::RowMajor> B1(a1,3,1);
  static tmv::BandMatrix<std::complex<T>,tmv::RowMajor> C1(ca1,3,1);
  b.push_back(B1.View());
  cb.push_back(C1.View());
  static tmv::BandMatrix<T,tmv::DiagMajor> B3(a1,3,1);
  static tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> C3(ca1,3,1);
  b.push_back(B3.View());
  cb.push_back(C3.View());
  static tmv::BandMatrix<T,tmv::DiagMajor> B5(B1,1,1);
  static tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> C5(C1,1,1);
  b.push_back(B5.View());
  cb.push_back(C5.View());
  static tmv::BandMatrix<T> B6(B1,3,0);
  static tmv::BandMatrix<std::complex<T> > C6(C1,3,0);
  b.push_back(B6.View());
  cb.push_back(C6.View());
#ifdef XTEST
  static tmv::BandMatrix<T> B4(a2,6,6);
  static tmv::BandMatrix<std::complex<T> > C4(ca2,6,6);
  b.push_back(B4.SubBandMatrix(0,2*N,0,N,3,3));
  cb.push_back(C4.SubBandMatrix(0,2*N,0,N,3,3));
  static tmv::BandMatrix<T> B4a(B4);
  static tmv::BandMatrix<std::complex<T> > C4a(C4);
  b.push_back(B4a.SubBandMatrix(0,N+2,0,N,4,4));
  cb.push_back(C4a.SubBandMatrix(0,N+2,0,N,4,4));
  static tmv::BandMatrix<T> B4b(B4);
  static tmv::BandMatrix<std::complex<T> > C4b(C4);
  b.push_back(B4b.SubBandMatrix(0,2*N,0,2*N,3,3,2,2));
  cb.push_back(C4b.SubBandMatrix(0,2*N,0,2*N,3,3,2,2));
  static tmv::BandMatrix<T,tmv::ColMajor> B2(a1,3,1);
  static tmv::BandMatrix<std::complex<T>,tmv::ColMajor> C2(ca1,3,1);
  b.push_back(B2.View());
  cb.push_back(C2.View());
  static tmv::BandMatrix<T,tmv::DiagMajor> B3b(a1,1,3);
  static tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> C3b(ca1,1,3);
  b.push_back(B3b.View());
  cb.push_back(C3b.View());
  static tmv::BandMatrix<T> B7(a1,0,3);
  static tmv::BandMatrix<std::complex<T> > C7(ca1,0,3);
  b.push_back(B7.View());
  cb.push_back(C7.View());
  static tmv::BandMatrix<T,tmv::DiagMajor> B8 = tmv::UpperBiDiagMatrix(v1,v2);
  static tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> C8 = 
    tmv::UpperBiDiagMatrix(cv1,cv2);
  b.push_back(B8.View());
  cb.push_back(C8.View());
  static tmv::BandMatrix<T,tmv::DiagMajor> B9 = tmv::LowerBiDiagMatrix(v2,v1);
  static tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> C9 =
    tmv::LowerBiDiagMatrix(cv2,cv1);
  b.push_back(B9.View());
  cb.push_back(C9.View());
  static tmv::BandMatrix<T,tmv::DiagMajor> B10 = tmv::TriDiagMatrix(v2,v1,v2);
  static tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> C10 =
    tmv::TriDiagMatrix(cv2,cv1,cv2);
  b.push_back(B10.View());
  cb.push_back(C10.View());
  static tmv::BandMatrix<T,tmv::DiagMajor> B11 = tmv::UpperBiDiagMatrix(v1,v1);
  static tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> C11 =
    tmv::UpperBiDiagMatrix(cv1,cv1);
  b.push_back(B11.View());
  cb.push_back(C11.View());
  static tmv::BandMatrix<T,tmv::DiagMajor> B12 = tmv::LowerBiDiagMatrix(v1,v1);
  static tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> C12 =
    tmv::LowerBiDiagMatrix(cv1,cv1);
  b.push_back(B12.View());
  cb.push_back(C12.View());
  static tmv::BandMatrix<T,tmv::DiagMajor> B13 = tmv::TriDiagMatrix(v1,v1,v2);
  static tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> C13 =
    tmv::TriDiagMatrix(cv1,cv1,cv2);
  b.push_back(B13.View());
  cb.push_back(C13.View());
  static tmv::BandMatrix<T,tmv::DiagMajor> B14 = tmv::TriDiagMatrix(v2,v1,v1);
  static tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> C14 =
    tmv::TriDiagMatrix(cv2,cv1,cv1);
  b.push_back(B14.View());
  cb.push_back(C14.View());
  static tmv::BandMatrix<T> B4c(B4);
  static tmv::BandMatrix<std::complex<T> > C4c(C4);
  b.push_back(B4c.SubBandMatrix(0,N,0,2*N,3,3));
  cb.push_back(C4c.SubBandMatrix(0,N,0,2*N,3,3));
  static tmv::BandMatrix<T> B4d(B4);
  static tmv::BandMatrix<std::complex<T> > C4d(C4);
  b.push_back(B4d.SubBandMatrix(0,N,0,N+2,4,4));
  cb.push_back(C4d.SubBandMatrix(0,N,0,N+2,4,4));
  b.push_back(BandMatrixViewOf(a1,3,N-1));
  cb.push_back(BandMatrixViewOf(ca1,3,N-1));
  static tmv::BandMatrix<T,tmv::DiagMajor> B15(a1,1,N-1);
  static tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> C15(ca1,1,N-1);
  b.push_back(B15.View());
  cb.push_back(C15.View());
  static tmv::BandMatrix<T,tmv::DiagMajor> B16(a1,3,N-2);
  static tmv::BandMatrix<std::complex<T>,tmv::DiagMajor> C16(ca1,3,N-2);
  b.push_back(B16.View());
  cb.push_back(C16.View());
  static tmv::BandMatrix<T> B17(a1,0,0);
  static tmv::BandMatrix<std::complex<T> > C17(ca1,0,0);
  b.push_back(B17.View());
  cb.push_back(C17.View());
#endif
}
