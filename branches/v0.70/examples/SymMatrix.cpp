#include "TMV.h"
// Note: extra include file for SymMatrix
#include "TMV_Sym.h"
// Also need to link with -ltmv_symband
#include <iostream>

int main() try 
{
  // Several ways to create and initialize matrices:
  
  // Create with uninitialized values
  tmv::SymMatrix<double> m1(5); 
  // In debug mode, the elements are all set to 888 to make it easier
  // to find errors related to not correctly initializing the matrix.
  std::cout<<"m1 = \n"<<m1<<std::endl;
  //!  m1 = 
  //!  5 5
  //!  (  888  888  888  888  888  )
  //!  (  888  888  888  888  888  )
  //!  (  888  888  888  888  888  )
  //!  (  888  888  888  888  888  )
  //!  (  888  888  888  888  888  )

  // Initialize with element access:
  for(size_t i=0;i<m1.nrows();i++) 
    for(size_t j=0;j<m1.ncols();j++) 
      if (i>=j) m1(i,j) = 2.*i-j*j+2.; 
  std::cout<<"m1 = \n"<<m1<<std::endl;
  //!  m1 = 
  //!  5 5
  //!  (  2  4  6  8  10  )
  //!  (  4  3  5  7  9  )
  //!  (  6  5  2  4  6  )
  //!  (  8  7  4  -1  1  )
  //!  (  10  9  6  1  -6  )

  tmv::SymMatrix<double> m2(5,2.); // Create with all 2's.
  std::cout<<"m2 = \n"<<m2<<std::endl;
  //!  m2 = 
  //!  5 5
  //!  (  2  2  2  2  2  )
  //!  (  2  2  2  2  2  )
  //!  (  2  2  2  2  2  )
  //!  (  2  2  2  2  2  )
  //!  (  2  2  2  2  2  )

  // Create from given elements:
  double mm[25] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
    16,17,18,19,20,21,22,23,24,25};
  tmv::SymMatrix<double> m3(5); // Default storage is Upper, ColMajor

  std::cout<<"m3 (Upper, ColMajor) = \n"<<m3<<std::endl;
  //!  m3 (Upper, ColMajor) = 
  //!  5 5
  //!  (  1  6  11  16  21  )
  //!  (  6  7  12  17  22  )
  //!  (  11  12  13  18  23  )
  //!  (  16  17  18  19  24  )
  //!  (  21  22  23  24  25  )
  tmv::SymMatrix<double,tmv::Upper,tmv::RowMajor> m4(5,mm);
  std::cout<<"m4 (Upper,RowMajor) = \n"<<m4<<std::endl;
  //!  m4 (Upper,RowMajor) = 
  //!  5 5
  //!  (  1  2  3  4  5  )
  //!  (  2  7  8  9  10  )
  //!  (  3  8  13  14  15  )
  //!  (  4  9  14  19  20  )
  //!  (  5  10  15  20  25  )
  tmv::SymMatrix<double,tmv::Lower,tmv::ColMajor> m5(5,mm);
  std::cout<<"m5 (Lower,ColMajor) = \n"<<m5<<std::endl;
  //!  m5 (Lower,ColMajor) = 
  //!  5 5
  //!  (  1  2  3  4  5  )
  //!  (  2  7  8  9  10  )
  //!  (  3  8  13  14  15  )
  //!  (  4  9  14  19  20  )
  //!  (  5  10  15  20  25  )
  tmv::SymMatrix<double,tmv::Lower,tmv::RowMajor> m6(5,mm);
  std::cout<<"m6 (Lower,RowMajor) = \n"<<m6<<std::endl;
  //!  m6 (Lower,RowMajor) = 
  //!  5 5
  //!  (  1  6  11  16  21  )
  //!  (  6  7  12  17  22  )
  //!  (  11  12  13  18  23  )
  //!  (  16  17  18  19  24  )
  //!  (  21  22  23  24  25  )

  // Make from the corresponding portion of a regular Matrix
  // Note - it copies from the correct triangle regardless of 
  // the storage order of the two matrices.
  tmv::Matrix<double,tmv::RowMajor> xm(5,5);
  xm <<
     2, 5, 1, 9, 8,
     1, 5, 7, 2, 0,
     3, 9, 6, 8, 4,
     4, 2, 1, 9, 0,
     9, 8, 3, 7, 5;
  tmv::SymMatrix<double,tmv::Lower,tmv::ColMajor> m7(xm);
  std::cout<<"m7 (Lower,ColMajor) = "<<m7<<std::endl;
  //!  m7 (Lower,ColMajor) = 5 5
  //!  (  2  1  3  4  9  )
  //!  (  1  5  9  2  8  )
  //!  (  3  9  6  1  3  )
  //!  (  4  2  1  9  7  )
  //!  (  9  8  3  7  5  )


  // Norms, etc. 

  std::cout<<"Norm1(m1) = "<<Norm1(m1)<<std::endl;
  //!  Norm1(m1) = 32
  std::cout<<"Norm2(m1) = "<<Norm2(m1)<<std::endl;
  //!  Norm2(m1) = 24.4968
  std::cout<<"NormInf(m1) = "<<NormInf(m1)<<std::endl;
  //!  NormInf(m1) = 32
  std::cout<<"NormF(m1) = "<<NormF(m1);
  std::cout<<" = "<<Norm(m1)<<std::endl;
  //!  NormF(m1) = 30.0333 = 30.0333
  std::cout<<"MaxAbsElement(m1) = "<<MaxAbsElement(m1)<<std::endl;
  //!  MaxAbsElement(m1) = 10
  std::cout<<"Trace(m1) = "<<Trace(m1)<<std::endl;
  //!  Trace(m1) = 0
  std::cout<<"Det(m1) = "<<Det(m1)<<std::endl;
  //!  Det(m1) = 4866


  // Views:

  std::cout<<"m1 = \n"<<m1<<std::endl;
  //!  m1 = 
  //!  5 5
  //!  (  2  4  6  8  10  )
  //!  (  4  3  5  7  9  )
  //!  (  6  5  2  4  6  )
  //!  (  8  7  4  -1  1  )
  //!  (  10  9  6  1  -6  )
  std::cout<<"m1.diag() = "<<m1.diag()<<std::endl;
  //!  m1.diag() = 5 ( 2  3  2  -1  -6 )
  std::cout<<"m1.diag(1) = "<<m1.diag(1)<<std::endl;
  //!  m1.diag(1) = 4 ( 4  5  4  1 )
  std::cout<<"m1.diag(-1) = "<<m1.diag(-1)<<std::endl;
  //!  m1.diag(-1) = 4 ( 4  5  4  1 )
  std::cout<<"m1.subSymMatrix(0,3) = \n"<<m1.subSymMatrix(0,3)<<std::endl;
  //!  m1.subSymMatrix(0,3) = 
  //!  3 3
  //!  (  2  4  6  )
  //!  (  4  3  5  )
  //!  (  6  5  2  )
  std::cout<<"m1.upperTri() = \n"<<m1.upperTri()<<std::endl;
  //!  m1.upperTri() = 
  //!  5 5
  //!  (  2  4  6  8  10  )
  //!  (  0  3  5  7  9  )
  //!  (  0  0  2  4  6  )
  //!  (  0  0  0  -1  1  )
  //!  (  0  0  0  0  -6  )
  std::cout<<"m1.lowerTri(tmv::UnitDiag) = \n"<<
    m1.lowerTri(tmv::UnitDiag)<<std::endl;
  //!  m1.lowerTri(tmv::UnitDiag) = 
  //!  5 5
  //!  (  1  0  0  0  0  )
  //!  (  4  1  0  0  0  )
  //!  (  6  5  1  0  0  )
  //!  (  8  7  4  1  0  )
  //!  (  10  9  6  1  1  )


  // Fortran-style indexing:

  tmv::SymMatrix<double,tmv::Upper,tmv::ColMajor,tmv::FortranStyle> fm1 = m1;
  std::cout<<"fm1 = m1 = \n"<<fm1<<std::endl;
  //!  fm1 = m1 = 
  //!  5 5
  //!  (  2  4  6  8  10  )
  //!  (  4  3  5  7  9  )
  //!  (  6  5  2  4  6  )
  //!  (  8  7  4  -1  1  )
  //!  (  10  9  6  1  -6  )
  std::cout<<"fm1(1,1) = "<<fm1(1,1)<<std::endl;
  //!  fm1(1,1) = 2
  std::cout<<"fm1(4,3) = "<<fm1(4,3)<<std::endl;
  //!  fm1(4,3) = 4
  std::cout<<"fm1.subSymMatrix(1,3) = \n"<<fm1.subSymMatrix(1,3)<<std::endl;
  //!  fm1.subSymMatrix(1,3) = 
  //!  3 3
  //!  (  2  4  6  )
  //!  (  4  3  5  )
  //!  (  6  5  2  )


  // Matrix arithmetic:

  tmv::SymMatrix<double> m1pm3 = m1 + m3;  
  std::cout<<"m1 + m3 = \n"<<m1pm3<<std::endl;
  //!  m1 + m3 = 
  //!  5 5
  //!  (  3  10  17  24  31  )
  //!  (  10  10  17  24  31  )
  //!  (  17  17  15  22  29  )
  //!  (  24  24  22  18  25  )
  //!  (  31  31  29  25  19  )
  // Works correctly even if matrices are stored in different order:
  tmv::SymMatrix<double> m3pm4 = m3 + m4; 
  std::cout<<"m3 + m4 = \n"<<m3pm4<<std::endl;
  //!  m3 + m4 = 
  //!  5 5
  //!  (  2  8  14  20  26  )
  //!  (  8  14  20  26  32  )
  //!  (  14  20  26  32  38  )
  //!  (  20  26  32  38  44  )
  //!  (  26  32  38  44  50  )
  tmv::SymMatrix<double> m4pm5 = m4 + m5; 
  std::cout<<"m4 + m5 = \n"<<m4pm5<<std::endl;
  //!  m4 + m5 = 
  //!  5 5
  //!  (  2  4  6  8  10  )
  //!  (  4  14  16  18  20  )
  //!  (  6  16  26  28  30  )
  //!  (  8  18  28  38  40  )
  //!  (  10  20  30  40  50  )

  m1 *= 2.;
  std::cout<<"m1 *= 2 = \n"<<m1<<std::endl;
  //!  m1 *= 2 = 
  //!  5 5
  //!  (  4  8  12  16  20  )
  //!  (  8  6  10  14  18  )
  //!  (  12  10  4  8  12  )
  //!  (  16  14  8  -2  2  )
  //!  (  20  18  12  2  -12  )

  m1 += m4;
  std::cout<<"m1 += m4 = \n"<<m1<<std::endl;
  //!  m1 += m4 = 
  //!  5 5
  //!  (  5  10  15  20  25  )
  //!  (  10  13  18  23  28  )
  //!  (  15  18  17  22  27  )
  //!  (  20  23  22  17  22  )
  //!  (  25  28  27  22  13  )

  // Vector outer product
  tmv::Vector<double> v = xm.col(0);
  tmv::SymMatrix<double> vv = v^v;
  std::cout<<"v = "<<v<<std::endl;
  //!  v = 5 ( 2  1  3  4  9 )
  std::cout<<"v^v = \n"<<vv<<std::endl;
  //!  v^v = 
  //!  5 5
  //!  (  4  2  6  8  18  )
  //!  (  2  1  3  4  9  )
  //!  (  6  3  9  12  27  )
  //!  (  8  4  12  16  36  )
  //!  (  18  9  27  36  81  )

  // SymMatrix * Vector product
  std::cout<<"m1 * v = "<<m1*v<<std::endl;
  //!  m1 * v = 5 ( 370  431  430  395  364 )
  std::cout<<"v * m1 = "<<v*m1<<std::endl;
  //!  v * m1 = 5 ( 370  431  430  395  364 )

  // SymMatrix * Matrix product
  tmv::Matrix<double> m1xm = m1 * xm; 
  std::cout<<"m1 * xm = \n"<<m1xm<<std::endl;
  //!  m1 * xm = 
  //!  5  5
  //!  (  370  450  260  540  225  )
  //!  (  431  547  316  663  292  )
  //!  (  430  578  346  694  323  )
  //!  (  395  623  396  709  358  )
  //!  (  364  656  444  786  373  )
  // Note: the product of two symmetrix matrices is not symmetric!:
  tmv::Matrix<double> m1m5 = m1 * m5; 
  std::cout<<"m1 * m5 = \n"<<m1m5<<std::endl;
  //!  m1 * m5 = 
  //!  5  5
  //!  (  275  630  945  1200  1375  )
  //!  (  322  742  1110  1406  1610  )
  //!  (  325  760  1123  1418  1625  )
  //!  (  310  750  1098  1358  1550  )
  //!  (  315  790  1153  1408  1575  )
  

  // Complex matrices:

  tmv::SymMatrix<std::complex<double> > cm1 = m1 * std::complex<double>(1,2);
  std::cout<<"cm1 = m1 * (1+2i) = \n"<<cm1<<std::endl;
  //!  cm1 = m1 * (1+2i) = 
  //!  5 5
  //!  (  (5,10)  (10,20)  (15,30)  (20,40)  (25,50)  )
  //!  (  (10,20)  (13,26)  (18,36)  (23,46)  (28,56)  )
  //!  (  (15,30)  (18,36)  (17,34)  (22,44)  (27,54)  )
  //!  (  (20,40)  (23,46)  (22,44)  (17,34)  (22,44)  )
  //!  (  (25,50)  (28,56)  (27,54)  (22,44)  (13,26)  )
  std::cout<<"cm1.conjugate() = \n"<<cm1.conjugate()<<std::endl;
  //!  cm1.conjugate() = 
  //!  5 5
  //!  (  (5,-10)  (10,-20)  (15,-30)  (20,-40)  (25,-50)  )
  //!  (  (10,-20)  (13,-26)  (18,-36)  (23,-46)  (28,-56)  )
  //!  (  (15,-30)  (18,-36)  (17,-34)  (22,-44)  (27,-54)  )
  //!  (  (20,-40)  (23,-46)  (22,-44)  (17,-34)  (22,-44)  )
  //!  (  (25,-50)  (28,-56)  (27,-54)  (22,-44)  (13,-26)  )
  std::cout<<"cm1.transpose() = \n"<<cm1.transpose()<<std::endl;
  //!  cm1.transpose() = 
  //!  5 5
  //!  (  (5,10)  (10,20)  (15,30)  (20,40)  (25,50)  )
  //!  (  (10,20)  (13,26)  (18,36)  (23,46)  (28,56)  )
  //!  (  (15,30)  (18,36)  (17,34)  (22,44)  (27,54)  )
  //!  (  (20,40)  (23,46)  (22,44)  (17,34)  (22,44)  )
  //!  (  (25,50)  (28,56)  (27,54)  (22,44)  (13,26)  )
  std::cout<<"cm1.adjoint() = \n"<<cm1.adjoint()<<std::endl;
  //!  cm1.adjoint() = 
  //!  5 5
  //!  (  (5,-10)  (10,-20)  (15,-30)  (20,-40)  (25,-50)  )
  //!  (  (10,-20)  (13,-26)  (18,-36)  (23,-46)  (28,-56)  )
  //!  (  (15,-30)  (18,-36)  (17,-34)  (22,-44)  (27,-54)  )
  //!  (  (20,-40)  (23,-46)  (22,-44)  (17,-34)  (22,-44)  )
  //!  (  (25,-50)  (28,-56)  (27,-54)  (22,-44)  (13,-26)  )
  std::cout<<"cm1.realPart() = \n"<<cm1.realPart()<<std::endl;
  //!  cm1.realPart() = 
  //!  5 5
  //!  (  5  10  15  20  25  )
  //!  (  10  13  18  23  28  )
  //!  (  15  18  17  22  27  )
  //!  (  20  23  22  17  22  )
  //!  (  25  28  27  22  13  )
  std::cout<<"cm1.imagPart() = \n"<<cm1.imagPart()<<std::endl;
  //!  cm1.imagPart() = 
  //!  5 5
  //!  (  10  20  30  40  50  )
  //!  (  20  26  36  46  56  )
  //!  (  30  36  34  44  54  )
  //!  (  40  46  44  34  44  )
  //!  (  50  56  54  44  26  )
  std::cout<<"Norm(cm1) = "<<Norm(cm1)<<std::endl;
  //!  Norm(cm1) = 227.035

  tmv::Vector<std::complex<double> > cv = v * std::complex<double>(-2,1);
  cm1 += cv^cv;
  std::cout<<"cm1 += cv^cv = \n"<<cm1<<std::endl;
  //!  cm1 += cv^cv = 
  //!  5 5
  //!  (  (17,-6)  (16,12)  (33,6)  (44,8)  (79,-22)  )
  //!  (  (16,12)  (16,22)  (27,24)  (35,30)  (55,20)  )
  //!  (  (33,6)  (27,24)  (44,-2)  (58,-4)  (108,-54)  )
  //!  (  (44,8)  (35,30)  (58,-4)  (65,-30)  (130,-100)  )
  //!  (  (79,-22)  (55,20)  (108,-54)  (130,-100)  (256,-298)  )
  cm1 += v^v;
  std::cout<<"cm1 += v^v = \n"<<cm1<<std::endl;
  //!  cm1 += v^v = 
  //!  5 5
  //!  (  (21,-6)  (18,12)  (39,6)  (52,8)  (97,-22)  )
  //!  (  (18,12)  (17,22)  (30,24)  (39,30)  (64,20)  )
  //!  (  (39,6)  (30,24)  (53,-2)  (70,-4)  (135,-54)  )
  //!  (  (52,8)  (39,30)  (70,-4)  (81,-30)  (166,-100)  )
  //!  (  (97,-22)  (64,20)  (135,-54)  (166,-100)  (337,-298)  )

  tmv::Matrix<std::complex<double> > cx(5,2);
  cx.col(0) = v;
  cx.col(1) = cv;
  cm1 -= cx*cx.transpose();
  std::cout<<"cm1 -= cx*cx.transpose() = \n"<<cm1<<std::endl;
  //!  cm1 -= cx*cx.transpose() = 
  //!  5 5
  //!  (  (5,10)  (10,20)  (15,30)  (20,40)  (25,50)  )
  //!  (  (10,20)  (13,26)  (18,36)  (23,46)  (28,56)  )
  //!  (  (15,30)  (18,36)  (17,34)  (22,44)  (27,54)  )
  //!  (  (20,40)  (23,46)  (22,44)  (17,34)  (22,44)  )
  //!  (  (25,50)  (28,56)  (27,54)  (22,44)  (13,26)  )

  tmv::HermMatrix<std::complex<double> > cm2 = m1;
  cm2.upperTri().offDiag() *= std::complex<double>(1,2);
  std::cout<<"cm2 = \n"<<cm2<<std::endl;
  //!  cm2 = 
  //!  5 5
  //!  (  (5,0)  (10,20)  (15,30)  (20,40)  (25,50)  )
  //!  (  (10,-20)  (13,0)  (18,36)  (23,46)  (28,56)  )
  //!  (  (15,-30)  (18,-36)  (17,0)  (22,44)  (27,54)  )
  //!  (  (20,-40)  (23,-46)  (22,-44)  (17,0)  (22,44)  )
  //!  (  (25,-50)  (28,-56)  (27,-54)  (22,-44)  (13,0)  )
  std::cout<<"cm2.conjugate() = \n"<<cm2.conjugate()<<std::endl;
  //!  cm2.conjugate() = 
  //!  5 5
  //!  (  (5,-0)  (10,-20)  (15,-30)  (20,-40)  (25,-50)  )
  //!  (  (10,20)  (13,-0)  (18,-36)  (23,-46)  (28,-56)  )
  //!  (  (15,30)  (18,36)  (17,-0)  (22,-44)  (27,-54)  )
  //!  (  (20,40)  (23,46)  (22,44)  (17,-0)  (22,-44)  )
  //!  (  (25,50)  (28,56)  (27,54)  (22,44)  (13,-0)  )
  std::cout<<"cm2.transpose() = \n"<<cm2.transpose()<<std::endl;
  //!  cm2.transpose() = 
  //!  5 5
  //!  (  (5,-0)  (10,-20)  (15,-30)  (20,-40)  (25,-50)  )
  //!  (  (10,20)  (13,-0)  (18,-36)  (23,-46)  (28,-56)  )
  //!  (  (15,30)  (18,36)  (17,-0)  (22,-44)  (27,-54)  )
  //!  (  (20,40)  (23,46)  (22,44)  (17,-0)  (22,-44)  )
  //!  (  (25,50)  (28,56)  (27,54)  (22,44)  (13,-0)  )
  std::cout<<"cm2.adjoint() = \n"<<cm2.adjoint()<<std::endl;
  //!  cm2.adjoint() = 
  //!  5 5
  //!  (  (5,0)  (10,20)  (15,30)  (20,40)  (25,50)  )
  //!  (  (10,-20)  (13,0)  (18,36)  (23,46)  (28,56)  )
  //!  (  (15,-30)  (18,-36)  (17,0)  (22,44)  (27,54)  )
  //!  (  (20,-40)  (23,-46)  (22,-44)  (17,0)  (22,44)  )
  //!  (  (25,-50)  (28,-56)  (27,-54)  (22,-44)  (13,0)  )
  std::cout<<"cm2.realPart() = \n"<<cm2.realPart()<<std::endl;
  //!  cm2.realPart() = 
  //!  5 5
  //!  (  5  10  15  20  25  )
  //!  (  10  13  18  23  28  )
  //!  (  15  18  17  22  27  )
  //!  (  20  23  22  17  22  )
  //!  (  25  28  27  22  13  )
  // imagPart is invalid for hermitian matrix, since the result 
  // would be anti-symmetric, which we don't have as a matrix type.
  std::cout<<"Norm(cm2) = "<<Norm(cm2)<<std::endl;
  //!  Norm(cm2) = 218.589

  cm2 += cv^cv.conjugate();
  std::cout<<"cm2 += cv^cv.conjugate() = \n"<<cm2<<std::endl;
  //!  cm2 += cv^cv.conjugate() = 
  //!  5 5
  //!  (  (25,0)  (20,20)  (45,30)  (60,40)  (115,50)  )
  //!  (  (20,-20)  (18,0)  (33,36)  (43,46)  (73,56)  )
  //!  (  (45,-30)  (33,-36)  (62,0)  (82,44)  (162,54)  )
  //!  (  (60,-40)  (43,-46)  (82,-44)  (97,0)  (202,44)  )
  //!  (  (115,-50)  (73,-56)  (162,-54)  (202,-44)  (418,0)  )
  cm2 += v^v;
  std::cout<<"cm2 += v^v = \n"<<cm2<<std::endl;
  //!  cm2 += v^v = 
  //!  5 5
  //!  (  (29,0)  (22,20)  (51,30)  (68,40)  (133,50)  )
  //!  (  (22,-20)  (19,0)  (36,36)  (47,46)  (82,56)  )
  //!  (  (51,-30)  (36,-36)  (71,0)  (94,44)  (189,54)  )
  //!  (  (68,-40)  (47,-46)  (94,-44)  (113,0)  (238,44)  )
  //!  (  (133,-50)  (82,-56)  (189,-54)  (238,-44)  (499,0)  )
  cm2 -= cx*cx.adjoint();
  std::cout<<"cm2 -= cx*cx.adjoint() = \n"<<cm2<<std::endl;
  //!  cm2 -= cx*cx.adjoint() = 
  //!  5 5
  //!  (  (5,0)  (10,20)  (15,30)  (20,40)  (25,50)  )
  //!  (  (10,-20)  (13,0)  (18,36)  (23,46)  (28,56)  )
  //!  (  (15,-30)  (18,-36)  (17,0)  (22,44)  (27,54)  )
  //!  (  (20,-40)  (23,-46)  (22,-44)  (17,0)  (22,44)  )
  //!  (  (25,-50)  (28,-56)  (27,-54)  (22,-44)  (13,0)  )


  // Can mix SymMatrix with other kinds of matrices:
  tmv::UpperTriMatrix<double> um(xm);
  std::cout<<"um + m1 = \n"<<um+m1<<std::endl;
  //!  um + m1 = 
  //!  5  5
  //!  (  7  15  16  29  33  )
  //!  (  10  18  25  25  28  )
  //!  (  15  18  23  30  31  )
  //!  (  20  23  22  26  22  )
  //!  (  25  28  27  22  18  )
  tmv::DiagMatrix<double> dm(xm);
  std::cout<<"dm * m1 = \n"<<dm*m1<<std::endl;
  //!  dm * m1 = 
  //!  5  5
  //!  (  10  20  30  40  50  )
  //!  (  50  65  90  115  140  )
  //!  (  90  108  102  132  162  )
  //!  (  180  207  198  153  198  )
  //!  (  125  140  135  110  65  )

  return 0;
} 
catch (tmv::Error& e) 
{
  std::cerr<<e<<std::endl;
  return 1;
}
