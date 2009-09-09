#include "TMV.h"
// Note: extra include file for BandMatrix
#include "TMV_Band.h"
// Also need to compile with -ltmv_symband
#include <iostream>

// Only need these, since I use UpperTriMatrix, etc. here.
// You don't need these if you are sticking to regular matrices
// and band matrices.
#include "TMV_Tri.h"
#include "TMV_Diag.h"

int main() try 
{

  // Several ways to create and initialize band matrices:
  
  std::cout<<"\nCreate / Initialize:\n\n";
  tmv::BandMatrix<double> m1(6,6,1,2); // Create with uninitialized values
  for(size_t i=0;i<m1.nrows();i++) 
    for(size_t j=0;j<m1.ncols();j++) 
      if (i<=j+m1.nlo() && j<=i+m1.nhi())
	m1(i,j) = 3.*i-j*j+7.; 
  std::cout<<"m1 = \n"<<m1<<std::endl;

  tmv::BandMatrix<double> m2(6,6,1,3,2.); // Create with all 2's.
  std::cout<<"m2 = \n"<<m2<<std::endl;

  // A BandMatrix can be non-square:
  tmv::BandMatrix<double> m3(6,8,1,3,2.);
  std::cout<<"m3 = \n"<<m3<<std::endl;

  // Create from given elements:
  double mm[21] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21};
  tmv::BandMatrix<double> m4(6,6,2,1,mm); // Default order is RowMajor
  std::cout<<"m4 (RowMajor) = \n"<<m4<<std::endl;
  tmv::BandMatrix<double,tmv::ColMajor> m5(6,6,2,1,mm);
  std::cout<<"m5 (ColMajor) = \n"<<m5<<std::endl;
  tmv::BandMatrix<double,tmv::DiagMajor> m6(6,6,2,1,mm);
  std::cout<<"m6 (DiagMajor) = \n"<<m6<<std::endl;

  // Can make from the banded portion of a regular Matrix:
  tmv::Matrix<double> xm(6,6);
  for(size_t i=0;i<xm.nrows();i++) 
    for(size_t j=0;j<xm.ncols();j++) 
      xm(i,j) = 5.*i-j*j+3.; 
  tmv::BandMatrix<double> m7(xm,3,2);
  std::cout<<"m7 = \n"<<m7<<std::endl;
  // Or from a wider BandMatrix:
  tmv::BandMatrix<double> m8(m7,3,0);
  std::cout<<"m8 = \n"<<m8<<std::endl;

  // Shortcuts to Bi- and Tri-diagonal matrices:
  tmv::Vector<double> v1(5,1.);
  tmv::Vector<double> v2(6,2.);
  tmv::Vector<double> v3(5,3.);
  tmv::BandMatrix<double> m9 = LowerBiDiagMatrix(v1,v2);
  tmv::BandMatrix<double> m10 = UpperBiDiagMatrix(v2,v3);
  tmv::BandMatrix<double> m11 = TriDiagMatrix(v1,v2,v3);
  std::cout<<"LowerBiDiagMatrix(v1,v2) = \n"<<m9<<std::endl;
  std::cout<<"UpperBiDiagMatrix(v2,v3) = \n"<<m10<<std::endl;
  std::cout<<"TriDiagMatrix(v1,v2,v3) = \n"<<m11<<std::endl;


  // Norms, etc. 
  // Note: All functions in this section may be written as f(m) or m.f()

  std::cout<<"\nNorms, etc.:\n\n";
  std::cout<<"Norm1(m1) = "<<Norm1(m1)<<std::endl;
  std::cout<<"Norm2(m1) = "<<Norm2(m1)<<std::endl;
  std::cout<<"NormInf(m1) = "<<NormInf(m1)<<std::endl;
  std::cout<<"NormF(m1) = "<<NormF(m1)<<" = "<<Norm(m1)<<std::endl;
  std::cout<<"MaxAbsElement(m1) = "<<MaxAbsElement(m1)<<std::endl;
  std::cout<<"Trace(m1) = "<<Trace(m1)<<std::endl;
  std::cout<<"Det(m1) = "<<Det(m1)<<std::endl;


  // Views:
  
  std::cout<<"\nViews:\n\n";
  std::cout<<"m1.diag() = "<<m1.diag()<<std::endl;
  std::cout<<"m1.diag(1) = "<<m1.diag(1)<<std::endl;
  std::cout<<"m1.diag(-1) = "<<m1.diag(-1)<<std::endl;
  std::cout<<"m1.SubBandMatrix(0,3,0,3,1,1) = \n"<<
    m1.SubBandMatrix(0,3,0,3,1,1)<<std::endl;
  std::cout<<"m1.Transpose() = \n"<<m1.Transpose()<<std::endl;
  // Rows, Cols shrink both dimensions of the matrix to include only
  // the portions that are in those rows or columns:
  std::cout<<"m1.Rows(0,4) = \n"<<m1.Rows(0,4)<<std::endl;
  std::cout<<"m1.Cols(1,4) = \n"<<m1.Cols(1,4)<<std::endl;
  std::cout<<"m1.Diags(0,2) = \n"<<m1.Diags(0,2)<<std::endl;
  std::cout<<"m1.Diags(-1,1) = \n"<<m1.Diags(-1,1)<<std::endl;


  // Fortran Indexing:

  std::cout<<"\nFortran Indexing:\n\n";
  tmv::BandMatrix<double,tmv::RowMajor,tmv::FortranStyle> fm1 = m1;

  std::cout<<"fm1 = m1 = \n"<<fm1<<std::endl;
  std::cout<<"fm1(1,1) = "<<fm1(1,1)<<std::endl;
  std::cout<<"fm1(4,3) = "<<fm1(4,3)<<std::endl;
  std::cout<<"fm1.SubBandMatrix(1,3,1,3,1,1) = \n"<<
    fm1.SubBandMatrix(1,3,1,3,1,1)<<std::endl;
  std::cout<<"fm1.Rows(1,4) = \n"<<fm1.Rows(1,4)<<std::endl;
  std::cout<<"fm1.Cols(2,4) = \n"<<fm1.Cols(2,4)<<std::endl;
  std::cout<<"fm1.Diags(0,1) = \n"<<fm1.Diags(0,1)<<std::endl;
  std::cout<<"fm1.Diags(-1,0) = \n"<<fm1.Diags(-1,0)<<std::endl;


  // Matrix arithmetic:

  std::cout<<"\nArithmetic:\n\n";
  tmv::BandMatrix<double> m1pm2 = m1 + m2;  // Add matrices
  std::cout<<"m1 + m2 = \n"<<m1pm2<<std::endl;
  // Works correctly even if matrices are stored in different order:
  tmv::BandMatrix<double> m5pm6 = m5 + m6; 
  std::cout<<"m5 + m6 = \n"<<m5pm6<<std::endl;
  // Also expands the number of off-diagonals appropriately as needed:
  tmv::BandMatrix<double> m2pm4 = m2 + m4; 
  std::cout<<"m2 + m4 = \n"<<m2pm4<<std::endl;
  

  m1 *= 2.;  // Multiply matrix by a scalar
  std::cout<<"m1 *= 2 = \n"<<m1<<std::endl;

  m2 += m1;  // Add one matrix to another
  std::cout<<"m2 += m1 = \n"<<m2<<std::endl;

  // Matrix * vector product
  tmv::Vector<double> v = xm.col(0);
  std::cout<<"v = "<<v<<std::endl;
  std::cout<<"m1 * v = "<<m1*v<<std::endl;
  std::cout<<"v * m1 = "<<v*m1<<std::endl;

  // Matrix * matrix product also expands bands appropriately:
  tmv::BandMatrix<double> m1m2 = m1 * m2; 
  std::cout<<"m1 * m2 = \n"<<m1m2<<std::endl;

  // Can mix BandMatrix with other kinds of matrices:
  std::cout<<"xm * m1 = \n"<<xm*m1<<std::endl;
  tmv::UpperTriMatrix<double> um(xm);
  std::cout<<"um + m1 = \n"<<um+m1<<std::endl;
  tmv::LowerTriMatrix<double> lm(xm);
  lm *= m8;
  std::cout<<"lm *= m8 = \n"<<lm<<std::endl;
  tmv::DiagMatrix<double> dm(xm);
  m1 *= dm;
  std::cout<<"m1 *= dm = \n"<<m1<<std::endl;
  
  return 0;
} 
catch (tmv::Error& e) 
{
  std::cerr<<e<<std::endl;
  return 1;
}
catch (std::exception& e) 
{
  std::cerr<<e.what()<<std::endl;
  return 1;
}
catch (...) 
{
  std::cerr<<"Unknown exception thrown\n";
  return 1;
}

