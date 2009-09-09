#include "TMV.h"
#include <iostream>

int main() try 
{

  // Several ways to create and initialize matrices:
  
  std::cout<<"\nCreate / Initialize:\n\n";
  tmv::Matrix<double> m1(4,3); // Create with uninitialized values
  for(size_t i=0;i<m1.nrows();i++) 
    for(size_t j=0;j<m1.ncols();j++) 
      m1(i,j) = 2.*i-3.*j+10.; 
  std::cout<<"m1 = \n"<<m1<<std::endl;

  tmv::Matrix<double> m2(4,3,2.); // Create with all 2's.
  std::cout<<"m2 = \n"<<m2<<std::endl;

  // Create from given elements:
  double mm[12] = {1,2,3,4, 5,6,7,8, 9,10,11,12};
  tmv::Matrix<double> m3(4,3,mm); // Default order is RowMajor
  std::cout<<"m3 (RowMajor) = \n"<<m3<<std::endl;
  tmv::Matrix<double,tmv::ColMajor> m4(4,3,mm);
  std::cout<<"m4 (ColMajor) = \n"<<m4<<std::endl;

  // Create from STL vector of vectors
  std::vector<std::vector<double> > mm2(3,std::vector<double>(3));
  for(size_t i=0;i<3;i++) 
    for(size_t j=0;j<3;j++) 
      mm2[i][j] = 2.*i+j-3.*i*j;
  tmv::Matrix<double> m5(mm2); 
  std::cout<<"m5 = \n"<<m5<<std::endl;


  // Norms, etc. 
  // Note: All functions in this section may be written as f(m) or m.f()

  std::cout<<"\nNorms, etc.:\n\n";
  std::cout<<"Norm1(m1) = "<<Norm1(m1)<<std::endl;
  std::cout<<"Norm2(m1) = "<<Norm2(m1)<<std::endl;
  std::cout<<"NormInf(m1) = "<<NormInf(m1)<<std::endl;
  std::cout<<"NormF(m1) = "<<NormF(m1);
  std::cout<<" = "<<Norm(m1)<<std::endl;
  std::cout<<"MaxAbsElement(m1) = "<<MaxAbsElement(m1)<<std::endl;
  std::cout<<"Trace(m5) = "<<Trace(m5)<<std::endl;
  std::cout<<"Det(m5) = "<<Det(m5)<<std::endl;


  // Views:
  
  std::cout<<"\nViews:\n\n";
  std::cout<<"m1.row(1) = "<<m1.row(1)<<std::endl;
  std::cout<<"m1.col(2) = "<<m1.col(2)<<std::endl;
  std::cout<<"m1.diag() = "<<m1.diag()<<std::endl;
  std::cout<<"m1.diag(1) = "<<m1.diag(1)<<std::endl;
  std::cout<<"m1.diag(-1) = "<<m1.diag(-1)<<std::endl;
  std::cout<<"m1.SubMatrix(2,4,0,2) = \n"<<m1.SubMatrix(2,4,0,2)<<std::endl;
  std::cout<<"m1.SubMatrix(0,4,1,3,2,1) = \n"<<m1.SubMatrix(0,4,1,3,2,1)<<std::endl;
  std::cout<<"m1.Transpose() = \n"<<m1.Transpose()<<std::endl;
  std::cout<<"m1.Cols(1,3) = \n"<<m1.Cols(1,3)<<std::endl;
  std::cout<<"m1.RowPair(3,0) = \n"<<m1.RowPair(3,0)<<std::endl;

  // Can use the views within expressions
  m1.Rows(0,3) += m5.Transpose();
  std::cout<<"m1.Rows(0,3) += m5.Transpose() => m1 = \n"<<m1<<std::endl;
  m1.row(0) *= 2.;
  std::cout<<"m1.row(0) *= 2 => m1 = \n"<<m1<<std::endl;


  // Fortran Indexing:

  std::cout<<"\nFortran Indexing:\n\n";
  tmv::Matrix<double,tmv::RowMajor,tmv::FortranStyle> fm1 = m1;

  std::cout<<"fm1 = m1 = \n"<<fm1<<std::endl;
  std::cout<<"fm1(1,1) = "<<fm1(1,1)<<std::endl;
  std::cout<<"fm1(4,3) = "<<fm1(4,3)<<std::endl;
  std::cout<<"fm1.row(1) = "<<fm1.row(1)<<std::endl;
  std::cout<<"fm1.col(3) = "<<fm1.col(3)<<std::endl;
  std::cout<<"fm1.SubMatrix(2,3,1,2) = \n"<<fm1.SubMatrix(2,3,1,2)<<std::endl;
  std::cout<<"fm1.SubMatrix(1,3,2,3,2,1) = \n"<<fm1.SubMatrix(1,3,2,3,2,1)<<std::endl;
  std::cout<<"fm1.Cols(1,2) = \n"<<fm1.Cols(1,2)<<std::endl;
  std::cout<<"fm1.RowPair(4,1) = \n"<<fm1.RowPair(4,1)<<std::endl;


  // Matrix arithmetic:

  std::cout<<"\nArithmetic:\n\n";
  tmv::Matrix<double> m1pm3 = m1 + m3;  // Add matrices
  std::cout<<"m1 + m3 = \n"<<m1pm3<<std::endl;
  // Works correctly even if matrices are stored in different order:
  tmv::Matrix<double> m3pm4 = m3 + m4; 
  std::cout<<"m3 + m4 = \n"<<m3pm4<<std::endl;

  m1 *= 2.;  // Multiply matrix by a scalar
  std::cout<<"m1 *= 2 = \n"<<m1<<std::endl;

  m1 += m4;  // Add one matrix to another
  std::cout<<"m1 += m4 = \n"<<m1<<std::endl;

  // Vector outer product
  tmv::Vector<double> v1 = m4.col(0);
  tmv::Vector<double> v2 = m4.row(1);
  tmv::Matrix<double> v1v2 = v1^v2;
  std::cout<<"v1 = "<<v1<<std::endl;
  std::cout<<"v2 = "<<v2<<std::endl;
  std::cout<<"v1^v2 = \n"<<v1v2<<std::endl;
  std::cout<<"ColVectorViewOf(v1)*RowVectorViewOf(v2) = \n"<<
    ColVectorViewOf(v1)*RowVectorViewOf(v2);

  // Matrix * vector product
  std::cout<<"m1 * v2 = "<<m1*v2<<std::endl;
  std::cout<<"v1 * m1 = "<<v1*m1<<std::endl;

  // Matrix * matrix product
  tmv::Matrix<double> m1m5 = m1 * m5; 
  std::cout<<"m1 * m5 = \n"<<m1m5<<std::endl;
  std::cout<<"m1.row(0) * m5.col(2) = "<<m1.row(0)*m5.col(2)<<std::endl;
  std::cout<<"(m1 * m5)(0,2) = "<<m1m5(0,2)<<std::endl;

  // Can handle aliases:
  std::cout<<"m1 + 3*m1-m2 = \n"<<m1+3.*m1-m2<<std::endl;
  m1 += 3.*m1-m2;
  std::cout<<"m1 += 3*m1-m2 = \n"<<m1<<std::endl;
  std::cout<<"m5 * m5 = \n"<<m5*m5<<std::endl;
  m5 *= m5;
  std::cout<<"m5 *= m5 = \n"<<m5<<std::endl;
  
  // Scalars can be treated as a multiple of identity matrix
  m5 += 32.;
  std::cout<<"m5 += 32 = \n"<<m5<<std::endl;

  
  // Complex matrices:

  std::cout<<"\nComplex matrices:\n\n";
  tmv::Matrix<std::complex<double> > cm1 = m1 * std::complex<double>(1,2);
  
  std::cout<<"cm1 = m1 * (1+2i) = \n"<<cm1<<std::endl;
  std::cout<<"cm1.Conjugate() = \n"<<cm1.Conjugate()<<std::endl;
  std::cout<<"cm1.Transpose() = \n"<<cm1.Transpose()<<std::endl;
  std::cout<<"cm1.Adjoint() = \n"<<cm1.Adjoint()<<std::endl;
  std::cout<<"cm1.Real() = \n"<<cm1.Real()<<std::endl;
  std::cout<<"cm1.Imag() = \n"<<cm1.Imag()<<std::endl;
  std::cout<<"Norm(cm1) = "<<Norm(cm1)<<std::endl;
  std::cout<<"cm1*cm1.Adjoint() = \n"<<cm1*cm1.Adjoint()<<std::endl;

  // Can mix real and complex in any way
  std::cout<<"cm1 - m1 = \n"<<cm1 - m1<<std::endl;
  std::cout<<"cm1 * m1.Transpose() * (1-2i) = \n"<<
    cm1 * m1.Transpose() * std::complex<double>(1,-2)<<std::endl;

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

