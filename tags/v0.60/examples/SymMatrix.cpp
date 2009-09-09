#include "TMV.h"
// Note: extra include file for SymMatrix
#include "TMV_Sym.h"
// Also need to compile with -ltmv_symband
#include <iostream>

// Only need these, since I use UpperTriMatrix, etc. here.
// You don't need these if you are sticking to regular matrices
// and symmetric matrices.
#include "TMV_Tri.h"
#include "TMV_Diag.h"

int main() try 
{

  // Several ways to create and initialize matrices:
  
  std::cout<<"\nCreate / Initialize:\n\n";
  tmv::SymMatrix<double> m1(5); // Create with uninitialized values
  for(size_t i=0;i<m1.nrows();i++) 
    for(size_t j=0;j<m1.ncols();j++) 
      if (i>=j) m1(i,j) = 2.*i-j*j+10.; 
  std::cout<<"m1 = \n"<<m1<<std::endl;

  tmv::SymMatrix<double> m2(5,2.); // Create with all 2's.
  std::cout<<"m2 = \n"<<m2<<std::endl;

  // Create from given elements:
  double mm[25] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
    16,17,18,19,20,21,22,23,24,25};
  tmv::SymMatrix<double> m3(5,mm); // Default storage is Upper, RowMajor
  std::cout<<"m3 (Upper, RowMajor) = \n"<<m3<<std::endl;
  tmv::SymMatrix<double,tmv::Upper,tmv::ColMajor> m4(5,mm);
  std::cout<<"m4 (Upper,ColMajor) = \n"<<m4<<std::endl;
  tmv::SymMatrix<double,tmv::Lower,tmv::RowMajor> m5(5,mm);
  std::cout<<"m5 (Lower,RowMajor) = \n"<<m5<<std::endl;
  tmv::SymMatrix<double,tmv::Lower,tmv::ColMajor> m6(5,mm);
  std::cout<<"m6 (Lower,ColMajor) = \n"<<m6<<std::endl;

  // Make from the corresponding portion of a regular Matrix
  tmv::Matrix<double> xm(5,5);
  for(size_t i=0;i<xm.nrows();i++) 
    for(size_t j=0;j<xm.ncols();j++) 
      xm(i,j) = 2.*i-j*j+10.; 
  tmv::SymMatrix<double,tmv::Lower,tmv::ColMajor> m7(xm);


  // Norms, etc. 
  // Note: All functions in this section may be written as f(m) or m.f()

  std::cout<<"\nNorms, etc.:\n\n";
  std::cout<<"Norm1(m1) = "<<Norm1(m1)<<std::endl;
  std::cout<<"Norm2(m1) = "<<Norm2(m1)<<std::endl;
  std::cout<<"NormInf(m1) = "<<NormInf(m1)<<std::endl;
  std::cout<<"NormF(m1) = "<<NormF(m1);
  std::cout<<" = "<<Norm(m1)<<std::endl;
  std::cout<<"MaxAbsElement(m1) = "<<MaxAbsElement(m1)<<std::endl;
  std::cout<<"Trace(m1) = "<<Trace(m1)<<std::endl;
  std::cout<<"Det(m1) = "<<Det(m1)<<std::endl;


  // Views:
  
  std::cout<<"\nViews:\n\n";
  std::cout<<"m1.diag() = "<<m1.diag()<<std::endl;
  std::cout<<"m1.diag(1) = "<<m1.diag(1)<<std::endl;
  std::cout<<"m1.diag(-1) = "<<m1.diag(-1)<<std::endl;
  std::cout<<"m1.SubSymMatrix(0,3) = \n"<<m1.SubSymMatrix(0,3)<<std::endl;
  std::cout<<"m1.UpperTri() = \n"<<m1.UpperTri()<<std::endl;
  std::cout<<"m1.LowerTri(tmv::UnitDiag) = \n"<<
    m1.LowerTri(tmv::UnitDiag)<<std::endl;


  // Fortran Indexing:

  std::cout<<"\nFortran Indexing:\n\n";
  tmv::SymMatrix<double,tmv::Upper,tmv::RowMajor,tmv::FortranStyle> fm1 = m1;

  std::cout<<"fm1 = m1 = \n"<<fm1<<std::endl;
  std::cout<<"fm1(1,1) = "<<fm1(1,1)<<std::endl;
  std::cout<<"fm1(4,3) = "<<fm1(4,3)<<std::endl;
  std::cout<<"fm1.SubSymMatrix(1,3) = \n"<<fm1.SubSymMatrix(1,3)<<std::endl;


  // Matrix arithmetic:

  std::cout<<"\nArithmetic:\n\n";
  tmv::SymMatrix<double> m1pm3 = m1 + m3;  // Add matrices
  std::cout<<"m1 + m3 = \n"<<m1pm3<<std::endl;
  // Works correctly even if matrices are stored in different order:
  tmv::SymMatrix<double> m3pm4 = m3 + m4; 
  std::cout<<"m3 + m4 = \n"<<m3pm4<<std::endl;
  tmv::SymMatrix<double> m4pm5 = m4 + m5; 
  std::cout<<"m4 + m5 = \n"<<m4pm5<<std::endl;

  m1 *= 2.;  // Multiply matrix by a scalar
  std::cout<<"m1 *= 2 = \n"<<m1<<std::endl;

  m1 += m4;  // Add one matrix to another
  std::cout<<"m1 += m4 = \n"<<m1<<std::endl;

  // Vector outer product - self outerproduct is symmetric
  tmv::Vector<double> v = xm.col(0);
  tmv::SymMatrix<double> vv = v^v;
  std::cout<<"v = "<<v<<std::endl;
  std::cout<<"v^v = \n"<<vv<<std::endl;

  // Matrix * vector product
  std::cout<<"m1 * v = "<<m1*v<<std::endl;
  std::cout<<"v * m1 = "<<v*m1<<std::endl;

  // Matrix * matrix product
  tmv::Matrix<double> m1xm = m1 * xm; 
  std::cout<<"m1 * xm = \n"<<m1xm<<std::endl;
  // Note: this next product is not symmetric!
  tmv::Matrix<double> m1m5 = m1 * m5; 
  std::cout<<"m1 * m5 = \n"<<m1m5<<std::endl;

  
  // Complex matrices:

  std::cout<<"\nComplex matrices:\n\n";
  tmv::SymMatrix<std::complex<double> > cm1 = m1 * std::complex<double>(1,2);
  
  std::cout<<"cm1 = m1 * (1+2i) = \n"<<cm1<<std::endl;
  std::cout<<"cm1.Conjugate() = \n"<<cm1.Conjugate()<<std::endl;
  std::cout<<"cm1.Transpose() = \n"<<cm1.Transpose()<<std::endl;
  std::cout<<"cm1.Adjoint() = \n"<<cm1.Adjoint()<<std::endl;
  std::cout<<"cm1.Real() = \n"<<cm1.Real()<<std::endl;
  std::cout<<"cm1.Imag() = \n"<<cm1.Imag()<<std::endl;
  std::cout<<"Norm(cm1) = "<<Norm(cm1)<<std::endl;
  tmv::Vector<std::complex<double> > cv = v * std::complex<double>(2,5);
  cm1 += cv^cv;
  std::cout<<"cm1 += cv^cv = \n"<<cm1<<std::endl;
  cm1 += v^v;
  std::cout<<"cm1 += v^v = \n"<<cm1<<std::endl;
  tmv::Matrix<std::complex<double> > cx(5,3);
  cx.col(0) = v;
  cx.col(1) = cv;
  cx.col(2) = 3.*cv;
  cm1 += cx*cx.Transpose();
  std::cout<<"cm1 += cx*cx.Transpose() = \n"<<cm1<<std::endl;


  tmv::HermMatrix<std::complex<double> > cm2 = m1;
  cm2.UpperTri().OffDiag() *= std::complex<double>(1,2);
  std::cout<<"cm2 = \n"<<cm2<<std::endl;
  std::cout<<"cm2.Conjugate() = \n"<<cm2.Conjugate()<<std::endl;
  std::cout<<"cm2.Transpose() = \n"<<cm2.Transpose()<<std::endl;
  std::cout<<"cm2.Adjoint() = \n"<<cm2.Adjoint()<<std::endl;
  std::cout<<"cm2.Real() = \n"<<cm2.Real()<<std::endl;
  std::cout<<"cm2.Imag() = \n"<<cm2.Imag()<<std::endl;
  std::cout<<"Norm(cm2) = "<<Norm(cm2)<<std::endl;
  cm2 += cv^cv.Conjugate();
  std::cout<<"cm2 += cv^cv.Conjugate() = \n"<<cm2<<std::endl;
  cm2 += v^v;
  std::cout<<"cm2 += v^v = \n"<<cm2<<std::endl;
  cm2 += cx*cx.Adjoint();
  std::cout<<"cm2 += cx*cx.Adjoint() = \n"<<cm2<<std::endl;

  // Can mix SymMatrix with other kinds of matrices:
  tmv::UpperTriMatrix<double> um(xm);
  std::cout<<"um + m1 = \n"<<um+m1<<std::endl;
  tmv::DiagMatrix<double> dm(xm);
  std::cout<<"dm * m1 = \n"<<dm*m1<<std::endl;
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

