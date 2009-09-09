#include "TMV.h"
#include "TMV_Diag.h"
#include <iostream>

int main() try 
{

  tmv::Matrix<double> A(4,4);
  for(size_t i=0;i<A.nrows();i++) 
    for(size_t j=0;j<A.ncols();j++) 
      A(i,j) = 6.*i-2*j*j+2.; 
  A.diag().AddToAll(5.);

  tmv::Vector<double> b(4);
  for(size_t i=0;i<b.size();i++) 
    b(i) = 3.+2.*i; 

  // Basic Ax=b solution:
  std::cout<<"A = \n"<<A<<std::endl;
  std::cout<<"b = \n"<<b<<std::endl;
  std::cout<<"Solve: Ax = b:\n";
  tmv::Vector<double> x = b/A; // Default: use LU decomposition
  std::cout<<"x = b/A = "<<x<<std::endl;
  std::cout<<"A*x = "<<A*x<<std::endl;

  // Can update A and then re-solve:
  A(0,0) = 20.;
  x = b/A;
  std::cout<<"A(0,0) = 20 => x = "<<x<<std::endl;
  std::cout<<"A*x = "<<A*x<<std::endl;

  // If matrix won't change, but you want to solve with multiple 
  // vectors, then it's faster to let TMV know this with SaveDiv().
  tmv::Matrix<double> A2 = A;
  A2.SaveDiv();
  x = b/A2;
  std::cout<<"\n\nx = "<<x<<std::endl;
  std::cout<<"A*x = "<<A2*x<<std::endl;
  b(0) = 10.;
  x = b/A2;  // This doesn't recalculate LU decomposition.
  std::cout<<"b(0) = 10 => x = "<<x<<std::endl;
  std::cout<<"A*x = "<<A2*x<<std::endl;

  // Matrix inverse:
  tmv::Matrix<double> A2inv = A2.Inverse();
  std::cout<<"Ainv = \n"<<A2inv<<std::endl;
  std::cout<<"Ainv*A = \n"<<A2inv*A2<<std::endl;
  std::cout<<"(Ainv*A).Write(std::cout,1.e-8) = \n";
  (A2inv*A2).Write(std::cout,1.e-8);

  // 1/x notation is treated as arithmetic:
  std::cout<<"1./A = \n"<<1./A2<<std::endl;
  std::cout<<"5./A = \n"<<5./A2<<std::endl;

  // Can also use Inverse() notation instead of /
  // x = b/A2 inlines to exactly the same thing as x = A2.Inverse() * b;
  // ie. accurate back-substitution methods are used rather than 
  // actually computing the inverse and then multiplying:
  x =  A2.Inverse() * b;
  std::cout<<"x = A.Inverse() * b = "<<x<<std::endl;

  // Division from the other side can either be done with this
  // Inverse() notation or with the % operator:
  // This is the solution to the equation x A = b, rather than A x = b.
  x = b * A2.Inverse();
  std::cout<<"x = b * A.Inverse() = "<<x<<std::endl;
  x = b % A2;
  std::cout<<"x = b % A = "<<x<<std::endl;
  
  tmv::Matrix<double> B(4,3);
  for(size_t i=0;i<B.nrows();i++) 
    for(size_t j=0;j<B.ncols();j++) 
      B(i,j) = 1.+2.*i+j*j;
  // Multiple right hand sides may be calculated at once if 
  // B is a matrix, rather than a vector:
  std::cout<<"B = \n"<<B<<std::endl;
  tmv::Matrix<double> X = B/A2;
  std::cout<<"X = B/A = \n"<<X<<std::endl;
  std::cout<<"AX = "<<A2*X<<std::endl;

  // And as always, you can mix complex and real objects:
  tmv::Vector<std::complex<double> > cb = b * std::complex<double>(3,-2);
  cb(1) = std::complex<double>(-1,8);
  cb(2) = std::complex<double>(1,1);
  std::cout<<"cb = "<<cb<<std::endl;
  std::cout<<"cb/A = "<<cb/A<<std::endl;
  tmv::Matrix<std::complex<double> > CA = A * std::complex<double>(5,-2);
  CA(1,1) = std::complex<double>(1,6);
  CA(2,3) = std::complex<double>(4,-1);
  std::cout<<"CA = "<<CA<<std::endl;
  std::cout<<"b/CA = "<<b/CA<<std::endl;
  std::cout<<"cb/CA = "<<cb/CA<<std::endl;

  // Least-squares solutions:
  // If A in the matrix equation A x = b has more rows than columns,
  // then there is, in general, no solution to the equation.
  // Instead, one is generally looking for the x that comes closest
  // to satisfying the equation, in a least-squares sense.
  // Specifically, the x for which Norm2(b-Ax) is minimized.
  // This is the solution produced by TMV for such matrices.
  //
  // Here I model a theoretical system for which each observation
  // is 5 + 6i - 3i^2 in the absense of measurement errors.
  // We observe the actual values which are not quite equal to 
  // that because of noise.  And the goal is to determine the 
  // coefficients of 1,i,i^2 (5,6,-3) from noisy observations:
  tmv::Vector<double> b3(10);
  tmv::Matrix<double> A3(10,3);
  double errors[10] = {0.01,-0.02,0.02,-0.02,0.00,-0.03,0.01,-0.02,0.03,0.02};
  double sigma = 0.02;
  for(int i=0;i<10;i++) {
    b3(i) = 5. + 6.*i - 3.*i*i + errors[i]; // Model of measurements
    A3(i,0) = 1.;  // Parameterization of the model...
    A3(i,1) = i;
    A3(i,2) = i*i;
  }
  std::cout<<"Solve Ax=b with:\n";
  std::cout<<"A*sigma = \n"<<A3<<std::endl;
  std::cout<<"b*sigma = \n"<<b3<<std::endl;
  A3 /= sigma;
  b3 /= sigma;

  tmv::Vector<double> x3 = b3/A3;  // Uses QR decomposition by default
  std::cout<<"\nx = "<<x3<<std::endl;
  std::cout<<"A*x => "<<A3*x3<<std::endl;
  std::cout<<"chisq = NormSq(A*x-b) = "<<NormSq(A3*x3-b3)<<std::endl;
  // (The expected value for this is 10-3 = 7.)

  // The covariance matrix for x is (A.Transpose() * A)^-1
  // This combination is easy to calculate from the QR decomposition
  // that TMV has used to do the division.  Therefore, we provide
  // it as an explicit function:
  tmv::Matrix<double> cov(3,3);
  A3.InverseATA(cov);
  std::cout<<"Cov(x) = \n"<<cov<<std::endl;
  
  // SVD can detect ill-conditioned matrices and correct for them.
  // For example, if you model the above observations with the 
  // 1,i,i^2 components as before, but add as well 6*i-5 as a component,
  // then that would be degenerate with 1 and i.  
  // SVD is able to detect this defect and deal with it appropriately:
  tmv::Matrix<double> A4(10,4);
  A4.Cols(0,3) = A3*sigma;
  for(int i=0;i<10;i++) A4(i,3) = 6.*i-5.;
  std::cout<<"Now A*sigma = \n"<<A4<<std::endl;
  A4 /= sigma;
  // This may or may not succeed, but if it does, the results will be 
  // unusable, typically with values around 1.e12 and such.
  try {
    tmv::Vector<double> x4 = b3/A4;
    std::cout<<"Successful x = b/A = \n"<<x4<<std::endl;
  } catch (tmv::Error& e) {
    std::cout<<"Tried x = b/A, but caught error: \n"<<e<<std::endl;
  }
    
  A4.DivideUsing(tmv::SV);
  A4.SetDiv();
  std::cout<<"Singular values for A are "<<A4.SVD().GetS()<<std::endl;
  std::cout<<"Using only the first "<<A4.SVD().GetKMax()<<" components\n";
  tmv::Vector<double> x4 = b3/A4;
  std::cout<<"SVD division yields: x = "<<x4<<std::endl;
  std::cout<<"chisq = "<<NormSq(A4*x4-b3);
  std::cout<<"      Norm(x) = "<<Norm(x4)<<std::endl;

  // QRP can also give useful results, but isn't quite as flexible 
  // as SVD:
  A4.DivideUsing(tmv::QRP);
  x4 = b3/A4;
  std::cout<<"QRP division yields: x = "<<x4<<std::endl;
  std::cout<<"chisq = "<<NormSq(A4*x4-b3);
  std::cout<<"      Norm(x) = "<<Norm(x4)<<std::endl;
  // Note that both methods give answers with an equal chisq, so
  // Ax is equally close to b for each but they have different 
  // specific choices for the degeneracy.
  // SVD will give the solution within this degeneracy freedom that
  // has the minimum Norm(x), but QRP will be faster, significantly so
  // for large matrices.

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

