#include "TMV.h"
#include <iostream>

int main() try 
{

  // Several ways to create and initialize vectors:
  
  std::cout<<"\nCreate / Initialize:\n\n";
  tmv::Vector<double> v1(6); // Create with uninitialized values
  double val = 1.;
  typedef tmv::Vector<double>::iterator vit;
  for(vit it = v1.begin(); it != v1.end(); it++) {
    *it = val;
    val *= 2.;
  }
  std::cout<<"v1 = "<<v1<<std::endl;

  tmv::Vector<double> v2(6,2.); // Create with all 2's.
  std::cout<<"v2 = "<<v2<<std::endl;

  double vv[6] = {1.1,8.,-15.,2.5,6.3,-12.};
  tmv::Vector<double> v3(6,vv); // Create with given elements;
  std::cout<<"v3 = "<<v3<<std::endl;

  tmv::Vector<double> v4(6); // Another uninitialized vector
  for(size_t i=0;i<v4.size();i++) {
    v4(i) = 2*i+10.;  // Could also use v4[i] instead
  }
  std::cout<<"v4 = "<<v4<<std::endl;


  // Norms, etc. 
  // Note: All functions in this section may be written as f(v) or v.f()

  std::cout<<"\nNorms, etc.:\n\n";
  std::cout<<"Norm1(v3) = "<<Norm1(v3);
  std::cout<<" = "<<v3.SumAbsElements()<<std::endl;

  std::cout<<"Norm2(v3) = "<<v3.Norm2();
  std::cout<<" = "<<v3.Norm()<<std::endl;

  std::cout<<"NormInf(v3) = "<<v3.NormInf();
  std::cout<<" = "<<MaxAbsElement(v3)<<std::endl;

  std::cout<<"SumElements(v3) = "<<SumElements(v3)<<std::endl;

  // Min/Max elements:
  size_t i1,i2,i3,i4;
  double x1 = v3.MinAbsElement(&i1);
  double x2 = v3.MaxAbsElement(&i2);
  double x3 = v3.MinElement(&i3);
  double x4 = v3.MaxElement(&i4);
  std::cout<<"|v3("<<i1<<")| = "<<x1<<" is the minimum absolute value\n";
  std::cout<<"|v3("<<i2<<")| = "<<x2<<" is the maximum absolute value\n";
  std::cout<<"v3("<<i3<<") = "<<x3<<" is the minimum value\n";
  std::cout<<"v3("<<i4<<") = "<<x4<<" is the maximum value\n";


  // Modifications:
  
  std::cout<<"\nModifications:\n\n";
  std::cout<<"v1 = "<<v1<<std::endl;
  v1.AddToAll(5.);
  std::cout<<"v1.AddToAll(5.) = \n"<<v1<<std::endl;

  v1.ReverseSelf();
  std::cout<<"v1.ReverseSelf() = \n"<<v1<<std::endl;

  v1.Zero();
  std::cout<<"v1.Zero() = \n"<<v1<<std::endl;

  v1.SetAllTo(20.);
  std::cout<<"v1.SetAllTo(20.) = \n"<<v1<<std::endl;

  v1.MakeBasis(2);
  std::cout<<"v1.MakeBasis(2) = \n"<<v1<<std::endl;


  // Views:
  
  std::cout<<"\nViews:\n\n";
  std::cout<<"v3.SubVector(0,3) = \n"<<v3.SubVector(0,3)<<std::endl;
  std::cout<<"v3.SubVector(0,6,2) = \n"<<v3.SubVector(0,6,2)<<std::endl;
  std::cout<<"v3.Reverse() = \n"<<v3.Reverse()<<std::endl;

  // Can use the views within expressions
  v3.Reverse() += v1;
  std::cout<<"v3.Reverse() += v1 => v3 = \n"<<v3<<std::endl;
  v3.SubVector(0,5) *= 2.;
  std::cout<<"v3.SubVector(0,5) *= 2 => v3 = \n"<<v3<<std::endl;

  // Fortran Indexing:

  std::cout<<"\nFortran Indexing:\n\n";
  tmv::Vector<double,tmv::FortranStyle> fv3 = v3;

  std::cout<<"fv3 = v3 = "<<fv3<<std::endl;
  std::cout<<"fv3(1) = "<<fv3(1)<<std::endl;
  std::cout<<"fv3(6) = "<<fv3(6)<<std::endl;
  std::cout<<"fv3.SubVector(1,3) = "<<fv3.SubVector(1,3)<<std::endl;
  std::cout<<"fv3.SubVector(1,5,2) = "<<fv3.SubVector(1,5,2)<<std::endl;
  std::cout<<"fv3.MakeBasis(2) = "<<fv3.MakeBasis(2)<<std::endl;
     // Note that this puts the 1 in the 2nd spot, rather than the 
     // 3rd spot that the CStyle version did.


  // Vector arithmetic:

  std::cout<<"\nArithmetic:\n\n";
  tmv::Vector<double> v3pv4 = v3 + v4;  // Add element by element
  std::cout<<"v3 + v4 = "<<v3pv4<<std::endl;

  double v3v4 = v3 * v4;  // Inner product
  std::cout<<"v3 * v4 = "<<v3v4<<std::endl;

  v3 *= 2.;  // Multiply vector by a scalar
  std::cout<<"v3 *= 2 = "<<v3<<std::endl;

  v3 += v4;  // Add one vector to another
  std::cout<<"v3 += v4 = "<<v3<<std::endl;

  // Get as complicated as you want:
  std::cout<<"(v1*v2) * v3 + ((-v4 + 4.*v1)*v2) * v2/20. = \n"<<
    (v1*v2) * v3 + ((-v4 + 4.*v1)*v2) * v2/20.<<std::endl;

  // Automatically checks for aliases:
  v3 = v4 - 3.*v3;
  std::cout<<"v3 = v4-3.*v3 => v3 = \n"<<v3<<std::endl;


  // Complex vectors:

  std::cout<<"\nComplex Vectors:\n\n";
  tmv::Vector<std::complex<double> > cv3 = v3 * std::complex<double>(1,2);
  
  std::cout<<"cv3 = v3 * (1+2i) = \n"<<cv3<<std::endl;
  std::cout<<"cv3.Conjugate() = \n"<<cv3.Conjugate()<<std::endl;
  std::cout<<"cv3.Real() = \n"<<cv3.Real()<<std::endl;
  std::cout<<"cv3.Imag() = \n"<<cv3.Imag()<<std::endl;
  std::cout<<"Norm(cv3) = "<<Norm(cv3)<<std::endl;
  std::cout<<"sqrt(cv3*cv3.Conjugate()) = "<<sqrt(cv3*cv3.Conjugate())<<std::endl;
  std::cout<<"cv3.MaxAbsElement() = "<<cv3.MaxAbsElement()<<std::endl;

  // Can mix real and complex in any way
  std::cout<<"cv3 - v3 = \n"<<cv3 - v3<<std::endl;
  std::cout<<"cv3 * v3 * (1-2i) = "<<
    cv3 * v3 * std::complex<double>(1,-2)<<std::endl;


  // Sorting:

  std::cout<<"\nSorting:\n\n";
  v4[0] = 2.;
  v4[1] = 5.3;
  v4[2] = -1.5;
  v4[3] = -7.;
  v4[4] = 0.5;
  v4[5] = -2.8;
  std::cout<<"v4 = "<<v4<<std::endl;
  size_t p[6];
  v4.Sort(p);
  std::cout<<"Sorted: v4 = "<<v4<<std::endl;
  v4.ReversePermute(p);
  std::cout<<"Sort undone: v4 = "<<v4<<std::endl;
  v4.Sort(); // Don't necessarily need p.
  std::cout<<"Resorted: v4 = "<<v4<<std::endl;
  
  // Can sort by other criteria:
  std::cout<<"v4.Sort(0,DESCEND) = \n"<<v4.Sort(0,tmv::DESCEND)<<std::endl;
  tmv::Vector<std::complex<double> > cv4(6);
  cv4[0] = std::complex<double>(-3.,4.);
  cv4[1] = std::complex<double>(1.,-1.);
  cv4[2] = std::complex<double>(-2.,0.);
  cv4[3] = std::complex<double>(-1.,-6.);
  cv4[4] = std::complex<double>(7.,5.);
  cv4[5] = std::complex<double>(3.,-1.);
  std::cout<<"cv4 = \n"<<cv4<<std::endl;
  std::cout<<"cv4.Sort(0,DESCEND,REAL_COMP) = \n"<<
    cv4.Sort(0,tmv::DESCEND,tmv::REAL_COMP)<<std::endl;
  std::cout<<"cv4.Sort(0,ASCEND,IMAG_COMP) = \n"<<
    cv4.Sort(0,tmv::ASCEND,tmv::IMAG_COMP)<<std::endl;
  std::cout<<"cv4.Sort(0,ASCEND,ABS_COMP) = \n"<<
    cv4.Sort(0,tmv::ASCEND,tmv::ABS_COMP)<<std::endl;
  std::cout<<"cv4.Sort(0,ASCEND,ARG_COMP) = \n"<<
    cv4.Sort(0,tmv::ASCEND,tmv::ARG_COMP)<<std::endl;

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

