#ifndef TMV_TEST_H
#define TMV_TEST_H

//#define LONGDOUBLE

#include <iostream>
#include <complex>
#include <vector>
#include <string>

using std::cerr;
using std::cout;
using std::endl;
using std::complex;
using std::abs;
using std::vector;
using std::ostream;
using std::string;

#define EPS (10*tmv::Epsilon<T>())

inline void myerror(const char* s, const char* s2="", const char* s3="")
{
  cout << "Error in test: "<< s << s2 << s3 << endl;
  exit(1);
}

inline void Assert(bool x,string s)
{
  if (!(x)) myerror(s.c_str()); 
#ifdef SHOWTESTS
  else cout<<"Passed: "<<s<<endl;
#endif
}

#endif
