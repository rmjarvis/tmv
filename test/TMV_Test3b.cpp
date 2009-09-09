// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#include <fstream>
#include "TMV_Test.h"
#include "TMV_Test3.h"

bool XXDEBUG1 = false;
bool XXDEBUG2 = false;
bool XXDEBUG3 = false;
bool XXDEBUG4 = false;
bool XXDEBUG5 = false;
bool XXDEBUG6 = false;
bool XXDEBUG7 = false;
bool XXDEBUG8 = false;
bool XXDEBUG9 = false;

bool showtests = false;
bool showacc = false;
bool showdiv = false;
bool showstartdone = false;
bool donorm2 = true;
bool symoprod = false;
bool dontthrow = false;
std::string lastsuccess = "";

//#include "TMV_Small.h"

int main() try {
  std::ofstream log("tmvtest3b.log");
  tmv::WriteWarningsTo(&log);

  //showacc=true;
  //showdiv=true;
  //showtests=true;
  //showstartdone=true;

#if 0
  tmv::SmallMatrix<double,2,2,tmv::RowMajor> a;
  a = tmv::ListInit, 4.5, 2., 1.2, 8.6;
  tmv::SmallVector<double,2> b;
  b = tmv::ListInit, 3.1, 0.9;
  tmv::SmallVector<double,2> c = b;
  tmv::Vector<double> v = b;
  double x = 5.;

  std::cout<<"a = "<<a<<std::endl;
  std::cout<<"b = "<<b<<std::endl;
  std::cout<<"c = "<<c<<std::endl;
  std::cout<<"v = "<<v<<std::endl;
  std::cout<<"x = "<<x<<std::endl;
 
  std::cout<<"(1) v+x*a*b = "<<v+x*tmv::Matrix<double>(a)*tmv::Vector<double>(b)<<std::endl;
  std::cout<<"(2) v+x*a*b = "<<v+tmv::Vector<double>(x*a*b)<<std::endl;
  std::cout<<"(3) v+x*a*b = "<<tmv::SmallVector<double,2>(v+tmv::SmallVector<double,2>(x*a*b))<<std::endl;
  c += x * a * b;
  std::cout<<"(4) c+=x*a*b = "<<c<<std::endl;
  v += x * a * b;
  std::cout<<"(5) v+=x*a*b = "<<v<<std::endl;
#endif

//#define SKIPREST

#ifndef SKIPREST

#ifdef INST_DOUBLE
  TestAllSmallMatrixA<double>();
#endif

#ifdef INST_FLOAT
  TestAllSmallMatrixA<float>();
#endif

#ifdef INST_LONGDOUBLE
  TestAllSmallMatrixA<long double>();
#endif 

#ifdef INST_INT
  TestAllSmallMatrixA<int>();
#endif 

#endif // SKIPREST

  return 0;
}
#ifndef NOTHROW
catch (tmv::Error& e) {
  std::cerr<<e<<std::endl;
  std::cerr<<"Last successful test was "<<lastsuccess<<std::endl;
  return 1;
}
#endif
catch (std::exception& e) {
  std::cerr<<e.what()<<std::endl;
  std::cerr<<"Last successful test was "<<lastsuccess<<std::endl;
  return 1;
}
catch (...) {
  std::cerr<<"Unknown exception thrown\n";
  std::cerr<<"Last successful test was "<<lastsuccess<<std::endl;
  return 1;
}


void PreAssert(std::string s)
{
  if (showtests) { 
    std::cout<<"Trying: "<<s;  
    std::cout.flush(); 
  } 
}

void DoAssert(bool x, std::string s)
{
  if (x) { 
    if (showtests) std::cout<<"  Passed"<<std::endl;
    lastsuccess = s; 
  } else { 
    if (showtests) std::cout<<"  Failed"<<std::endl;
    if (dontthrow) std::cout<<"Failed test: "<<s<<std::endl;  
    else
#ifdef NOTHROW
    { std::cerr<<"Error in test: "<<s<<std::endl; exit(1); }
#else
    throw tmv::Error("Error in test: ",s);  
#endif
  } 
}

