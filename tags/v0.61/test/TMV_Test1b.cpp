#ifdef NDEBUG
#undef NDEBUG
#endif

#include "TMV_Test.h"
#include "TMV_Test1.h"

bool XXDEBUG1 = false;
bool XXDEBUG2 = false;
bool XXDEBUG3 = false;
bool XXDEBUG4 = false;
bool XXDEBUG5 = false;
bool XXDEBUG6 = false;
bool XXDEBUG7 = false;
bool XXDEBUG8 = false;

bool showtests = false;
bool showacc = false;
bool showdiv = false;
bool showstartdone = false;
bool donorm2 = true;
bool symoprod = false;
bool dontthrow = false;
std::string lastsuccess = "";

int main() try {

  //showacc=true;
  //showdiv=true;
  //showtests=true;
  //showstartdone=true;
//#define SKIPREST

#ifndef SKIPREST

#ifdef INST_DOUBLE
  TestDiagMatrix<double>();
  TestDiagDiv<double>();
#endif // DOUBLE

#ifdef INST_FLOAT
  TestDiagMatrix<float>();
  TestDiagDiv<float>();
#endif // FLOAT

#ifdef INST_LONGDOUBLE
  TestDiagMatrix<long double>();
  TestDiagDiv<long double>();
#endif // LONGDOUBLE

#ifdef INST_INT
  TestDiagMatrix<int>();
#endif  // INT

#endif // SKIPREST

  return 0;
}
catch (tmv::Error& e) {
  std::cerr<<e<<std::endl;
  std::cerr<<"Last successful test was "<<lastsuccess<<std::endl;
  return 1;
}
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
    else throw tmv::Error("Error in test: ",s.c_str());  
  } 
}

