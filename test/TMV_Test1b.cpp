
#define TESTDIV
#define TESTDIAG
#define TESTTRI
#define TESTMEM

#ifdef TESTMEM
#define MEMDEBUG
#include "MemDebug.h"
AllocList* allocList=0;
#endif

#include "TMV_Test.h"
#include "TMV_Test1.h"

bool showtests = false;
bool showacc = false;
bool showdiv = false;
bool showstartdone = false;
bool donorm2 = false;
bool symoprod = false;
bool dontthrow = false;
std::string lastsuccess = "";

int main() try {

#ifdef TESTMEM
  atexit(&DumpUnfreed);
#endif

#ifdef XTEST
  donorm2 = true;
#endif

  //showacc=true;
  //showdiv=true;
  //showtests=true;
  //showstartdone=true;
//#define SKIPREST

#ifndef SKIPREST

#ifdef INST_DOUBLE
#ifdef TESTDIAG
  TestDiagMatrix<double>();
  TestDiagDiv<double>();
#endif
#endif // DOUBLE

#ifdef INST_FLOAT
#ifdef TESTDIAG
  TestDiagMatrix<float>();
  TestDiagDiv<float>();
#endif
#endif // FLOAT

#ifdef INST_LONGDOUBLE
#ifdef TESTDIAG
  TestDiagMatrix<long double>();
  TestDiagDiv<long double>();
#endif
#endif // LONGDOUBLE

#ifdef INST_INT
#ifdef TESTDIAG
  TestDiagMatrix<int>();
#endif
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
    else throw tmv::Error(std::string("Error in test: ") + s);  
  } 
}

