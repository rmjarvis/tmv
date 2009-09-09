
#define TESTDIV

#include "TMV_Test.h"

bool showtests = false;
bool showacc = false;
bool showdiv = false;
bool doallarith = false;
bool showstartdone = false;
bool donorm2 = true;

int main() try {

//#define SKIPREST

#ifndef SKIPREST
  TestAllVector<double>();
  TestAllMatrix<double>();
#ifdef TESTDIV
  TestAllMatrixDiv<double>();
#endif

#ifndef NOFLOAT
  TestAllVector<float>();
  TestAllMatrix<float>();
#ifdef TESTDIV
  TestAllMatrixDiv<float>();
#endif
#endif

#ifdef LONGDOUBLE
  TestAllVector<long double>();
  TestAllMatrix<long double>();
#ifdef TESTDIV
  TestAllMatrixDiv<long double>();
#endif
#endif 

#endif // SKIPREST

  return 0;
}
catch (tmv::Error& e) {
  cout<<e<<endl;
  exit(1);
}

