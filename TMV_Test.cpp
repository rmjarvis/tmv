
#define TESTDIV
#define TESTDIAG
#define TESTTRI
#define TESTBAND
#define TESTSYM
#define TESTMEM

#ifdef TESTMEM
#define MEMDEBUG
#include <valarray>
// With gcc, it seems that valarray needs to be included before MemDebug
#include "MemDebug.h"
AllocList* allocList=0;
#endif

#include "TMV_Test.h"

#ifdef TESTDIAG
#include "TMV_Diag.h"
#endif

#ifdef TESTTRI
#include "TMV_Tri.h"
#endif

#ifdef TESTBAND
#include "TMV_Band.h"
#endif

#ifdef TESTSYM
#include "TMV_Sym.h"
#endif

bool showtests = false;
bool showacc = false;
bool showdiv = false;
bool doallarith = false;
bool showstartdone = false;
bool donorm2 = true;

int main() {

  try {

#ifdef TESTMEM
    atexit(&DumpUnfreed);
#endif

#ifdef XTEST
    doallarith = true;
    donorm2 = true;
#endif

    //showacc=true;
    //showdiv=true;
    //showtests=true;
    //showstartdone=true;
//#define SKIPREST

#ifndef SKIPREST
    TestAllVector<double>();
    TestAllMatrix<double>();
#ifdef TESTDIAG
    TestDiagMatrix<double>();
    TestDiagDiv<double>();
#endif
#ifdef TESTTRI
    TestTriMatrix<double>();
    TestAllTriDiv<double>();
#endif
#ifdef TESTDIV
    TestAllMatrixDiv<double>();
#endif
#ifdef TESTBAND
    TestBandMatrix<double>();
    TestAllBandDiv<double>();
#endif
#ifdef TESTSYM
    TestSymMatrix<double>();
    TestAllSymDiv<double>();
#endif

#ifndef NOFLOAT
    TestAllVector<float>();
    TestAllMatrix<float>();
#ifdef TESTDIV
    TestAllMatrixDiv<float>();
#endif
#ifdef TESTDIAG
    TestDiagMatrix<float>();
    TestDiagDiv<float>();
#endif
#ifdef TESTTRI
    TestTriMatrix<float>();
    TestAllTriDiv<float>();
#endif
#ifdef TESTBAND
    TestBandMatrix<float>();
    TestAllBandDiv<float>();
#endif
#ifdef TESTSYM
    TestSymMatrix<float>();
    TestAllSymDiv<float>();
#endif
#endif

#ifdef LONGDOUBLE
    TestAllVector<long double>();
    TestAllMatrix<long double>();
#ifdef TESTDIV
    TestAllMatrixDiv<long double>();
#endif
#ifdef TESTDIAG
    TestDiagMatrix<long double>();
    TestDiagDiv<long double>();
#endif
#ifdef TESTTRI
    TestTriMatrix<long double>();
    TestAllTriDiv<long double>();
#endif
#ifdef TESTBAND
    TestBandMatrix<long double>();
    TestAllBandDiv<long double>();
#endif
#ifdef TESTSYM
    TestSymMatrix<long double>();
    TestAllSymDiv<long double>();
#endif
#endif 

#endif // SKIPREST

    return 0;
  }
  catch (tmv::Error& e) {
    std::cerr<<e<<std::endl;
    return 1;
  }
  catch (std::exception& e) {
    std::cerr<<e.what()<<std::endl;
    return 1;
  }
  catch (...) {
    std::cerr<<"Unknown exception thrown\n";
    return 1;
  }

}
