
#include <fstream>
#include "TMV_Test.h"
#include "TMV_Test_2.h"

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
bool symoprod = true;
bool dontthrow = false;
std::string lastsuccess = "";

int main() try 
{
    std::ofstream log("tmvtest2.log");
    tmv::WriteWarningsTo(&log);

    //showtests=true;
    //showacc=true;
    //showdiv=true;
    //showstartdone=true;

#if 1

#ifdef TEST_DOUBLE
    TestBandMatrix<double>();
    TestSymMatrix<double>();
    TestSymBandMatrix<double>();
    TestAllBandDiv<double>();
    TestAllSymDiv<double>();
    TestAllSymBandDiv<double>();
#endif

#ifdef TEST_FLOAT
    TestBandMatrix<float>();
    TestSymMatrix<float>();
    TestSymBandMatrix<float>();
    TestAllBandDiv<float>();
    TestAllSymDiv<float>();
    TestAllSymBandDiv<float>();
#endif

#ifdef TEST_LONGDOUBLE
    TestBandMatrix<long double>();
    TestSymMatrix<long double>();
    TestSymBandMatrix<long double>();
    TestAllBandDiv<long double>();
    TestAllSymDiv<long double>();
    TestAllSymBandDiv<long double>();
#endif 

#ifdef TEST_INT
    TestBandMatrix<int>();
    TestSymMatrix<int>();
    TestSymBandMatrix<int>();
#endif 

#endif

    return 0;
}
#if 1
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
} catch (...) {
    std::cerr<<"Unknown exception thrown\n";
    std::cerr<<"Last successful test was "<<lastsuccess<<std::endl;
    return 1;
}
#else
catch (double) {}
#endif


