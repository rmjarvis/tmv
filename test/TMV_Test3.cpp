
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

int main() try {
    std::ofstream log("tmvtest3.log");
    tmv::WriteWarningsTo(&log);

    //showacc=true;
    //showdiv=true;
    //showtests=true;
    //showstartdone=true;

#if 1

#ifdef TEST_DOUBLE
    TestSmallVector<double>();
    TestSmallMatrix<double>();
    TestSmallMatrixA<double>();
    TestSmallMatrixB<double>();
    TestSmallMatrixDiv<double>();
    TestSmallMatrixDivA<double>();
    TestSmallMatrixDivB<double>();
#endif

#ifdef TEST_FLOAT
    TestSmallVector<float>();
    TestSmallMatrix<float>();
    TestSmallMatrixA<float>();
    TestSmallMatrixB<float>();
    TestSmallMatrixDiv<float>();
    TestSmallMatrixDivA<float>();
    TestSmallMatrixDivB<float>();
#endif

#ifdef TEST_LONGDOUBLE
    TestSmallVector<long double>();
    TestSmallMatrix<long double>();
    TestSmallMatrixA<long double>();
    TestSmallMatrixB<long double>();
    TestSmallMatrixDiv<long double>();
    TestSmallMatrixDivA<long double>();
    TestSmallMatrixDivB<long double>();
#endif 

#ifdef TEST_INT
    TestSmallVector<int>();
    TestSmallMatrix<int>();
#endif 

#endif

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

