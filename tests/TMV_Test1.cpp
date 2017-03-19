
#include <fstream>
#include "TMV_Test.h"
#include "TMV_Test_1.h"

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

int main() try {
    std::ofstream log("tmvtest1.log");
    tmv::WriteWarningsTo(&log);

    //showacc=true;
    //showdiv=true;
    //showtests=true;
    //showstartdone=true;

#if 1

#ifdef TEST_DOUBLE
    TestVector<double>();
    TestMatrix<double>();
    TestPermutation<double>();
    TestDiagMatrix<double>();
    TestDiagDiv<double>();
    TestTriMatrix<double>();
    TestTriDiv<double>();
    TestMatrixDiv<double>();
    TestMatrixDet<double>();
#endif // DOUBLE

#ifdef TEST_FLOAT
    TestVector<float>();
    TestMatrix<float>();
    TestPermutation<float>();
    TestDiagMatrix<float>();
    TestDiagDiv<float>();
    TestTriMatrix<float>();
    TestTriDiv<float>();
    TestMatrixDiv<float>();
    TestMatrixDet<float>();
#endif // FLOAT

#ifdef TEST_INT
    TestVector<int>();
    TestMatrix<int>();
    TestPermutation<int>();
    TestDiagMatrix<int>();
    TestTriMatrix<int>();
    TestMatrixDet<int>();
#endif  // INT

#ifdef TEST_LONGDOUBLE
    TestVector<long double>();
    TestMatrix<long double>();
    TestPermutation<long double>();
    TestDiagMatrix<long double>();
    TestDiagDiv<long double>();
    TestTriMatrix<long double>();
    TestTriDiv<long double>();
    TestMatrixDiv<long double>();
    TestMatrixDet<long double>();
#endif // LONGDOUBLE

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
}
catch (...) {
    std::cerr<<"Unknown exception thrown\n";
    std::cerr<<"Last successful test was "<<lastsuccess<<std::endl;
    return 1;
}
#else
catch (double) {}
#endif

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

