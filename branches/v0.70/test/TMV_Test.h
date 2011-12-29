#ifndef TMV_TEST_H
#define TMV_TEST_H

#ifndef XTEST
#define XTEST 0
#endif

#include <iostream>
#include <cmath>
#include "tmv/TMV_Base.h"

#define EPS (10*tmv::TMV_Epsilon<T>())

#ifndef NO_TEST_DOUBLE
#define TEST_DOUBLE
#endif

#ifndef NO_TEST_FLOAT
#define TEST_FLOAT
#endif

#ifndef NO_TEST_COMPLEX
#define TEST_COMPLEX
#endif

extern bool showtests;
extern bool showacc;
extern bool showdiv;
extern bool donorm2;
extern bool showstartdone;
extern bool aliasok;
extern bool symoprod;
extern bool dontthrow;
extern std::string lastsuccess;

inline void PreAssert(std::string s)
{
    if (showtests) { 
        std::cout<<"Trying: "<<s;  
        std::cout.flush(); 
    } 
}

inline void DoAssert(bool x, std::string s)
{
    if (x) { 
        if (showtests) std::cout<<"  Passed"<<std::endl;
        lastsuccess = s; 
    } else { 
        if (showtests) std::cout<<"  Failed"<<std::endl;
        if (dontthrow) std::cout<<"Failed test: "<<s<<std::endl;  
        else {
#ifdef NOTHROW
            std::cerr<<"Error in test: "<<s<<std::endl; exit(1); 
#else
            throw tmv::Error("Error in test: ",s);  
#endif
        }
    } 
}

#define Assert(x,s) \
    do {  \
        PreAssert(s);  \
        DoAssert(x,s); \
    } while (false)

template <class M1, class M2, class T>
inline bool Equal(const M1& a, const M2& b, T eps)
{ 
    T normdiff = Norm(a-b);
    if (showtests) std::cout<<"  "<<normdiff<<" <=? "<<eps<<"  ";
    return normdiff <= eps; 
}
template <class X1, class X2, class T>
inline bool Equal2(const X1& a, const X2& b, T eps)
{
    T absdiff = tmv::TMV_ABS2(a-b);
    if (showtests) std::cout<<"  "<<absdiff<<" <=? "<<eps<<"  ";
    return absdiff <= eps;
}

template <class M1, class M2>
inline bool Equal(const M1& a, const M2& b, int )
{ return a == b; }
template <class X1, class X2>
inline bool Equal2(const X1& a, const X2& b, int )
{ return a == b; }

extern bool XXDEBUG1;
extern bool XXDEBUG2;
extern bool XXDEBUG3;
extern bool XXDEBUG4;
extern bool XXDEBUG5;
extern bool XXDEBUG6;
extern bool XXDEBUG7;
extern bool XXDEBUG8;
extern bool XXDEBUG9;

#ifdef TMV_TEST_DEPRECATED
#define subVector SubVector
#define subMatrix SubMatrix
#define subDiagMatrix SubDiagMatrix
#define subTriMatrix SubTriMatrix
#define subBandMatrix SubBandMatrix
#define subSymMatrix SubSymMatrix
#define subSymBandMatrix SubSymBandMatrix

#define reverse Reverse
#define view View
#define conjugate Conjugate
#define realPart Real
#define imagPart Imag
#define flatten Flatten

#define transpose Transpose
#define adjoint Adjoint
#define upperTri UpperTri
#define lowerTri LowerTri
#define offDiag OffDiag
#define upperBand UpperBand
#define lowerBand LowerBand
#define colPair ColPair
#define rowPair RowPair
#define colRange Cols
#define rowRange Rows

#define logDet LogDet
#define maxElement MaxElement
#define maxAbsElement MaxAbsElement
#define minElement MinElement
#define minAbsElement MinAbsElement
#define isSquare IsSquare
#define doNorm2 DoNorm2
#define doCondition DoCondition
#define condition Condition
#define writeCompact WriteCompact

#define inverse Inverse
#define makeInverse Inverse
#define makeInverseATA InverseATA

#define saveDiv SaveDiv
#define setDiv SetDiv
#define unsetDiv UnSetDiv
#define divideUsing DivideUsing
#define checkDecomp CheckDecomp
#define isTrans IsTrans

#define getS GetS
#define getV GetV
#define getU GetU
#define getL GetL
#define getP GetP
#define getD GetD
#define getQ GetQ
#define getR GetR

#define setZero Zero
#define clip Clip
#define setAllTo SetAllTo
#define addToAll AddToAll
#define conjugateSelf ConjugateSelf
#define makeBasis MakeBasis
#define swap Swap
#define permute Permute
#define reversePermute ReversePermute

#define transposeSelf TransposeSelf
#define swapRows SwapRows
#define swapCols SwapCols
#define swapRowsCols SwapRowsCols
#define permuteRows PermuteRows
#define permuteCols PermuteCols
#define permuteRowsCols PermuteRowsCols
#define reversePermuteRows ReversePermuteRows
#define reversePermuteCols ReversePermuteCols
#define reversePermuteRowsCols ReversePermuteRowsCols

#define sort Sort
#define Ascend ASCEND
#define Descend DESCEND
#define RealComp REAL_COMP
#define ImagComp IMAG_COMP
#define AbsComp ABS_COMP
#define ArgComp ARG_COMP

#define PolarDecompose Polar_Decompose

#define lud LUD
#define qrd QRD
#define qrpd QRPD
#define svd SVD
#define symsvd SymSVD
#define chd CHD

#endif

#endif
