
#define TMV_EXTRA_DEBUG 

#include "TMV.h"
// Note: extra include file for BandMatrix
#include "TMV_Band.h"
// Also need to link with -ltmv_symband
#include <iostream>

int main() try 
{
    // Several ways to create and initialize band matrices:

    // Create with uninitialized values
    tmv::BandMatrix<double> m1(6,6,1,2);
    for(int i=0;i<m1.nrows();i++) 
        for(int j=0;j<m1.ncols();j++) 
            if (i<=j+m1.nlo() && j<=i+m1.nhi())
                m1(i,j) = 3.*i-j*j+7.; 
    std::cout<<"m1 =\n"<<m1;
    //! m1 =
    //! 6  6  
    //! ( 7  6  3  0  0  0 )
    //! ( 10  9  6  1  0  0 )
    //! ( 0  12  9  4  -3  0 )
    //! ( 0  0  12  7  0  -9 )
    //! ( 0  0  0  10  3  -6 )
    //! ( 0  0  0  0  6  -3 )

    // Create with all 2's.
    tmv::BandMatrix<double> m2(6,6,1,3,2.);
    std::cout<<"m2 =\n"<<m2;
    //! m2 =
    //! 6  6  
    //! ( 2  2  2  2  0  0 )
    //! ( 2  2  2  2  2  0 )
    //! ( 0  2  2  2  2  2 )
    //! ( 0  0  2  2  2  2 )
    //! ( 0  0  0  2  2  2 )
    //! ( 0  0  0  0  2  2 )

    // A BandMatrix can be non-square:
    tmv::BandMatrix<double> m3(6,8,1,3,2.);
    std::cout<<"m3 =\n"<<m3;
    //! m3 =
    //! 6  8  
    //! ( 2  2  2  2  0  0  0  0 )
    //! ( 2  2  2  2  2  0  0  0 )
    //! ( 0  2  2  2  2  2  0  0 )
    //! ( 0  0  2  2  2  2  2  0 )
    //! ( 0  0  0  2  2  2  2  2 )
    //! ( 0  0  0  0  2  2  2  2 )

    // Create from given elements:
    double mm[20] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
    tmv::BandMatrix<double,tmv::ColMajor> m4(6,6,2,1);
    std::copy(mm,mm+20,m4.colmajor_begin());
    std::cout<<"m4 (ColMajor) =\n"<<m4;
    //! m4 (ColMajor) =
    //! 6  6  
    //! ( 1  4  0  0  0  0 )
    //! ( 2  5  8  0  0  0 )
    //! ( 3  6  9  12  0  0 )
    //! ( 0  7  10  13  16  0 )
    //! ( 0  0  11  14  17  19 )
    //! ( 0  0  0  15  18  20 )
    tmv::BandMatrix<double,tmv::RowMajor> m5(6,6,2,1);
    std::copy(mm,mm+20,m5.rowmajor_begin());
    std::cout<<"m5 (RowMajor) =\n"<<m5;
    //! m5 (RowMajor) =
    //! 6  6  
    //! ( 1  2  0  0  0  0 )
    //! ( 3  4  5  0  0  0 )
    //! ( 6  7  8  9  0  0 )
    //! ( 0  10  11  12  13  0 )
    //! ( 0  0  14  15  16  17 )
    //! ( 0  0  0  18  19  20 )
    tmv::BandMatrix<double,tmv::DiagMajor> m6(6,6,2,1);
    std::copy(mm,mm+20,m6.diagmajor_begin());
    std::cout<<"m6 (DiagMajor) =\n"<<m6;
    //! m6 (DiagMajor) =
    //! 6  6  
    //! ( 10  16  0  0  0  0 )
    //! ( 5  11  17  0  0  0 )
    //! ( 1  6  12  18  0  0 )
    //! ( 0  2  7  13  19  0 )
    //! ( 0  0  3  8  14  20 )
    //! ( 0  0  0  4  9  15 )

    // Can make from the banded portion of a regular Matrix:
    tmv::Matrix<double> xm(6,6);
    for(int i=0;i<xm.nrows();i++) 
        for(int j=0;j<xm.ncols();j++) 
            xm(i,j) = 5.*i-j*j+3.; 
    tmv::BandMatrix<double> m7(xm,3,2);
    std::cout<<"m7 =\n"<<m7;
    //! m7 =
    //! 6  6  
    //! ( 3  2  -1  0  0  0 )
    //! ( 8  7  4  -1  0  0 )
    //! ( 13  12  9  4  -3  0 )
    //! ( 18  17  14  9  2  -7 )
    //! ( 0  22  19  14  7  -2 )
    //! ( 0  0  24  19  12  3 )
    // Or from a wider BandMatrix:
    tmv::BandMatrix<double> m8(m7,3,0);
    std::cout<<"m8 =\n"<<m8;
    //! m8 =
    //! 6  6  
    //! ( 3  0  0  0  0  0 )
    //! ( 8  7  0  0  0  0 )
    //! ( 13  12  9  0  0  0 )
    //! ( 18  17  14  9  0  0 )
    //! ( 0  22  19  14  7  0 )
    //! ( 0  0  24  19  12  3 )

    // Shortcuts to Bi- and Tri-diagonal matrices:
    tmv::Vector<double> v1(5,1.);
    tmv::Vector<double> v2(6,2.);
    tmv::Vector<double> v3(5,3.);
    tmv::BandMatrix<double> m9 = LowerBiDiagMatrix(v1,v2);
    tmv::BandMatrix<double> m10 = UpperBiDiagMatrix(v2,v3);
    tmv::BandMatrix<double> m11 = TriDiagMatrix(v1,v2,v3);
    std::cout<<"LowerBiDiagMatrix(v1,v2) =\n"<<m9;
    //! LowerBiDiagMatrix(v1,v2) =
    //! 6  6  
    //! ( 2  0  0  0  0  0 )
    //! ( 1  2  0  0  0  0 )
    //! ( 0  1  2  0  0  0 )
    //! ( 0  0  1  2  0  0 )
    //! ( 0  0  0  1  2  0 )
    //! ( 0  0  0  0  1  2 )
    std::cout<<"UpperBiDiagMatrix(v2,v3) =\n"<<m10;
    //! UpperBiDiagMatrix(v2,v3) =
    //! 6  6  
    //! ( 2  3  0  0  0  0 )
    //! ( 0  2  3  0  0  0 )
    //! ( 0  0  2  3  0  0 )
    //! ( 0  0  0  2  3  0 )
    //! ( 0  0  0  0  2  3 )
    //! ( 0  0  0  0  0  2 )
    std::cout<<"TriDiagMatrix(v1,v2,v3) =\n"<<m11;
    //! TriDiagMatrix(v1,v2,v3) =
    //! 6  6  
    //! ( 2  3  0  0  0  0 )
    //! ( 1  2  3  0  0  0 )
    //! ( 0  1  2  3  0  0 )
    //! ( 0  0  1  2  3  0 )
    //! ( 0  0  0  1  2  3 )
    //! ( 0  0  0  0  1  2 )


    // Norms, etc. 

    std::cout<<"Norm1(m1) = "<<Norm1(m1)<<std::endl;
    //! Norm1(m1) = 30
    std::cout<<"Norm2(m1) = "<<Norm2(m1)<<std::endl;
    //! Norm2(m1) = 24.0314
    std::cout<<"NormInf(m1) = "<<NormInf(m1)<<std::endl;
    //! NormInf(m1) = 28
    std::cout<<"NormF(m1) = "<<NormF(m1)<<" = "<<Norm(m1)<<std::endl;
    //! NormF(m1) = 32.0312 = 32.0312
    std::cout<<"MaxAbsElement(m1) = "<<MaxAbsElement(m1)<<std::endl;
    //! MaxAbsElement(m1) = 12
    std::cout<<"Trace(m1) = "<<Trace(m1)<<std::endl;
    //! Trace(m1) = 32
    std::cout<<"Det(m1) = "<<Det(m1)<<std::endl;
    //! Det(m1) = 67635


    // Views:

    std::cout<<"m1 =\n"<<m1;
    //! m1 =
    //! 6  6  
    //! ( 7  6  3  0  0  0 )
    //! ( 10  9  6  1  0  0 )
    //! ( 0  12  9  4  -3  0 )
    //! ( 0  0  12  7  0  -9 )
    //! ( 0  0  0  10  3  -6 )
    //! ( 0  0  0  0  6  -3 )
    std::cout<<"m1.diag() = "<<m1.diag()<<std::endl;
    //! m1.diag() = 6  ( 7  9  9  7  3  -3 )
    std::cout<<"m1.diag(1) = "<<m1.diag(1)<<std::endl;
    //! m1.diag(1) = 5  ( 6  6  4  0  -6 )
    std::cout<<"m1.diag(-1) = "<<m1.diag(-1)<<std::endl;
    //! m1.diag(-1) = 5  ( 10  12  12  10  6 )
    std::cout<<"m1.subBandMatrix(0,3,0,3,1,1) =\n"<<
        m1.subBandMatrix(0,3,0,3,1,1);
    //! m1.subBandMatrix(0,3,0,3,1,1) =
    //! 3  3  
    //! ( 7  6  0 )
    //! ( 10  9  6 )
    //! ( 0  12  9 )
    std::cout<<"m1.transpose() =\n"<<m1.transpose();
    //! m1.transpose() =
    //! 6  6  
    //! ( 7  10  0  0  0  0 )
    //! ( 6  9  12  0  0  0 )
    //! ( 3  6  9  12  0  0 )
    //! ( 0  1  4  7  10  0 )
    //! ( 0  0  -3  0  3  6 )
    //! ( 0  0  0  -9  -6  -3 )

    // rowRange, colRange shrink both dimensions of the matrix to include only
    // the portions that are in those rows or columns:
    std::cout<<"m1.rowRange(0,4) =\n"<<m1.rowRange(0,4);
    //! m1.rowRange(0,4) =
    //! 4  6  
    //! ( 7  6  3  0  0  0 )
    //! ( 10  9  6  1  0  0 )
    //! ( 0  12  9  4  -3  0 )
    //! ( 0  0  12  7  0  -9 )
    std::cout<<"m1.colRange(1,4) =\n"<<m1.colRange(1,4);
    //! m1.colRange(1,4) =
    //! 5  3  
    //! ( 6  3  0 )
    //! ( 9  6  1 )
    //! ( 12  9  4 )
    //! ( 0  12  7 )
    //! ( 0  0  10 )
    std::cout<<"m1.diagRange(0,2) =\n"<<m1.diagRange(0,2);
    //! m1.diagRange(0,2) =
    //! 6  6  
    //! ( 7  6  0  0  0  0 )
    //! ( 0  9  6  0  0  0 )
    //! ( 0  0  9  4  0  0 )
    //! ( 0  0  0  7  0  0 )
    //! ( 0  0  0  0  3  -6 )
    //! ( 0  0  0  0  0  -3 )
    std::cout<<"m1.diagRange(-1,1) =\n"<<m1.diagRange(-1,1);
    //! m1.diagRange(-1,1) =
    //! 6  6  
    //! ( 7  0  0  0  0  0 )
    //! ( 10  9  0  0  0  0 )
    //! ( 0  12  9  0  0  0 )
    //! ( 0  0  12  7  0  0 )
    //! ( 0  0  0  10  3  0 )
    //! ( 0  0  0  0  6  -3 )


    // Fortran Indexing:

    tmv::BandMatrix<double,tmv::FortranStyle> fm1 = m1;
    std::cout<<"fm1 = m1 =\n"<<fm1;
    //! fm1 = m1 =
    //! 6  6  
    //! ( 7  6  3  0  0  0 )
    //! ( 10  9  6  1  0  0 )
    //! ( 0  12  9  4  -3  0 )
    //! ( 0  0  12  7  0  -9 )
    //! ( 0  0  0  10  3  -6 )
    //! ( 0  0  0  0  6  -3 )
    std::cout<<"fm1(1,1) = "<<fm1(1,1)<<std::endl;
    //! fm1(1,1) = 7
    std::cout<<"fm1(4,3) = "<<fm1(4,3)<<std::endl;
    //! fm1(4,3) = 12
    std::cout<<"fm1.subBandMatrix(1,3,1,3,1,1) =\n"<<
        fm1.subBandMatrix(1,3,1,3,1,1);
    //! fm1.subBandMatrix(1,3,1,3,1,1) =
    //! 3  3  
    //! ( 7  6  0 )
    //! ( 10  9  6 )
    //! ( 0  12  9 )
    std::cout<<"fm1.rowRange(1,4) =\n"<<fm1.rowRange(1,4);
    //! fm1.rowRange(1,4) =
    //! 4  6  
    //! ( 7  6  3  0  0  0 )
    //! ( 10  9  6  1  0  0 )
    //! ( 0  12  9  4  -3  0 )
    //! ( 0  0  12  7  0  -9 )
    std::cout<<"fm1.colRange(2,4) =\n"<<fm1.colRange(2,4);
    //! fm1.colRange(2,4) =
    //! 5  3  
    //! ( 6  3  0 )
    //! ( 9  6  1 )
    //! ( 12  9  4 )
    //! ( 0  12  7 )
    //! ( 0  0  10 )
    std::cout<<"fm1.diagRange(0,1) =\n"<<fm1.diagRange(0,1);
    //! fm1.diagRange(0,1) =
    //! 6  6  
    //! ( 7  6  0  0  0  0 )
    //! ( 0  9  6  0  0  0 )
    //! ( 0  0  9  4  0  0 )
    //! ( 0  0  0  7  0  0 )
    //! ( 0  0  0  0  3  -6 )
    //! ( 0  0  0  0  0  -3 )
    std::cout<<"fm1.diagRange(-1,0) =\n"<<fm1.diagRange(-1,0);
    //! fm1.diagRange(-1,0) =
    //! 6  6  
    //! ( 7  0  0  0  0  0 )
    //! ( 10  9  0  0  0  0 )
    //! ( 0  12  9  0  0  0 )
    //! ( 0  0  12  7  0  0 )
    //! ( 0  0  0  10  3  0 )
    //! ( 0  0  0  0  6  -3 )


    // Matrix arithmetic:

    tmv::BandMatrix<double> m1pm2 = m1 + m2;
    std::cout<<"m1 + m2 =\n"<<m1pm2;
    //! m1 + m2 =
    //! 6  6  
    //! ( 9  8  5  2  0  0 )
    //! ( 12  11  8  3  2  0 )
    //! ( 0  14  11  6  -1  2 )
    //! ( 0  0  14  9  2  -7 )
    //! ( 0  0  0  12  5  -4 )
    //! ( 0  0  0  0  8  -1 )
    // Works correctly even if matrices are stored in different order:
    tmv::BandMatrix<double> m5pm6 = m5 + m6; 
    std::cout<<"m5 + m6 =\n"<<m5pm6;
    //! m5 + m6 =
    //! 6  6  
    //! ( 11  18  0  0  0  0 )
    //! ( 8  15  22  0  0  0 )
    //! ( 7  13  20  27  0  0 )
    //! ( 0  12  18  25  32  0 )
    //! ( 0  0  17  23  30  37 )
    //! ( 0  0  0  22  28  35 )
    // Also expands the number of off-diagonals appropriately as needed:
    tmv::BandMatrix<double> m2pm4 = m2 + m4; 
    std::cout<<"m2 + m4 =\n"<<m2pm4;
    //! m2 + m4 =
    //! 6  6  
    //! ( 3  6  2  2  0  0 )
    //! ( 4  7  10  2  2  0 )
    //! ( 3  8  11  14  2  2 )
    //! ( 0  7  12  15  18  2 )
    //! ( 0  0  11  16  19  21 )
    //! ( 0  0  0  15  20  22 )

    m1 *= 2.;
    std::cout<<"m1 *= 2 =\n"<<m1;
    //! m1 *= 2 =
    //! 6  6  
    //! ( 14  12  6  0  0  0 )
    //! ( 20  18  12  2  0  0 )
    //! ( 0  24  18  8  -6  0 )
    //! ( 0  0  24  14  0  -18 )
    //! ( 0  0  0  20  6  -12 )
    //! ( 0  0  0  0  12  -6 )

    m2 += m1;
    std::cout<<"m2 += m1 =\n"<<m2;
    //! m2 += m1 =
    //! 6  6  
    //! ( 16  14  8  2  0  0 )
    //! ( 22  20  14  4  2  0 )
    //! ( 0  26  20  10  -4  2 )
    //! ( 0  0  26  16  2  -16 )
    //! ( 0  0  0  22  8  -10 )
    //! ( 0  0  0  0  14  -4 )

    tmv::Vector<double> v = xm.col(0);
    std::cout<<"v = "<<v<<std::endl;
    //! v = 6  ( 3  8  13  18  23  28 )
    std::cout<<"m1 * v = "<<m1*v<<std::endl;
    //! m1 * v = 6  ( 216  396  432  60  162  108 )
    std::cout<<"v * m1 = "<<v*m1<<std::endl;
    //! v * m1 = 6  ( 202  492  780  832  396  -768 )

    // Matrix * matrix product also expands bands appropriately:
    tmv::BandMatrix<double> m1m2 = m1 * m2; 
    std::cout<<"m1 * m2 =\n"<<m1m2;
    //! m1 * m2 =
    //! 6  6  
    //! ( 488  592  400  136  0  12 )
    //! ( 716  952  704  264  -8  -8 )
    //! ( 528  948  904  272  -56  -32 )
    //! ( 0  624  844  464  -320  -104 )
    //! ( 0  0  520  452  -80  -332 )
    //! ( 0  0  0  264  12  -96 )

    // Can mix BandMatrix with other kinds of matrices:
    std::cout<<"xm * m1 =\n"<<xm*m1;
    //! xm * m1 =
    //! 6  6  
    //! ( 82  48  -120  -348  -336  396 )
    //! ( 252  318  180  -128  -276  216 )
    //! ( 422  588  480  92  -216  36 )
    //! ( 592  858  780  312  -156  -144 )
    //! ( 762  1128  1080  532  -96  -324 )
    //! ( 932  1398  1380  752  -36  -504 )
    tmv::UpperTriMatrix<double> um(xm);
    std::cout<<"um + m1 =\n"<<um+m1;
    //! um + m1 =
    //! 6  6  
    //! ( 17  14  5  -6  -13  -22 )
    //! ( 20  25  16  1  -8  -17 )
    //! ( 0  24  27  12  -9  -12 )
    //! ( 0  0  24  23  2  -25 )
    //! ( 0  0  0  20  13  -14 )
    //! ( 0  0  0  0  12  -3 )
    tmv::LowerTriMatrix<double> lm(xm);
    lm *= m8;
    std::cout<<"lm *= m8 =\n"<<lm;
    //! lm *= m8 =
    //! 6  6  
    //! ( 9  0  0  0  0  0 )
    //! ( 80  49  0  0  0  0 )
    //! ( 252  192  81  0  0  0 )
    //! ( 534  440  252  81  0  0 )
    //! ( 744  774  500  224  49  0 )
    //! ( 954  1064  782  396  120  9 )
    tmv::DiagMatrix<double> dm(xm);
    m1 *= dm;
    std::cout<<"m1 *= dm =\n"<<m1;
    //! m1 *= dm =
    //! 6  6  
    //! ( 42  84  54  0  0  0 )
    //! ( 60  126  108  18  0  0 )
    //! ( 0  168  162  72  -42  0 )
    //! ( 0  0  216  126  0  -54 )
    //! ( 0  0  0  180  42  -36 )
    //! ( 0  0  0  0  84  -18 )
    return 0;
} catch (tmv::Error& e) {
    std::cerr<<e<<std::endl;
    return 1;
}
