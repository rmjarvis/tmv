
#define TMV_EXTRA_DEBUG // To get the 888's below.

#include "TMV.h"
// Note: extra include file for SymMatrix
#include "TMV_Sym.h"
// Also need to link with -ltmv_symband
#include <iostream>

int main() try 
{
    // Several ways to create and initialize matrices:

    // Create with uninitialized values
    tmv::SymMatrix<double> m1(5); 
    // In debug mode, the elements are all set to 888 to make it easier
    // to find errors related to not correctly initializing the matrix.
    std::cout<<"m1 =\n"<<m1;
    //! m1 =
    //! 5  5  
    //! ( 888  888  888  888  888 )
    //! ( 888  888  888  888  888 )
    //! ( 888  888  888  888  888 )
    //! ( 888  888  888  888  888 )
    //! ( 888  888  888  888  888 )

    // Initialize with element access:
    for(int i=0;i<m1.nrows();i++) 
        for(int j=0;j<m1.ncols();j++) 
            if (i>=j) m1(i,j) = 2.*i-j*j+2.; 
    std::cout<<"m1 =\n"<<m1;
    //! m1 =
    //! 5  5  
    //! ( 2  4  6  8  10 )
    //! ( 4  3  5  7  9 )
    //! ( 6  5  2  4  6 )
    //! ( 8  7  4  -1  1 )
    //! ( 10  9  6  1  -6 )

    tmv::SymMatrix<double> m2(5,2.); // Create with all 2's.
    std::cout<<"m2 =\n"<<m2;
    //! m2 =
    //! 5  5  
    //! ( 2  2  2  2  2 )
    //! ( 2  2  2  2  2 )
    //! ( 2  2  2  2  2 )
    //! ( 2  2  2  2  2 )
    //! ( 2  2  2  2  2 )

    // Initialize from given elements:
    double mm[15] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    tmv::SymMatrix<double,tmv::Upper|tmv::ColMajor> m3(5); 
    // Use colmajor iterator of the upperTri portion:
    std::copy(mm,mm+15,m3.upperTri().colmajor_begin());
    std::cout<<"m3 (Upper, ColMajor) =\n"<<m3;
    //! m3 (Upper, ColMajor) =
    //! 5  5  
    //! ( 1  2  4  7  11 )
    //! ( 2  3  5  8  12 )
    //! ( 4  5  6  9  13 )
    //! ( 7  8  9  10  14 )
    //! ( 11  12  13  14  15 )
    
    // The iteration does not have to match the internal storage of 
    // the SymMatrix.  It will just be slightly less efficient.
    tmv::SymMatrix<double,tmv::Lower|tmv::RowMajor> m4(5);
    std::copy(mm,mm+15,m4.upperTri().rowmajor_begin());
    std::cout<<"m4 (Lower,RowMajor) initialized with (Upper,RowMajor) =\n"<<m4;
    //! m4 (Lower,RowMajor) initialized with (Upper,RowMajor) =
    //! 5  5  
    //! ( 1  2  3  4  5 )
    //! ( 2  6  7  8  9 )
    //! ( 3  7  10  11  12 )
    //! ( 4  8  11  13  14 )
    //! ( 5  9  12  14  15 )
    
    // Initialiion from comma separated list is not possible, 
    // but you can initialize either triangle:
    tmv::SymMatrix<double> m5(5); // Default is Lower, ColMajor
    m5.lowerTri() << 
        4, 
        8, 10,
        1, -4, -7,
        9,  3,  3, 12, 
        5,  0,  1,  8, 9;
    std::cout<<"m5 =\n"<<m5;
    //! m5 =
    //! 5  5  
    //! ( 4  8  1  9  5 )
    //! ( 8  10  -4  3  0 )
    //! ( 1  -4  -7  3  1 )
    //! ( 9  3  3  12  8 )
    //! ( 5  0  1  8  9 )
    
    tmv::SymMatrix<double> m6(5);
    m6.upperTri() <<
        4,  8,  1,  9,  5,
           10, -4,  3,  0,
               -7,  3,  1,
                   12,  8,
                        9;

    std::cout<<"m6 =\n"<<m6;
    //! m6 =
    //! 5  5  
    //! ( 4  8  1  9  5 )
    //! ( 8  10  -4  3  0 )
    //! ( 1  -4  -7  3  1 )
    //! ( 9  3  3  12  8 )
    //! ( 5  0  1  8  9 )

    // Make from the corresponding portion of a regular Matrix
    // Note - it copies from the correct triangle regardless of 
    // the storage order of the two matrices.
    tmv::Matrix<double,tmv::RowMajor> xm(5,5);
    xm <<
        2, 5, 1, 9, 8,
        1, 5, 7, 2, 0,
        3, 9, 6, 8, 4,
        4, 2, 1, 9, 0,
        9, 8, 3, 7, 5;
    tmv::SymMatrix<double,tmv::Lower|tmv::ColMajor> m7(xm);
    std::cout<<"m7 =\n"<<m7;
    //! m7 =
    //! 5  5  
    //! ( 2  1  3  4  9 )
    //! ( 1  5  9  2  8 )
    //! ( 3  9  6  1  3 )
    //! ( 4  2  1  9  7 )
    //! ( 9  8  3  7  5 )


    // Norms, etc. 

    std::cout<<"Norm1(m1) = "<<Norm1(m1)<<std::endl;
    //! Norm1(m1) = 32
    std::cout<<"Norm2(m1) = "<<Norm2(m1)<<std::endl;
    //! Norm2(m1) = 24.4968
    std::cout<<"NormInf(m1) = "<<NormInf(m1)<<std::endl;
    //! NormInf(m1) = 32
    std::cout<<"NormF(m1) = "<<NormF(m1);
    std::cout<<" = "<<Norm(m1)<<std::endl;
    //! NormF(m1) = 30.0333 = 30.0333
    std::cout<<"MaxAbsElement(m1) = "<<MaxAbsElement(m1)<<std::endl;
    //! MaxAbsElement(m1) = 10
    std::cout<<"Trace(m1) = "<<Trace(m1)<<std::endl;
    //! Trace(m1) = 0
    std::cout<<"Det(m1) = "<<Det(m1)<<std::endl;
    //! Det(m1) = 4866


    // Views:

    std::cout<<"m1 =\n"<<m1;
    //! m1 =
    //! 5  5  
    //! ( 2  4  6  8  10 )
    //! ( 4  3  5  7  9 )
    //! ( 6  5  2  4  6 )
    //! ( 8  7  4  -1  1 )
    //! ( 10  9  6  1  -6 )
    std::cout<<"m1.diag() = "<<m1.diag()<<std::endl;
    //! m1.diag() = 5  ( 2  3  2  -1  -6 )
    std::cout<<"m1.diag(1) = "<<m1.diag(1)<<std::endl;
    //! m1.diag(1) = 4  ( 4  5  4  1 )
    std::cout<<"m1.diag(-1) = "<<m1.diag(-1)<<std::endl;
    //! m1.diag(-1) = 4  ( 4  5  4  1 )
    std::cout<<"m1.subSymMatrix(0,3) =\n"<<m1.subSymMatrix(0,3);
    //! m1.subSymMatrix(0,3) =
    //! 3  3  
    //! ( 2  4  6 )
    //! ( 4  3  5 )
    //! ( 6  5  2 )
    std::cout<<"m1.upperTri() =\n"<<m1.upperTri();
    //! m1.upperTri() =
    //! 5  5  
    //! ( 2  4  6  8  10 )
    //! ( 0  3  5  7  9 )
    //! ( 0  0  2  4  6 )
    //! ( 0  0  0  -1  1 )
    //! ( 0  0  0  0  -6 )
    std::cout<<"m1.lowerTri(tmv::UnitDiag) =\n"<<m1.lowerTri(tmv::UnitDiag);
    //! m1.lowerTri(tmv::UnitDiag) =
    //! 5  5  
    //! ( 1  0  0  0  0 )
    //! ( 4  1  0  0  0 )
    //! ( 6  5  1  0  0 )
    //! ( 8  7  4  1  0 )
    //! ( 10  9  6  1  1 )


    // Fortran-style indexing:

    tmv::SymMatrix<double,tmv::FortranStyle> fm1 = m1;
    std::cout<<"fm1 = m1 =\n"<<fm1;
    //! fm1 = m1 =
    //! 5  5  
    //! ( 2  4  6  8  10 )
    //! ( 4  3  5  7  9 )
    //! ( 6  5  2  4  6 )
    //! ( 8  7  4  -1  1 )
    //! ( 10  9  6  1  -6 )
    std::cout<<"fm1(1,1) = "<<fm1(1,1)<<std::endl;
    //! fm1(1,1) = 2
    std::cout<<"fm1(4,3) = "<<fm1(4,3)<<std::endl;
    //! fm1(4,3) = 4
    std::cout<<"fm1.subSymMatrix(1,3) =\n"<<fm1.subSymMatrix(1,3);
    //! fm1.subSymMatrix(1,3) =
    //! 3  3  
    //! ( 2  4  6 )
    //! ( 4  3  5 )
    //! ( 6  5  2 )


    // Matrix arithmetic:

    tmv::SymMatrix<double> m1pm3 = m1 + m3;  
    std::cout<<"m1 + m3 =\n"<<m1pm3;
    //! m1 + m3 =
    //! 5  5  
    //! ( 3  6  10  15  21 )
    //! ( 6  6  10  15  21 )
    //! ( 10  10  8  13  19 )
    //! ( 15  15  13  9  15 )
    //! ( 21  21  19  15  9 )
    
    // Works correctly even if matrices are stored in different order:
    tmv::SymMatrix<double> m3pm4 = m3 + m4; 
    std::cout<<"m3 + m4 =\n"<<m3pm4;
    //! m3 + m4 =
    //! 5  5  
    //! ( 2  4  7  11  16 )
    //! ( 4  9  12  16  21 )
    //! ( 7  12  16  20  25 )
    //! ( 11  16  20  23  28 )
    //! ( 16  21  25  28  30 )
    
    tmv::SymMatrix<double> m4pm5 = m4 + m5; 
    std::cout<<"m4 + m5 =\n"<<m4pm5;
    //! m4 + m5 =
    //! 5  5  
    //! ( 5  10  4  13  10 )
    //! ( 10  16  3  11  9 )
    //! ( 4  3  3  14  13 )
    //! ( 13  11  14  25  22 )
    //! ( 10  9  13  22  24 )

    m1 *= 2.;
    std::cout<<"m1 *= 2 =\n"<<m1;
    //! m1 *= 2 =
    //! 5  5  
    //! ( 4  8  12  16  20 )
    //! ( 8  6  10  14  18 )
    //! ( 12  10  4  8  12 )
    //! ( 16  14  8  -2  2 )
    //! ( 20  18  12  2  -12 )

    m1 += m4;
    std::cout<<"m1 += m4 =\n"<<m1;
    //! m1 += m4 =
    //! 5  5  
    //! ( 5  10  15  20  25 )
    //! ( 10  12  17  22  27 )
    //! ( 15  17  14  19  24 )
    //! ( 20  22  19  11  16 )
    //! ( 25  27  24  16  3 )

    // Vector outer product
    tmv::Vector<double> v = xm.col(0);
    tmv::SymMatrix<double> vv = v^v;
    std::cout<<"v = "<<v<<std::endl;
    //! v = 5  ( 2  1  3  4  9 )
    std::cout<<"v^v =\n"<<vv;
    //! v^v =
    //! 5  5  
    //! ( 4  2  6  8  18 )
    //! ( 2  1  3  4  9 )
    //! ( 6  3  9  12  27 )
    //! ( 8  4  12  16  36 )
    //! ( 18  9  27  36  81 )

    // SymMatrix * Vector product
    std::cout<<"m1 * v = "<<m1*v<<std::endl;
    //! m1 * v = 5  ( 370  414  381  307  240 )
    std::cout<<"v * m1 = "<<v*m1<<std::endl;
    // Note: This is only the same answer because m1 is symmetric.
    //! v * m1 = 5  ( 370  414  381  307  240 )

    // SymMatrix * Matrix product
    tmv::Matrix<double> m1xm = m1 * xm; 
    std::cout<<"m1 * xm =\n"<<m1xm;
    //! m1 * xm =
    //! 5  5  
    //! ( 370  450  260  540  225 )
    //! ( 414  523  299  637  283 )
    //! ( 381  516  309  620  296 )
    //! ( 307  531  347  587  316 )
    //! ( 240  532  383  636  311 )
    
    // Note: the product of two symmetrix matrices is not symmetric!:
    tmv::Matrix<double> m1m5 = m1 * m5; 
    std::cout<<"m1 * m5 =\n"<<m1m5;
    //! m1 * m5 =
    //! 5  5  
    //! ( 420  140  -55  560  425 )
    //! ( 486  198  -64  657  486 )
    //! ( 501  291  -70  648  457 )
    //! ( 454  337  -152  563  351 )
    //! ( 499  422  -200  594  304 )


    // Complex matrices:

    tmv::SymMatrix<std::complex<double> > cm1 = m1 * std::complex<double>(1,2);
    std::cout<<"cm1 = m1 * (1+2i) =\n"<<cm1;
    //! cm1 = m1 * (1+2i) =
    //! 5  5  
    //! ( (5,10)  (10,20)  (15,30)  (20,40)  (25,50) )
    //! ( (10,20)  (12,24)  (17,34)  (22,44)  (27,54) )
    //! ( (15,30)  (17,34)  (14,28)  (19,38)  (24,48) )
    //! ( (20,40)  (22,44)  (19,38)  (11,22)  (16,32) )
    //! ( (25,50)  (27,54)  (24,48)  (16,32)  (3,6) )
    std::cout<<"cm1.conjugate() =\n"<<cm1.conjugate();
    //! cm1.conjugate() =
    //! 5  5  
    //! ( (5,-10)  (10,-20)  (15,-30)  (20,-40)  (25,-50) )
    //! ( (10,-20)  (12,-24)  (17,-34)  (22,-44)  (27,-54) )
    //! ( (15,-30)  (17,-34)  (14,-28)  (19,-38)  (24,-48) )
    //! ( (20,-40)  (22,-44)  (19,-38)  (11,-22)  (16,-32) )
    //! ( (25,-50)  (27,-54)  (24,-48)  (16,-32)  (3,-6) )
    std::cout<<"cm1.transpose() =\n"<<cm1.transpose();
    //! cm1.transpose() =
    //! 5  5  
    //! ( (5,10)  (10,20)  (15,30)  (20,40)  (25,50) )
    //! ( (10,20)  (12,24)  (17,34)  (22,44)  (27,54) )
    //! ( (15,30)  (17,34)  (14,28)  (19,38)  (24,48) )
    //! ( (20,40)  (22,44)  (19,38)  (11,22)  (16,32) )
    //! ( (25,50)  (27,54)  (24,48)  (16,32)  (3,6) )
    std::cout<<"cm1.adjoint() =\n"<<cm1.adjoint();
    //! cm1.adjoint() =
    //! 5  5  
    //! ( (5,-10)  (10,-20)  (15,-30)  (20,-40)  (25,-50) )
    //! ( (10,-20)  (12,-24)  (17,-34)  (22,-44)  (27,-54) )
    //! ( (15,-30)  (17,-34)  (14,-28)  (19,-38)  (24,-48) )
    //! ( (20,-40)  (22,-44)  (19,-38)  (11,-22)  (16,-32) )
    //! ( (25,-50)  (27,-54)  (24,-48)  (16,-32)  (3,-6) )
    std::cout<<"cm1.realPart() =\n"<<cm1.realPart();
    //! cm1.realPart() =
    //! 5  5  
    //! ( 5  10  15  20  25 )
    //! ( 10  12  17  22  27 )
    //! ( 15  17  14  19  24 )
    //! ( 20  22  19  11  16 )
    //! ( 25  27  24  16  3 )
    std::cout<<"cm1.imagPart() =\n"<<cm1.imagPart();
    //! cm1.imagPart() =
    //! 5  5  
    //! ( 10  20  30  40  50 )
    //! ( 20  24  34  44  54 )
    //! ( 30  34  28  38  48 )
    //! ( 40  44  38  22  32 )
    //! ( 50  54  48  32  6 )
    std::cout<<"Norm(cm1) = "<<Norm(cm1)<<std::endl;
    //! Norm(cm1) = 207.183

    tmv::Vector<std::complex<double> > cv = v * std::complex<double>(-2,1);
    cm1 += cv^cv;
    std::cout<<"cm1 += cv^cv =\n"<<cm1;
    //! cm1 += cv^cv =
    //! 5  5  
    //! ( (17,-6)  (16,12)  (33,6)  (44,8)  (79,-22) )
    //! ( (16,12)  (15,20)  (26,22)  (34,28)  (54,18) )
    //! ( (33,6)  (26,22)  (41,-8)  (55,-10)  (105,-60) )
    //! ( (44,8)  (34,28)  (55,-10)  (59,-42)  (124,-112) )
    //! ( (79,-22)  (54,18)  (105,-60)  (124,-112)  (246,-318) )
    cm1 += v^v;
    std::cout<<"cm1 += v^v =\n"<<cm1;
    //! cm1 += v^v =
    //! 5  5  
    //! ( (21,-6)  (18,12)  (39,6)  (52,8)  (97,-22) )
    //! ( (18,12)  (16,20)  (29,22)  (38,28)  (63,18) )
    //! ( (39,6)  (29,22)  (50,-8)  (67,-10)  (132,-60) )
    //! ( (52,8)  (38,28)  (67,-10)  (75,-42)  (160,-112) )
    //! ( (97,-22)  (63,18)  (132,-60)  (160,-112)  (327,-318) )

    tmv::Matrix<std::complex<double> > cx(5,2);
    cx.col(0) = v;
    cx.col(1) = cv;
    cm1 -= cx*cx.transpose();
    std::cout<<"cm1 -= cx*cx.transpose() =\n"<<cm1;
    //! cm1 -= cx*cx.transpose() =
    //! 5  5  
    //! ( (5,10)  (10,20)  (15,30)  (20,40)  (25,50) )
    //! ( (10,20)  (12,24)  (17,34)  (22,44)  (27,54) )
    //! ( (15,30)  (17,34)  (14,28)  (19,38)  (24,48) )
    //! ( (20,40)  (22,44)  (19,38)  (11,22)  (16,32) )
    //! ( (25,50)  (27,54)  (24,48)  (16,32)  (3,6) )

    tmv::HermMatrix<std::complex<double> > cm2 = m1;
    cm2.upperTri().offDiag() *= std::complex<double>(1,2);
    std::cout<<"cm2 =\n"<<cm2;
    //! cm2 =
    //! 5  5  
    //! ( (5,-0)  (10,20)  (15,30)  (20,40)  (25,50) )
    //! ( (10,-20)  (12,-0)  (17,34)  (22,44)  (27,54) )
    //! ( (15,-30)  (17,-34)  (14,-0)  (19,38)  (24,48) )
    //! ( (20,-40)  (22,-44)  (19,-38)  (11,-0)  (16,32) )
    //! ( (25,-50)  (27,-54)  (24,-48)  (16,-32)  (3,-0) )
    std::cout<<"cm2.conjugate() =\n"<<cm2.conjugate();
    //! cm2.conjugate() =
    //! 5  5  
    //! ( (5,0)  (10,-20)  (15,-30)  (20,-40)  (25,-50) )
    //! ( (10,20)  (12,0)  (17,-34)  (22,-44)  (27,-54) )
    //! ( (15,30)  (17,34)  (14,0)  (19,-38)  (24,-48) )
    //! ( (20,40)  (22,44)  (19,38)  (11,0)  (16,-32) )
    //! ( (25,50)  (27,54)  (24,48)  (16,32)  (3,0) )
    std::cout<<"cm2.transpose() =\n"<<cm2.transpose();
    //! cm2.transpose() =
    //! 5  5  
    //! ( (5,-0)  (10,-20)  (15,-30)  (20,-40)  (25,-50) )
    //! ( (10,20)  (12,-0)  (17,-34)  (22,-44)  (27,-54) )
    //! ( (15,30)  (17,34)  (14,-0)  (19,-38)  (24,-48) )
    //! ( (20,40)  (22,44)  (19,38)  (11,-0)  (16,-32) )
    //! ( (25,50)  (27,54)  (24,48)  (16,32)  (3,-0) )
    std::cout<<"cm2.adjoint() =\n"<<cm2.adjoint();
    //! cm2.adjoint() =
    //! 5  5  
    //! ( (5,0)  (10,20)  (15,30)  (20,40)  (25,50) )
    //! ( (10,-20)  (12,0)  (17,34)  (22,44)  (27,54) )
    //! ( (15,-30)  (17,-34)  (14,0)  (19,38)  (24,48) )
    //! ( (20,-40)  (22,-44)  (19,-38)  (11,0)  (16,32) )
    //! ( (25,-50)  (27,-54)  (24,-48)  (16,-32)  (3,0) )
    std::cout<<"cm2.realPart() =\n"<<cm2.realPart();
    //! cm2.realPart() =
    //! 5  5  
    //! ( 5  10  15  20  25 )
    //! ( 10  12  17  22  27 )
    //! ( 15  17  14  19  24 )
    //! ( 20  22  19  11  16 )
    //! ( 25  27  24  16  3 )
    // Note: imagPart is invalid for hermitian matrix, since the result 
    // would be anti-symmetric, which we don't have as a matrix type.
    std::cout<<"Norm(cm2) = "<<Norm(cm2)<<std::endl;
    //! Norm(cm2) = 202.349

    cm2 += cv^cv.conjugate();
    std::cout<<"cm2 += cv^cv.conjugate() =\n"<<cm2;
    //! cm2 += cv^cv.conjugate() =
    //! 5  5  
    //! ( (25,0)  (20,20)  (45,30)  (60,40)  (115,50) )
    //! ( (20,-20)  (17,0)  (32,34)  (42,44)  (72,54) )
    //! ( (45,-30)  (32,-34)  (59,0)  (79,38)  (159,48) )
    //! ( (60,-40)  (42,-44)  (79,-38)  (91,0)  (196,32) )
    //! ( (115,-50)  (72,-54)  (159,-48)  (196,-32)  (408,0) )
    cm2 += v^v;
    std::cout<<"cm2 += v^v =\n"<<cm2;
    //! cm2 += v^v =
    //! 5  5  
    //! ( (29,0)  (22,20)  (51,30)  (68,40)  (133,50) )
    //! ( (22,-20)  (18,0)  (35,34)  (46,44)  (81,54) )
    //! ( (51,-30)  (35,-34)  (68,0)  (91,38)  (186,48) )
    //! ( (68,-40)  (46,-44)  (91,-38)  (107,0)  (232,32) )
    //! ( (133,-50)  (81,-54)  (186,-48)  (232,-32)  (489,0) )
    cm2 -= cx*cx.adjoint();
    std::cout<<"cm2 -= cx*cx.adjoint() =\n"<<cm2;
    //! cm2 -= cx*cx.adjoint() =
    //! 5  5  
    //! ( (5,0)  (10,20)  (15,30)  (20,40)  (25,50) )
    //! ( (10,-20)  (12,0)  (17,34)  (22,44)  (27,54) )
    //! ( (15,-30)  (17,-34)  (14,0)  (19,38)  (24,48) )
    //! ( (20,-40)  (22,-44)  (19,-38)  (11,0)  (16,32) )
    //! ( (25,-50)  (27,-54)  (24,-48)  (16,-32)  (3,0) )


    // Can mix SymMatrix with other kinds of matrices:
    tmv::UpperTriMatrix<double> um(xm);
    std::cout<<"um + m1 =\n"<<um+m1;
    //! um + m1 =
    //! 5  5  
    //! ( 7  15  16  29  33 )
    //! ( 10  17  24  24  27 )
    //! ( 15  17  20  27  28 )
    //! ( 20  22  19  20  16 )
    //! ( 25  27  24  16  8 )
    tmv::DiagMatrix<double> dm(xm);
    std::cout<<"dm * m1 =\n"<<dm*m1;
    //! dm * m1 =
    //! 5  5  
    //! ( 10  20  30  40  50 )
    //! ( 50  60  85  110  135 )
    //! ( 90  102  84  114  144 )
    //! ( 180  198  171  99  144 )
    //! ( 125  135  120  80  15 )

    return 0;
} catch (tmv::Error& e) {
    std::cerr<<e<<std::endl;
    return 1;
}
