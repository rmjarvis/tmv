
#define TMV_EXTRA_DEBUG 

#include "TMV.h"
#include "TMV_Diag.h"
#include <iostream>

int main() try 
{
    tmv::Matrix<double> A(4,4);
    for(int i=0;i<A.nrows();i++) 
        for(int j=0;j<A.ncols();j++) 
            A(i,j) = 6.*i-2*j*j+2.; 
    A.diag().addToAll(5.);

    tmv::Vector<double> b(4);
    for(int i=0;i<b.size();i++) 
        b(i) = 3.+2.*i; 

    // Basic Ax=b solution:
    std::cout<<"A =\n"<<A;
    //! A =
    //! 4  4  
    //! ( 7  0  -6  -16 )
    //! ( 8  11  0  -10 )
    //! ( 14  12  11  -4 )
    //! ( 20  18  12  7 )
    std::cout<<"b = "<<b<<std::endl;
    //! b = 4  ( 3  5  7  9 )

    // Solve: Ax = b
    tmv::Vector<double> x = b/A; // Default: use LU decomposition
    std::cout<<"x = b/A = "<<x<<std::endl;
    //! x = b/A = 4  ( 0.294545  0.170909  0.0472727  -0.0763636 )
    std::cout<<"Check: A*x = "<<A*x<<std::endl;
    //! Check: A*x = 4  ( 3  5  7  9 )

    // Can update A and then re-solve:
    A(0,0) = 20.;
    x = b/A;
    std::cout<<"Now x = b/A = "<<x<<std::endl;
    //! Now x = b/A = 4  ( 0.133458  0.29877  0.117091  -0.064587 )
    std::cout<<"A*x = "<<A*x<<std::endl;
    //! A*x = 4  ( 3  5  7  9 )

    // If the matrix won't change, but you want to solve with multiple 
    // vectors, then it's faster to let TMV know this with saveDiv().
    tmv::Matrix<double> A2 = A;
    A2.saveDiv();
    x = b/A2;
    std::cout<<"x1 = "<<x<<std::endl;
    //! x1 = 4  ( 0.133458  0.29877  0.117091  -0.064587 )
    std::cout<<"A*x = "<<A2*x<<std::endl;
    //! A*x = 4  ( 3  5  7  9 )
    x = b.reverse()/A2; // Fast, since doesn't recalculate LU decomposition.
    std::cout<<"x2 = "<<x<<std::endl;
    //! x2 = 4  ( 0.136753  0.207381  -0.0775483  -0.362478 )
    std::cout<<"A*x2 = "<<A2*x<<std::endl;
    //! A*x2 = 4  ( 9  7  5  3 )
    A2.row(0) *= 2.;
    x = b/A2;  // Wrong, since doesn't recalculate LU decomposition.
    std::cout<<"Wrong x = "<<x<<std::endl;
    //! Wrong x = 4  ( 0.133458  0.29877  0.117091  -0.064587 )
    std::cout<<"A*x = "<<A2*x<<std::endl;
    //! A*x = 4  ( 6  5  7  9 )
    // If the matrix does change when saveDiv() is set, 
    // you can manually recalculate the decomposition:
    A2.resetDiv();
    x = b/A2;  // Now it is correct.
    std::cout<<"x = "<<x<<std::endl;
    //! x = 4  ( 0.0703537  0.348858  0.144442  -0.0599736 )
    std::cout<<"A*x = "<<A2*x<<std::endl;
    //! A*x = 4  ( 3  5  7  9 )

    // Matrix inverse:
    tmv::Matrix<double> A2inv = A2.inverse();
    std::cout<<"Ainv =\n"<<A2inv;
    //! Ainv =
    //! 4  4  
    //! ( 0.0210347  -0.0395431  -0.012522  0.0325132 )
    //! ( -0.016696  0.0966608  -0.0527241  0.0316344 )
    //! ( -0.00911687  -0.0643234  0.139631  -0.0537786 )
    //! ( -0.00153779  -0.0253076  -0.0680141  0.0608084 )
    std::cout<<"Ainv*A =\n"<<A2inv*A2;
    //! Ainv*A =
    //! 4  4  
    //! ( 1  0  5.55112e-17  1.66533e-16 )
    //! ( 1.11022e-16  1  -1.66533e-16  -4.996e-16 )
    //! ( 0  2.22045e-16  1  2.22045e-16 )
    //! ( 0  -2.22045e-16  -1.11022e-16  1 )
    
    // This is a case where it can be useful to see the 
    // matrix elements that are larger than some threshold value:
    std::cout<<"Ainv*A =\n" << tmv::ThreshIO(1.e-8) << (A2inv*A2);
    //! Ainv*A =
    //! 4  4  
    //! ( 1  0  0  0 )
    //! ( 0  1  0  0 )
    //! ( 0  0  1  0 )
    //! ( 0  0  0  1 )

    // 1/x notation is treated as arithmetic:
    // (But the 1 has to be the same type as the elements of the matrix.
    /// In this case, double.  With a float matrix, use 1.F instead.)
    std::cout<<"1./A =\n"<<1./A2;
    //! 1./A =
    //! 4  4  
    //! ( 0.0210347  -0.0395431  -0.012522  0.0325132 )
    //! ( -0.016696  0.0966608  -0.0527241  0.0316344 )
    //! ( -0.00911687  -0.0643234  0.139631  -0.0537786 )
    //! ( -0.00153779  -0.0253076  -0.0680141  0.0608084 )
    std::cout<<"5./A =\n"<<5./A2;
    //! 5./A =
    //! 4  4  
    //! ( 0.105174  -0.197715  -0.0626098  0.162566 )
    //! ( -0.0834798  0.483304  -0.26362  0.158172 )
    //! ( -0.0455844  -0.321617  0.698155  -0.268893 )
    //! ( -0.00768893  -0.126538  -0.34007  0.304042 )

    // Can also use inverse() notation instead of /
    // x = b/A2 inlines to exactly the same thing as x = A2.inverse() * b;
    // ie. accurate back-substitution methods are used rather than 
    // actually computing the inverse and then multiplying:
    x =  A2.inverse() * b;
    std::cout<<"x = A.inverse() * b = "<<x<<std::endl;
    //! x = A.inverse() * b = 4  ( 0.0703537  0.348858  0.144442  -0.0599736 )

    // Division from the other side can either be done with this
    // inverse() notation or with the % operator:
    // This is the solution to the equation x A = b, rather than A x = b.
    x = b * A2.inverse();
    std::cout<<"x = b * A.inverse() = "<<x<<std::endl;
    //! x = b * A.inverse() = 4  ( -0.0980338  -0.313357  0.0641037  0.426538 )
    x = b % A2;
    std::cout<<"x = b % A = "<<x<<std::endl;
    //! x = b % A = 4  ( -0.0980338  -0.313357  0.0641037  0.426538 )
    std::cout<<"Check: x*A = "<<x*A2<<std::endl;
    //! Check: x*A = 4  ( 3  5  7  9 )

    tmv::Matrix<double> B(4,3);
    for(int i=0;i<B.nrows();i++) 
        for(int j=0;j<B.ncols();j++) 
            B(i,j) = 1.+2.*i+j*j;
    // Multiple right hand sides may be calculated at once if 
    // B is a matrix, rather than a vector:
    std::cout<<"B =\n"<<B;
    //! B =
    //! 4  3  
    //! ( 1  2  5 )
    //! ( 3  4  7 )
    //! ( 5  6  9 )
    //! ( 7  8  11 )
    tmv::Matrix<double> X = B/A2;
    std::cout<<"X = B/A =\n"<<X;
    //! X = B/A =
    //! 4  3  
    //! ( 0.067388  0.0688708  0.0733194 )
    //! ( 0.231107  0.289982  0.466608 )
    //! ( 0.119618  0.13203  0.169266 )
    //! ( 0.0081283  -0.0259227  -0.128076 )
    std::cout<<"AX = "<<A2*X;
    //! AX = 4  3  
    //! ( 1  2  5 )
    //! ( 3  4  7 )
    //! ( 5  6  9 )
    //! ( 7  8  11 )

    // And as always, you can mix complex and real objects:
    tmv::Vector<std::complex<double> > cb = b * std::complex<double>(3,-2);
    cb(1) = std::complex<double>(-1,8);
    cb(2) = std::complex<double>(1,1);
    std::cout<<"cb = "<<cb<<std::endl;
    //! cb = 4  ( (9,-6)  (-1,8)  (1,1)  (27,-18) )
    tmv::Vector<std::complex<double> > cx = cb/A;
    std::cout<<"cx = cb/A =\n"<<cx<<std::endl;
    //! cx = cb/A =
    //! 4  ( (1.2835,-1.16652)  (0.404218,0.351494)  (-1.41217,0.70246)  (1.57144,-1.34657) )
    std::cout<<"A*cx = "<<A*cx<<std::endl;
    //! A*cx = 4  ( (9,-6)  (-1,8)  (1,1)  (27,-18) )
    tmv::Matrix<std::complex<double> > CA = A * std::complex<double>(5,-2);
    CA(1,1) = std::complex<double>(1,6);
    CA(2,3) = std::complex<double>(4,-1);
    std::cout<<"CA = "<<CA;
    //! CA = 4  4  
    //! ( (100,-40)  (0,0)  (-30,12)  (-80,32) )
    //! ( (40,-16)  (1,6)  (0,0)  (-50,20) )
    //! ( (70,-28)  (60,-24)  (55,-22)  (4,-1) )
    //! ( (100,-40)  (90,-36)  (60,-24)  (35,-14) )
    std::cout<<"cx = b/CA =\n"<<(cx=b/CA)<<std::endl;
    //! cx = b/CA =
    //! 4  ( (-0.0858093,0.219988)  (0.270489,-0.423628)  (-0.0665359,0.214597)  (-0.114638,0.18158) )
    std::cout<<"CA*cx = "<<tmv::ThreshIO(1.e-8)<<CA*cx<<std::endl;
    //! CA*cx = 4  ( (3,0)  (5,0)  (7,0)  (9,0) )
    std::cout<<"cb/CA =\n"<<(cx=cb/CA)<<std::endl;
    //! cb/CA =
    //! 4  ( (0.0374442,0.698211)  (0.792852,-1.67098)  (-0.916414,0.915387)  (0.267616,0.555355) )
    std::cout<<"CA*cx = "<<CA*cx<<std::endl;
    //! CA*cx = 4  ( (9,-6)  (-1,8)  (1,1)  (27,-18) )


    // Least-squares solutions:
    // If A in the matrix equation A x = b has more rows than columns,
    // then there is, in general, no solution to the equation.
    // Instead, one is generally looking for the x that comes closest
    // to satisfying the equation, in a least-squares sense.
    // Specifically, the x for which Norm2(b-Ax) is minimized.
    // This is the solution produced by TMV for such matrices.
    //
    // Here I model a theoretical system for which each observation
    // is 5 + 6i - 3i^2 in the absense of measurement errors.
    // We observe the actual values which are not quite equal to 
    // that because of noise.  And the goal is to determine the 
    // coefficients of 1,i,i^2 (5,6,-3) from noisy observations:
    tmv::Vector<double> b3(10);
    tmv::Matrix<double> A3(10,3);
    double errors[10] = {0.01,-0.02,0.02,-0.02,0.00,-0.03,0.01,-0.02,0.03,0.02};
    for(int i=0;i<10;i++) {
        b3(i) = 5. + 6.*i - 3.*i*i + errors[i]; // Model of measurements
        A3(i,0) = 1.;  // Parameterization of the model...
        A3(i,1) = i;
        A3(i,2) = i*i;
    }
    double sigma = 0.02; // sigma = estimate of rms errors
    A3 /= sigma;  
    b3 /= sigma;

    tmv::Vector<double> x3 = b3/A3;  // Uses QR decomposition by default
    std::cout<<"x = "<<x3<<std::endl;
    //! x = 3  ( 5.00773  5.98989  -2.99867 )
    std::cout<<"A*x => \n"<<A3*x3<<std::endl;
    //! A*x => 
    //! 10  ( 250.386  399.947  249.64  -200.534  -950.576  -2000.48  -3350.26  -4999.91  -6949.42  -9198.8 )
    std::cout<<"chisq = NormSq(A*x-b) = "<<NormSq(A3*x3-b3)<<std::endl;
    //! chisq = NormSq(A*x-b) = 6.99811
    // The expected value for this is 10 observations minus 3 degrees
    // of freedom = 7.

    // The covariance matrix for x is (A.Transpose() * A)^-1
    // This combination is easy to calculate from the QR decomposition
    // that TMV has used to do the division.  Therefore, we provide
    // it as an explicit function:
    tmv::Matrix<double> cov(3,3);
    A3.makeInverseATA(cov);
    std::cout<<"Cov(x) =\n"<<cov;
    //! Cov(x) =
    //! 3  3  
    //! ( 0.000247273  -0.000103636  9.09091e-06 )
    //! ( -0.000103636  6.62121e-05  -6.81818e-06 )
    //! ( 9.09091e-06  -6.81818e-06  7.57576e-07 )

    // The singular value decomposition can detect ill-conditioned matrices
    // and correct for them.
    // For example, if you model the above observations with the 
    // 1,i,i^2 components as before, but add as well 6*i-5 as a component,
    // then that would be degenerate with 1 and i.  
    // SVD is able to detect this defect and deal with it appropriately:
    tmv::Matrix<double> A4(10,4);
    A4.colRange(0,3) = A3*sigma;
    for(int i=0;i<10;i++) A4(i,3) = 6.*i-5.;
    std::cout<<"Now A*sigma =\n"<<A4;
    //! Now A*sigma =
    //! 10  4  
    //! ( 1  0  0  -5 )
    //! ( 1  1  1  1 )
    //! ( 1  2  4  7 )
    //! ( 1  3  9  13 )
    //! ( 1  4  16  19 )
    //! ( 1  5  25  25 )
    //! ( 1  6  36  31 )
    //! ( 1  7  49  37 )
    //! ( 1  8  64  43 )
    //! ( 1  9  81  49 )
    A4 /= sigma;
    try {
        // This may or may not succeed, but if it does, the results will be 
        // unusable, typically with values around 1.e13 and such.
        tmv::Vector<double> x4 = b3/A4;
        std::cout<<"Unstable x = b/A =\n"<<x4<<std::endl;
        //! Unstable x = b/A =
        //! 4  ( -1.31399e+13  1.57679e+13  -2.99808  -2.62798e+12 )
        std::cout<<"A*x =\n"<<A4*x4<<std::endl;
        //! A*x =
        //! 10  ( 251.125  400.203  249.5  -201.25  -951.5  -2001  -3351.5  -5000  -6950  -9199 )
        std::cout<<"chisq = "<<NormSq(A4*x4-b3);
        std::cout<<"      Norm(x) = "<<Norm(x4)<<std::endl;
        //! chisq = 13.9006      Norm(x) = 2.06927e+13
    } catch (tmv::Error& e) {
        std::cout<<"Tried x = b/A, but caught error: \n"<<e<<std::endl;
    }

    // So instead, tell TMV to use SVD for division rather than QR.
    A4.divideUsing(tmv::SV);
    std::cout<<"Singular values for A are "<<A4.svd().getS().diag()<<std::endl;
    //! Singular values for A are 4  ( 7618.67  733.074  116.312  5.15828e-14 )
    std::cout<<"Using only the first "<<A4.svd().getKMax()<<" components\n";
    //! Using only the first 3 components
    tmv::Vector<double> x4 = b3/A4;
    std::cout<<"SVD division yields: x = "<<x4<<std::endl;
    //! SVD division yields: x = 4  ( 5.88681  4.93498  -2.99867  0.175817 )
    std::cout<<"chisq = "<<NormSq(A4*x4-b3);
    std::cout<<"      Norm(x) = "<<Norm(x4)<<std::endl;
    //! chisq = 6.99811      Norm(x) = 8.24813

    // QRP can also give useful results, but isn't quite as flexible 
    // as SVD:
    A4.divideUsing(tmv::QRP);
    x4 = b3/A4;
    std::cout<<"QRP division yields: x = "<<x4<<std::endl;
    std::cout<<"chisq = "<<NormSq(A4*x4-b3);
    std::cout<<"      Norm(x) = "<<Norm(x4)<<std::endl;
    //! QRP division yields: x = 4  ( 5.00773  5.98989  -2.99867  0 )
    //! chisq = 6.99811      Norm(x) = 8.3635

    // Note that both methods give answers with an equal chisq, so
    // Ax is equally close to b for each but they have different 
    // specific choices for the degeneracy.
    // SVD will give the solution within this degeneracy freedom that
    // has the minimum Norm(x), but QRP will be faster -- 
    // significantly so for large matrices.

    return 0;
} catch (tmv::Error& e) {
    std::cerr<<e<<std::endl;
    return 1;
}
