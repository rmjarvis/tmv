
#define TMV_EXTRA_DEBUG // To get the 888's below.

#include "TMV.h"
#include <iostream>

int main() try 
{
    // Several ways to create and initialize vectors:

    // Create with uninitialized values
    tmv::Vector<double> v1(6); 
    // In debug mode, all are initialized to 888 to make it easier to 
    // notice when you fail to initialize correctly.
    std::cout<<"v1 = "<<v1<<std::endl;
    //! v1 = 6  ( 888  888  888  888  888  888 )

    // Initialize with STL-compliant iterator
    srand(51572);
    std::generate(v1.begin(),v1.end(),rand);
    v1 /= double(RAND_MAX);
    std::cout<<"v1 = "<<v1<<std::endl;
    //! v1 = 6  ( 0.403622  0.66681  0.0776762  0.504604  0.871968  0.163157 )

    // Create with all 2's.
    tmv::Vector<double> v2(6,2.); 
    std::cout<<"v2 = "<<v2<<std::endl;
    //! v2 = 6  ( 2  2  2  2  2  2 )

    // Create with given elements from C-array
    double vv[6] = {1.1,8.,-15.,2.5,6.3,-12.};
    tmv::Vector<double> v3(6); 
    std::copy(vv,vv+6,v3.begin());
    std::cout<<"v3 = "<<v3<<std::endl;
    //! v3 = 6  ( 1.1  8  -15  2.5  6.3  -12 )

    // Initialize with direct access of each element
    tmv::Vector<double> v4(6); 
    for(int i=0;i<v4.size();i++) {
        v4(i) = 2*i+10.;  // Could also use v4[i] instead of v4(i)
    }
    std::cout<<"v4 = "<<v4<<std::endl;
    //! v4 = 6  ( 10  12  14  16  18  20 )

    // Initialize with comma-delimited list
    v4 << 1.2, 9., 12, 2.5, -7.4, 14;
    std::cout<<"v4 = "<<v4<<std::endl;
    //! v4 = 6  ( 1.2  9  12  2.5  -7.4  14 )
    // If the list is the wrong size, a run time error occurs:
    try {
        v4 << 1.2, 9., 12, 2.5;
    } catch (tmv::Error& e) {
        std::cout<<"Caught e = "<<e;
        //! Caught e = TMV Error: Reading from list initialization.
        //! Reached end of list, but expecting 2 more elements.
    }
    try {
        v4 << 1.2, 9., 12, 2.5, -7.4, 14, 99;
    } catch (tmv::Error& e) {
        std::cout<<"Caught e = "<<e;
        //! Caught e = TMV Error: Reading from list initialization.
        //! List has more elements than expected.
    }

    // Norms, etc. 
    std::cout<<"Norm1(v3) = "<<Norm1(v3);
    std::cout<<" = "<<v3.sumAbsElements()<<std::endl;
    //! Norm1(v3) = 44.9 = 44.9

    std::cout<<"Norm2(v3) = "<<Norm2(v3);
    std::cout<<" = "<<v3.norm()<<std::endl;
    //! Norm2(v3) = 21.9123 = 21.9123

    std::cout<<"NormInf(v3) = "<<NormInf(v3);
    std::cout<<" = "<<v3.maxAbsElement()<<std::endl;
    //! NormInf(v3) = 15 = 15

    std::cout<<"SumElements(v3) = "<<SumElements(v3)<<std::endl;
    //! SumElements(v3) = -9.1

    // Min/Max elements:
    int i1,i2,i3,i4;
    double x1 = v3.minAbsElement(&i1);
    double x2 = v3.maxAbsElement(&i2);
    double x3 = v3.minElement(&i3);
    double x4 = v3.maxElement(&i4);
    std::cout<<"|v3("<<i1<<")| = "<<x1<<" is the minimum absolute value\n";
    //! |v3(0)| = 1.1 is the minimum absolute value
    std::cout<<"|v3("<<i2<<")| = "<<x2<<" is the maximum absolute value\n";
    //! |v3(2)| = 15 is the maximum absolute value
    std::cout<<"v3("<<i3<<") = "<<x3<<" is the minimum value\n";
    //! v3(2) = -15 is the minimum value
    std::cout<<"v3("<<i4<<") = "<<x4<<" is the maximum value\n";
    //! v3(1) = 8 is the maximum value


    // Modifications:

    v1 << 1, 2, 4, 8, 16, 32;
    std::cout<<"v1 = "<<v1<<std::endl;
    //! v1 = 6  ( 1  2  4  8  16  32 )

    v1.addToAll(5.);
    std::cout<<"v1.addToAll(5.) = "<<v1<<std::endl;
    //! v1.addToAll(5.) = 6  ( 6  7  9  13  21  37 )

    v1.reverseSelf();
    std::cout<<"v1.reverseSelf() = "<<v1<<std::endl;
    //! v1.reverseSelf() = 6  ( 37  21  13  9  7  6 )

    v1.setZero();
    std::cout<<"v1.setZero() = "<<v1<<std::endl;
    //! v1.setZero() = 6  ( 0  0  0  0  0  0 )

    v1.setAllTo(20.);
    std::cout<<"v1.setAllTo(20.) = "<<v1<<std::endl;
    //! v1.setAllTo(20.) = 6  ( 20  20  20  20  20  20 )

    v1.makeBasis(2);
    std::cout<<"v1.makeBasis(2) = "<<v1<<std::endl;
    //! v1.makeBasis(2) = 6  ( 0  0  1  0  0  0 )


    // Views:

    std::cout<<"v3 = "<<v3<<std::endl;
    //! v3 = 6  ( 1.1  8  -15  2.5  6.3  -12 )
    std::cout<<"v3.subVector(0,3) = "<<v3.subVector(0,3)<<std::endl;
    //! v3.subVector(0,3) = 3  ( 1.1  8  -15 )
    std::cout<<"v3.subVector(0,6,2) = "<<v3.subVector(0,6,2)<<std::endl;
    //! v3.subVector(0,6,2) = 3  ( 1.1  -15  6.3 )
    std::cout<<"v3.reverse() = "<<v3.reverse()<<std::endl;
    //! v3.reverse() = 6  ( -12  6.3  2.5  -15  8  1.1 )

    // Views can be initialized with << too.
    v3.reverse() << 1.2, 9., 12, 2.5, -7.4, 14;
    std::cout<<"v3 = "<<v3<<std::endl;
    //! v3 = 6  ( 14  -7.4  2.5  12  9  1.2 )
    v3.subVector(0,6,2) << 18, 22, 33;
    std::cout<<"v3 = "<<v3<<std::endl;
    //! v3 = 6  ( 18  -7.4  22  12  33  1.2 )

    // Can use the views within expressions
    v3.reverse() += v4;
    std::cout<<"v3.reverse() += v4 => v3 = "<<v3<<std::endl;
    //! v3.reverse() += v4 => v3 = 6  ( 32  -14.8  24.5  24  42  2.4 )

    v3.subVector(0,3) *= 2.;
    std::cout<<"v3.subVector(0,3) *= 2 => v3 = "<<v3<<std::endl;
    //! v3.subVector(0,3) *= 2 => v3 = 6  ( 64  -29.6  49  24  42  2.4 )


    // Fortran Indexing:

    tmv::Vector<double,tmv::FortranStyle> fv3 = v3;

    std::cout<<"fv3 = v3 = "<<fv3<<std::endl;
    //! fv3 = v3 = 6  ( 64  -29.6  49  24  42  2.4 )
    std::cout<<"fv3(1) = "<<fv3(1)<<std::endl;
    //! fv3(1) = 64
    std::cout<<"fv3(6) = "<<fv3(6)<<std::endl;
    //! fv3(6) = 2.4
    std::cout<<"fv3.subVector(1,3) = "<<fv3.subVector(1,3)<<std::endl;
    //! fv3.subVector(1,3) = 3  ( 64  -29.6  49 )
    std::cout<<"fv3.makeBasis(3) = "<<fv3.makeBasis(3)<<std::endl;
    //! fv3.makeBasis(3) = 6  ( 0  0  1  0  0  0 )


    // Vector arithmetic:

    tmv::Vector<double> v3pv4 = v3 + v4;  
    std::cout<<"v3 + v4 = "<<v3pv4<<std::endl;
    //! v3 + v4 = 6  ( 65.2  -20.6  61  26.5  34.6  16.4 )

    // Inner product
    double v3v4 = v3 * v4; 
    std::cout<<"v3 * v4 = "<<v3v4<<std::endl;
    //! v3 * v4 = 181.2

    v3 *= 2.;  
    std::cout<<"v3 *= 2 = "<<v3<<std::endl;
    //! v3 *= 2 = 6  ( 128  -59.2  98  48  84  4.8 )

    v3 += v4; 
    std::cout<<"v3 += v4 = "<<v3<<std::endl;
    //! v3 += v4 = 6  ( 129.2  -50.2  110  50.5  76.6  18.8 )

    // Get as complicated as you want:
    std::cout<<"(v1*v2) * v3 + ((-v4 + 4.*v1)*v2) * v2/20. =\n"<<
        (v1*v2) * v3 + ((-v4 + 4.*v1)*v2) * v2/20.<<std::endl;
    //! (v1*v2) * v3 + ((-v4 + 4.*v1)*v2) * v2/20. =
    //! 6  ( 252.94  -105.86  214.54  95.54  147.74  32.14 )

    // Automatically checks for aliases:
    v3 = v4 - 3.*v3;
    std::cout<<"v3 = v4-3.*v3 => v3 = "<<v3<<std::endl;
    //! v3 = v4-3.*v3 => v3 = 6  ( -386.4  159.6  -318  -149  -237.2  -42.4 )


    // Complex vectors:

    tmv::Vector<std::complex<double> > cv4 = v4 * std::complex<double>(1,2);
    std::cout<<"cv4 = v4 * (1+2i) =\n"<<cv4<<std::endl;
    //! cv4 = v4 * (1+2i) =
    //! 6  ( (1.2,2.4)  (9,18)  (12,24)  (2.5,5)  (-7.4,-14.8)  (14,28) )

    std::cout<<"cv4.conjugate() =\n"<<cv4.conjugate()<<std::endl;
    //! cv4.conjugate() =
    //! 6  ( (1.2,-2.4)  (9,-18)  (12,-24)  (2.5,-5)  (-7.4,14.8)  (14,-28) )
    std::cout<<"cv4.realPart() = "<<cv4.realPart()<<std::endl;
    //! cv4.realPart() = 6  ( 1.2  9  12  2.5  -7.4  14 )
    std::cout<<"cv4.imagPart() = "<<cv4.imagPart()<<std::endl;
    //! cv4.imagPart() = 6  ( 2.4  18  24  5  -14.8  28 )
    std::cout<<"Norm(cv4) = "<<Norm(cv4)<<std::endl;
    //! Norm(cv4) = 49.1655
    std::cout<<"sqrt(cv4*cv4.conjugate()) = "<<
        sqrt(cv4*cv4.conjugate())<<std::endl;
    //! sqrt(cv4*cv4.conjugate()) = (49.1655,0)
    std::cout<<"cv4.maxAbsElement() = "<<cv4.maxAbsElement()<<std::endl;
    //! cv4.maxAbsElement() = 31.305

    // Can mix real and complex in any combination
    std::cout<<"cv4 - v4 = "<<cv4 - v4<<std::endl;
    //! cv4 - v4 = 6  ( (0,2.4)  (0,18)  (0,24)  (0,5)  (0,-14.8)  (0,28) )
    std::cout<<"cv4 * v4 * (1-2i) = "<<
        cv4 * v4 * std::complex<double>(1,-2)<<std::endl;
    //! cv4 * v4 * (1-2i) = (2417.25,0)


    // Sorting:

    v4 << 2, 5.3, -1.5, -7, 0.5, -2.8;
    std::cout<<"v4 = "<<v4<<std::endl;
    //! v4 = 6  ( 2  5.3  -1.5  -7  0.5  -2.8 )
    tmv::Permutation P;
    v4.sort(P);
    std::cout<<"Sorted: v4 = "<<v4<<std::endl;
    //! Sorted: v4 = 6  ( -7  -2.8  -1.5  0.5  2  5.3 )
    v4 = P.inverse() * v4;
    std::cout<<"Sort undone: v4 = "<<v4<<std::endl;
    //! Sort undone: v4 = 6  ( 2  5.3  -1.5  -7  0.5  -2.8 )
    v4.sort(); // Don't necessarily need P.
    std::cout<<"Resorted: v4 = "<<v4<<std::endl;
    //! Resorted: v4 = 6  ( -7  -2.8  -1.5  0.5  2  5.3 )

    // Can sort by other criteria:
    std::cout<<"v4.sort(Descend) = "<<v4.sort(tmv::Descend)<<std::endl;
    //! v4.sort(Descend) = 6  ( 5.3  2  0.5  -1.5  -2.8  -7 )

    cv4.realPart() << -3,  1, -2, -1,  7,  3;
    cv4.imagPart() <<  4, -1,  0, -6,  5, -1;
    // (I find this complex initialization to be more readable than a list 
    //  filled with values that look like complex<double>(-3,4), ...)
    std::cout<<"cv4 = "<<cv4<<std::endl;
    //! cv4 = 6  ( (-3,4)  (1,-1)  (-2,0)  (-1,-6)  (7,5)  (3,-1) )

    std::cout<<"cv4.sort(Descend,RealComp) =\n"<<
        cv4.sort(tmv::Descend,tmv::RealComp)<<std::endl;
    //! cv4.sort(Descend,RealComp) =
    //! 6  ( (7,5)  (3,-1)  (1,-1)  (-1,-6)  (-2,0)  (-3,4) )
    std::cout<<"cv4.sort(Ascend,ImagComp) =\n"<<
        cv4.sort(tmv::Ascend,tmv::ImagComp)<<std::endl;
    //! cv4.sort(Ascend,ImagComp) =
    //! 6  ( (-1,-6)  (3,-1)  (1,-1)  (-2,0)  (-3,4)  (7,5) )
    std::cout<<"cv4.sort(Ascend,AbsComp) =\n"<<
        cv4.sort(tmv::Ascend,tmv::AbsComp)<<std::endl;
    //! cv4.sort(Ascend,AbsComp) =
    //! 6  ( (1,-1)  (-2,0)  (3,-1)  (-3,4)  (-1,-6)  (7,5) )
    std::cout<<"cv4.sort(Ascend,ArgComp) =\n"<<
        cv4.sort(tmv::Ascend,tmv::ArgComp)<<std::endl;
    //! cv4.sort(Ascend,ArgComp) =
    //! 6  ( (-1,-6)  (1,-1)  (3,-1)  (7,5)  (-3,4)  (-2,0) )

    // The default component for complex vectors is RealComp:
    std::cout<<"cv4.sort() =\n"<< cv4.sort()<<std::endl;
    //! cv4.sort() =
    //! 6  ( (-3,4)  (-2,0)  (-1,-6)  (1,-1)  (3,-1)  (7,5) )

    return 0;
} catch (tmv::Error& e) {
    std::cerr<<e<<std::endl;
    return 1;
}
