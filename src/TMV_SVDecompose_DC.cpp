///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//#define XDEBUG

#include "TMV_SVDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Givens.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::endl;

#ifdef XDEBUG
#ifdef _OPENMP

#include <sstream>
struct ThreadSafeWriter
{
    ThreadSafeWriter() {}
    ~ThreadSafeWriter()
    {
#pragma omp critical
        {
            std::cout<<s.str();
        }
    }
    template <class T>
    ThreadSafeWriter& operator<<(const T& x)
    { s << x; return *this; }

    // The next bit is to get dbgcout<<std::endl; working.
    // See: http://stackoverflow.com/questions/1134388/stdendl-is-of-unknown-type-when-overloading-operator

    typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
    typedef CoutType& (*StandardEndLine)(CoutType&);
    ThreadSafeWriter& operator<<(StandardEndLine)
    {
#pragma omp critical
        {
            std::cout<<s.str()<<std::endl;
        }
        s.clear();
        s.str(std::string());
        return *this;
    }

    std::stringstream s;
};

#include <omp.h>

#define dbgcout ThreadSafeWriter()

#else

#define dbgcout std::cout
//#define dbgcout if (false) std::cout

#endif

#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixArith.h"
#define THRESH 1.e-5
//#define TESTUV  // Should only use this for full USVt decompositions
using std::cerr;
#else
#define dbgcout if (false) std::cout
#endif

namespace tmv {

#define RT TMV_RealType(T)

#ifdef TMV_BLOCKSIZE
#define DC_LIMIT TMV_BLOCKSIZE/2
#else
#define DC_LIMIT 32
#endif

    // Note about OpenMP here.
    // Both the divide step (recursing on SV_DecomposeFromBidiagonal_DC) 
    // and the conquer step (the loop in FindDCSingularValues) are
    // pontentially parallelizable.  However, recursive functions do not
    // lend themselves to being parallelized by OPENMP very easily.
    // So I think I need to rewrite the divide step in a way that does
    // not rely on recursion before I apply OPENMP pragmas to it.
    // So for now, we only have recursion in the FindDCSingularValues function.

    template <class T> 
    static T FindDCSingularValue(
        const ptrdiff_t k, const ptrdiff_t N, const T rho, const T* D, const T* z, 
        const T* zsq, const T normsqz, T* diff, T* sum)
    {
        dbgcout<<"start FindDCV: k = "<<k<<endl;
        // Find solutions to the equation:
        // f(s) = rho + Sum_i=1..N z_i^2/(D_i^2-s^2)
        //
        // diff outputs the values of D(i)-s more accurately than
        // a direct calculation can usually do.
        //
        // For this, we follow the advice of Ren-Cang Li's 1993 paper
        // Solving Secular Equations Stably and Efficiently.
        // I won't give the derivations of any of this.  For more details
        // you can see his paper, which I found on the net at:
        // http://www.netlib.org/lapack/lawnspdf/lawn89.pdf
        //
        const T eps = TMV_Epsilon<T>();
        const ptrdiff_t TMV_MAXITER = 100;

        dbgcout<<"k = "<<k<<endl;
        if (k<N-1) {
            dbgcout<<"root bracketed by "<<D[k]<<"  "<<D[k+1]<<endl;
            dbgcout<<"   range = "<<D[k+1]-D[k]<<endl;
            dbgcout<<"s^2 bracketed by "<<D[k]*D[k]<<"  "<<D[k+1]*D[k+1]<<endl;
            dbgcout<<"   range = "<<D[k+1]*D[k+1]-D[k]*D[k]<<endl;
        } else {
            dbgcout<<"root bracketed by "<<D[k]<<"  "<<
                TMV_SQRT(D[k]*D[k]+normsqz)<<endl;
            dbgcout<<"   range = "<<TMV_SQRT(D[k]*D[k]+normsqz)-D[k]*D[k]<<endl;
            dbgcout<<"s^2 bracketed by "<<D[k]*D[k]<<"  "<<
                D[k]*D[k]+normsqz<<endl;
            dbgcout<<"   range = "<<normsqz<<endl;
        }

        // Most of the routine is the same for k=N-1 as for the normal
        // case, but there are some lines where the normal case uses
        // k, but when k==N-1, I need k-1.  So I use kk for these.
        ptrdiff_t kk = k < N-1 ? k : k-1;

        // First we need an initial guess for s_k 
        dbgcout<<"D[kk+1]-D[kk] = "<<D[kk+1]-D[kk]<<endl;
        dbgcout<<"eps*D["<<kk+1<<"] = "<<eps*D[kk+1]<<endl;
        T delta = (D[kk+1] - D[kk]) * (D[kk+1] + D[kk]);
        // It is known that D_k < s_k < D_k+1
        // Let s^2 = (D_k^2 + D_k+1^2)/2
        // Define t,tau: s = D_k + t, and s^2 = D_k^2 + tau
        T tau = k < N-1 ? delta/T(2) : normsqz/(T(2)*rho);
        T t = tau / (D[k] + TMV_SQRT(D[k]*D[k]+tau));
        T s = D[k] + t;
        dbgcout<<"Initial s = "<<s<<" = "<<D[k]<<" + "<<t<<std::endl;

        for(ptrdiff_t j=0;j<N;j++) sum[j] = D[j]+s; 
        for(ptrdiff_t j=0;j<N;j++) diff[j] = D[j]-D[k]; 
        for(ptrdiff_t j=0;j<N;j++) diff[j] -= t; 

        // Let c = f(s) - z_k^2/(D_k^2-s^2) - z_k+1^2/(D_k+1^2-s^2)
        // i.e. c is f(s) without the two terms for the poles that
        // bracket the root we are seraching for.
        // Since the largest terms are normally the ones for the poles nearer
        // to s, we perform the calculation from each end and proceed towards k.
        // Note: for k==N-1, the two poles do not bracket the root, but 
        // we proceed similarly anyway.
        T psi(0);
        for(ptrdiff_t j=0;j<kk;j++) {
            psi += (zsq[j] / sum[j]) / diff[j];
            //dbgcout<<"psi += "<<zsq[j]<<" / ("<<diff[j]<<" * "<<sum[k]<<") = "<<psi<<endl;
        }
        T phi(0);
        for(ptrdiff_t j=N-1;j>kk+1;j--) {
            phi += (zsq[j] / sum[j]) / diff[j];
            //dbgcout<<"phi += "<<zsq[j]<<" / ("<<diff[j]<<" * "<<sum[k]<<") = "<<phi<<endl;
        }
        T c = rho + psi + phi;
        dbgcout<<"c = "<<c<<endl;

        // Finish calculating f = f(s)
        T f = c + (zsq[kk] / sum[kk]) / diff[kk] + 
            (zsq[kk+1] / sum[kk+1]) / diff[kk+1];
        dbgcout<<"f = c + "<<zsq[kk]<<" / ("<<diff[kk]<<" * "<<sum[kk]<<")\n";
        dbgcout<<"      + "<<zsq[kk+1]<<" / ("<<
            diff[kk+1]<<" * "<<sum[kk+1]<<")\n";
        dbgcout<<"f("<<tau<<") = "<<f<<endl;

        T lowerbound, upperbound; // the allowed range for solution tau
        ptrdiff_t k1; // The index of the pole closest to the root
        if (k < N-1) {
            dbgcout<<"k < N-1\n";
            if (f >= T(0)) { // D_k < s_k < s
                dbgcout<<"f >= 0\n";
                k1 = k;
                lowerbound = T(0);
                upperbound = tau;
                dbgcout<<"bounds = 0,tau = "<<tau<<endl;
                T a = c + zsq[k]/delta + zsq[k+1]/delta;
                T b = zsq[k]/delta;
                T d = a*a-T(4)*b*c;
                if (d < T(0)) d = T(0); // Then rounding error.  Use d = 0
                else d = TMV_SQRT(d);
                // Our initial estimate for s_k^2 is D_k^2 + tau
                tau = (a > 0) ? T(2)*b/(a+d) : (a-d)/(T(2)*c);
                tau *= delta;
                if (tau > upperbound) tau = upperbound;
            } else {
                dbgcout<<"f < 0\n";
                k1 = k+1;
                lowerbound = -tau;
                upperbound = T(0);
                dbgcout<<"bounds = -tau,0 = "<<-tau<<endl;
                T a = -c + zsq[k]/delta + zsq[k+1]/delta;
                T b = -zsq[k+1]/delta;
                T d = a*a-T(4)*b*c;
                if (d < T(0)) d = T(0); 
                else d = TMV_SQRT(d);
                // Our initial estimate for s_k^2 is D_k+1^2 + tau 
                // (where tau < 0)
                tau = (a > 0) ? T(2)*b/(a+d) : (a-d)/(T(2)*c);
                tau *= delta;
                if (tau < lowerbound) tau = lowerbound;
            }
        } else { // k == N-1
            dbgcout<<"k == N-1\n";
            k1 = k;
            do {
                if (f >= T(0)) { // D_k < s_k < s
                    dbgcout<<"f >= 0\n";
                    lowerbound = 0;
                    upperbound = tau;
                    dbgcout<<"bounds = 0,tau = "<<tau<<endl;
                } else { // s < s_k < sqrt(D_k^2 + |z|^2)
                    dbgcout<<"f < 0\n";
                    lowerbound = tau;
                    upperbound = normsqz/rho;
                    dbgcout<<"bounds = tau,2tau = "<<tau<<"  "<< normsqz/rho<<endl;
                    if (c <= T(0)) {
                        // In this case there is no positive solution.  Best guess
                        // is just to start at upper bound.
                        tau = upperbound;
                        break;
                    }
                }
                T a = -c + zsq[k-1]/delta + zsq[k]/delta;
                T b = -zsq[k]/delta;
                T d = a*a-T(4)*b*c;
                if (d < T(0)) d = T(0); 
                else d = TMV_SQRT(d);
                tau = (a < 0) ? T(2)*b/(a-d) : (a+d)/(T(2)*c);
                tau *= delta;
                if (tau > upperbound) tau = upperbound;
                // The other times we check this are basically for rounding 
                // error reasons.
                // But this time, if f < 0, the solution to the quadratic
                // can really be much > upperbound.  So the check is especially
                // important in this case.
            } while (false);
        }
        dbgcout<<"Initial tau = "<<tau<<endl;
        dbgcout<<"Allowed bounds for tau = "<<
            lowerbound<<"  "<<upperbound<<endl;
        TMVAssert(tau >= lowerbound && tau <= upperbound);
        t = tau / ( TMV_SQRT(D[k1]*D[k1] + tau) + D[k1] );

        s = D[k1] + t;
        dbgcout<<"First refinement: s = "<<s<<" = "<<D[k1]<<" + "<<t<<std::endl;
        // We should be able to write diff[j] = (D[j]-D[k1)-t,
        // but some compilers do optimizations that lose accuracy 
        // and can end up with diff = 0.
        // So need to do this in two steps.
        for(ptrdiff_t j=0;j<N;j++) sum[j] = D[j]+s;
        for(ptrdiff_t j=0;j<N;j++) diff[j] = D[j]-D[k1];
        for(ptrdiff_t j=0;j<N;j++) diff[j] -= t; 

        // Define f(s) = rho + psi(s) + z_k1^2/(D_k1^2-s^2) + phi(s)
        // psi(s) = Sum_k=1..k1-1 z_k^2/(D_k^2-s^2)
        // phi(s) = Sum_k=k1+1..N z_k^2/(D_k^2-s^2)

        T e, dpsi, dphi, fk, dfk, k1term, dk1term, df;
        // e is used for an error bound on f, which we will accumulate
        // as we go.  The full formula is :
        // e = 2rho + Sum_j=1..k (k-j+6) z_j^2/abs(D_j^2-s^2) 
        //       + Sum_j=k+1..N (j-k+5) z_j^2/(D_j^2-s^2)
        //       + abs(f(s)) + abs(tau)*abs(f'(s))

        // psix, phix store the value of psi,phi without the last term added.
        // These are useful in the threepoles section below.
        T psix=0, phix=0, dpsix=0, dphix=0;

        T f1 = f;
        for(int iter = 0; iter < TMV_MAXITER; iter++) {
            // Calculate psi, phi
            psi = phi = e = dpsi = dphi = T(0);
            for(ptrdiff_t j=0;j<k1;j++) {
                T temp = (z[j] / sum[j]) / diff[j];
                psix = psi; dpsix = dpsi;
                psi += z[j] * temp;
                dpsi += temp * temp;
                e -= psi;
            }
            for(ptrdiff_t j=N-1;j>k1;j--) {
                T temp = (z[j] / sum[j]) / diff[j];
                phix = phi; dphix = dphi;
                phi += z[j] * temp;
                dphi += temp * temp;
                e += phi;
            }

            // Finish calculating f
            fk = rho + psi + phi; // == f(s) without z_k1^2/(D_k1^2-s^2) term
            dfk = dpsi + dphi;
            T temp = (z[k1] / sum[k1]) / diff[k1];
            k1term = temp*z[k1];
            dk1term = temp*temp;
            f = fk + k1term;
            df = dfk + dk1term;
            e += T(2)*rho + T(6)*(phi-psi + TMV_ABS(k1term)) +
                TMV_ABS(f) + TMV_ABS(tau)*df;
            dbgcout<<"rho = "<<rho<<endl;
            dbgcout<<"psi = "<<psi<<endl;
            dbgcout<<"phi = "<<phi<<endl;
            dbgcout<<"fk = "<<fk<<endl;
            dbgcout<<"k1 = "<<k1<<"  k = "<<k<<endl;
            dbgcout<<"sum = "<<sum[k1]<<endl;
            dbgcout<<"diff = "<<diff[k1]<<endl;
            dbgcout<<"z = "<<z[k1]<<endl;
            dbgcout<<"temp = "<<temp<<endl;
            dbgcout<<"k1term = "<<k1term<<endl;
            dbgcout<<"f = "<<f<<endl;
            dbgcout<<"df = "<<df<<endl;

            if (k == N-1) break;
            //if (TMV_ABS(f) <= eps*e) break; 
            // Do it this way instead so nan will break out too.
            if (!(TMV_ABS(f) > eps*e)) break; 
            if (TMV_ABS(f) <= RT(1.e3)*TMV_ABS(f1) && f*f1 < T(0)) break;

            // If we get to here, then the initial guess was not very good.  
            // This can mean that the function is not very well described
            // by the nearby poles.
            // So use bisection until we at least do better than the 
            // midpoint guess we started with.
            // Then we can move on to the fancy, quick-converging techniques 
            // below.
            // This is quite rare, but when the above calculation gives a very 
            // bad initial guess, this adjustment can be critical.
            dbgcout<<"f = "<<f<<", f1 = "<<f1<<endl;
            T dt; // = d(tau)
            if (f*f1 < T(0)) {
                // Then we are too close to the pole.
                // Try doubling the distance from the pole.
                if (f < T(0)) {
                    dbgcout<<"Too close to left pole\n";
                    dt = std::pow(2.,iter)*(tau-lowerbound);
                    if (tau+dt > upperbound) dt = (upperbound - tau)/T(2);
                    lowerbound = tau;
                } else {
                    dbgcout<<"Too close to right pole\n";
                    dt = std::pow(2.,iter)*(tau-upperbound);
                    if (tau+dt < lowerbound) dt = (lowerbound - tau)/T(2);
                    upperbound = tau;
                }
            } else {
                dbgcout<<"Use bisection to update initial tau:\n";
                if (f < T(0)) { lowerbound = tau; dt = (upperbound-tau)/T(2); }
                else { upperbound = tau; dt = (lowerbound-tau)/T(2); }
            }
            dbgcout<<"dt = "<<dt<<endl;
            tau += dt;
            dbgcout<<"tau -> "<<tau<<endl;
            dbgcout<<"bounds -> "<<lowerbound<<"  "<<upperbound<<endl;

            // change dt from d(tau) to d(t)
            dt /= s + TMV_SQRT(s*s+dt);
            t += dt;
            s += dt;
            for(ptrdiff_t j=0;j<N;j++) sum[j] = D[j]+s;
            if (TMV_ABS(s) < TMV_ABS(dt)) {
                dbgcout<<"Need to redo the diff calculations\n";
                for(ptrdiff_t j=0;j<N;j++) diff[j] = D[j]-s; 
            } else  {
                for(ptrdiff_t j=0;j<N;j++) diff[j] -= dt; 
            }

            if (iter == TMV_MAXITER-1) {
                TMV_Warning(
                    "Warning - Unable to find appropriate initial "
                    "guess in FindDCSingularValue\n"
                    "Unsuccessful after 100 iterations.\n");
            }
        }

        // Use 3 poles when fk<0 for k1=k or fk>0 for k1=k+1
        bool threepoles = (((fk < T(0)) == (k1 == k)) && k1>0 && k1<N-1);
        bool fixedweight = true;
        if (threepoles) dbgcout<<"USE THREE POLES\n";

        bool last = false;
        for(int iter = 0; iter < TMV_MAXITER; iter++) {
            dbgcout<<"k = "<<k<<", loop 2: iter = "<<iter<<std::endl;
            dbgcout<<"Main iter = "<<iter<<", f = "<<f<<
                ", eps*e = "<<eps*e<<endl;
            dbgcout<<"f("<<tau<<") = "<<f<<"  s = "<<s<<
                ", s^2 = "<<D[k1]*D[k1]+tau<<endl;
            dbgcout<<"df = "<<df<<", eta_newt = -f/df = "<<-f/df<<std::endl;

            T eta; // = the change to be added to s

            if (!(TMV_ABS(f) > eps*e)) {
                last = true;
                eta = T(0); // This will trigger a Newton step below.
            } else if (!threepoles) {
                // The normal case: only 2 poles used for approx
                T d1 = diff[kk] * sum[kk];
                T d2 = diff[kk+1] * sum[kk+1];
                dbgcout<<"d1,d2 = "<<d1<<", "<<d2<<endl;
                if (fixedweight) 
                    if (kk==k1) c = fk - d2*dfk;
                    else c = fk - d1*dfk;
                else c = fk - d1*dpsi - d2*dphi;
                // c eta^2 - a eta + b = 0
                T d1x = d1/d2;
                T a = (d1x+T(1))*f - d1*df;
                T b = d1x*f;
                dbgcout<<"c,a,b = "<<c<<", "<<a<<", "<<b<<endl;
                if (c == T(0)) {
                    if (a == T(0)) {
                        // Then the above has rounding error.
                        // This version can't have a == 0:
                        if (fixedweight)
                            // For k1 == k case:
                            // a = (d1+d2)*f - d1*d2*df
                            // with f = d2*dfk + k1term  (i.e. c==0)
                            // a = (d1+d2)*d2*dfk + (d1+d2)*k1term 
                            //       - d1*d2*dfk - d1*d2*dk1term
                            //   = d2^2*dfk + (d1+d2)*k1term - d1*d2*(k1term/d1)
                            //   = d2^2*dfk + d1 zsq(k1)/d1
                            // For k1 == k+1 case, swap d1,d2
                            a = zsq[k1]/d2 + (k1==kk ? d2 : d1*d1x)*dfk;
                        else
                            // a = (d1+d2)*f - d1*d2*df
                            // with f = d1*dpsi + d2*dphi + k1term
                            // a = (d1+d2)*d1*dpsi + (d1+d2)*d2*dphi - 
                            //      - d1*d2*(dpsi+dphi) 
                            //      + (d1+d2) k1term - d1*d2*dk1term
                            //   = d1*d1*dpsi + d2*d2*dphi + d1 k1term
                            a = d1x*d1*dpsi + d2*dphi + 
                                (k1==kk ? d1x : T(1))*k1term;
                    }
                    eta = b/a;
                    eta *= d2;
                    dbgcout<<"c=0  eta = "<<eta<<endl;
                    if (k==N-1)
                        TMVAssert(eta > -tau);
                    else
                        TMVAssert(eta > d1 && eta < d2);
                } else {
                    T d = a*a-T(4)*b*c;
                    dbgcout<<"d = "<<d<<endl;
                    if (d < T(0)) d = T(0); 
                    else d = TMV_SQRT(d);
                    if (k == N-1) {
                        // Then we want the solution with other sign for d
                        
                        // Actually, now that I am scaling a,b by d2 and d2^2,
                        // this changes things, since d2 is negative.
                        // So now we actually choose the same solution as
                        // in the normal case.  The only difference in this
                        // section now is the bounds test.
                        eta = a <= 0 ? (a-d)/(T(2)*c) : T(2)*b/(a+d);
                        eta *= d2;
                        dbgcout<<"eta = "<<eta<<endl;
                        if (eta <= -tau) // Then rounding errors - use bisect
                            eta = -tau/RT(2);
                        TMVAssert(eta > -tau);
                    } else {
                        eta = a <= 0 ? (a-d)/(T(2)*c) : T(2)*b/(a+d);
                        eta *= d2;
                        dbgcout<<"eta = "<<eta<<endl;
                        if (eta <= d1) // Then rounding errors - use bisect
                            eta = d1/RT(2);
                        if (eta >= d2) // Then rounding errors - use bisect
                            eta = d2/RT(2);
                        dbgcout<<"eta => "<<eta<<endl;
                        TMVAssert(eta > d1 && eta < d2);
                    }
                }
            } else {
                // More complicated case using 3 poles. 
                dbgcout<<"threepoles section\n";
                T d1 = diff[k1-1] * sum[k1-1];
                T d2 = diff[k1+1] * sum[k1+1];
                dbgcout<<"d1,d2 = "<<d1<<" "<<d2<<endl;
                if (fixedweight) 
                    if (k1 == k) c = (rho+psix+phi) - d2*(dpsix+dphi);
                    else  c = (rho+psi+phix) - d1*(dpsi+dphix);
                else 
                    if (k1 == k) c = (rho+psix+phi) - d1*dpsix - d2*dphi;
                    else c = (rho+psi+phix) - d1*dpsi - d2*dphix;
                dbgcout<<"c = "<<c<<endl;
                // f(eta) is modeled as:
                // f(eta) = c + s1/(d1-eta) + z^2[k1]/(-tau-eta) + s2/(d2-eta)
                // given that f(0) = f
                // and f'(0) = df
                // These lead to: c + s1/d1 + s2/d2 = fk
                // and  c (d1+d2) + s1 + s2 = fk (d1+d2) - dfk d1 d2
                //      
                // We solve f(eta) = 0 as:
                // c (d1-eta)(tau+eta)(d2-eta) + s1(tau+eta)(d2-eta)
                //       - z^2[k1](d1-eta)(d2-eta) + s2(d1-eta)(tau+eta) = 0
                // using the above equations, we remove s1,s2:
                // c eta^3 + (c tau - fk (d1+d2) + dfk d1 d2 + tau k1term) eta^2
                //    + (fk d1 d2 - tau f (d1+d2) + tau dfk d1 d2) eta
                //    + tau d1 d2 f = 0
                // So, cubic equation is 
                // h(eta) = c eta^3 - a eta^2 + b eta + g = 0
                // with a = {fk(d1+d2) - dfk d1 d2} + z^2[k1] - c tau
                //      b = d1 d2 (fk + tau dfk) - f (d1+d2) tau 
                //      b = d1 d2 fk - tau ( 
                //            {fk(d1+d2) - dfk d1 d2} + (d1+d2) k1term)
                //      g = d1 d2 tau f
                // We need to be careful with the terms in {} above 
                // (in a and b), since the terms can nearly cancel.
                // Do better with phix,psix:
                // For k1 == k:
                //      {} = (rho + psix + zsq[k1-1]/d1 + phi)(d1+d2) 
                //           - (dpsix + zsq[k1-1]/d1^2 + dphi)(d1 d2)
                //         = (rho+psix+phi)(d1+d2) - (dpsix+dphi)(d1 d2) 
                //           + zsq[k1-1]
                //
                // Also, to help avoid problems with underflow, we effectively
                // measure eta in units of d2.
                // So a /= d2, b /= d2^2, g /= d2^3
                T d1x = d1/d2;
                T taux = tau/d2;
                T temp = (
                    k1==k ? 
                    (d1x+T(1))*(rho+psix+phi) - d1*(dpsix+dphi) + zsq[k1-1]/d2 :
                    (d1x+T(1))*(rho+psi+phix) - d1*(dpsi+dphix) + zsq[k1+1]/d2);
                T a = temp + zsq[k1]/d2 - c*taux;
                T b = d1x*fk - taux*(temp + (d1x+T(1))*k1term);
                T g = d1x*taux*f; 
                T mineta = (lowerbound - tau)/d2;
                T maxeta = (upperbound - tau)/d2;
                dbgcout<<"temp = "<<temp<<endl;
                dbgcout<<"a = "<<a<<endl;
                dbgcout<<"b = "<<b<<endl;
                dbgcout<<"g = "<<g<<endl;
                dbgcout<<"mineta, maxeta = "<<mineta<<", "<<maxeta<<endl;

                // Bounds on eta are:
                // lb < tau + eta < ub 
                // lb - tau < eta < ub - tau
                // h(0) = g > 0 regardless of the sign of f.
                // if f > 0, then deta needs to be negative
                // if f < 0, then deta needs to be positive
                // This dictates which bound is the one to use for eta2
                // In below iteration, we maintain h(eta2) * h(eta) < 0
                T eta2 = (f>0) ? mineta : maxeta;
                dbgcout<<"eta2 = "<<eta2<<endl;
                T h2 = ((c*eta2 - a)*eta2 + b)*eta2 + g;

                eta = T(0);
                T h = g; // = h(0)

                if (h*h2 > T(0)) {
                    // Normally h and h2 are opposite sign.
                    // if they are the same sign, then that means that the 
                    // threepoles approximation is the same sign over the
                    // entire allowed bounds for eta.
                    // Since threepoles isn't a good approximation, use 
                    // Newton instead.
                    dbgcout<<"h(0) = g = "<<h<<endl;
                    dbgcout<<"h(eta2) = h2 = "<<h2<<endl;
                    dbgcout<<"Same sign, so use Newton\n";
                    eta = (-f/df)/d2;
                } else for(int iter3 = 0; iter3 < TMV_MAXITER; iter3++) {
                    dbgcout<<"iter3 = "<<iter3<<", eta = "<<eta<<
                        "  h = "<<h<<endl;
                    dbgcout<<"eta2 = "<<eta2<<"  h2 = "<<h2<<endl;
                    // We do something sort of like a Newton step, but use
                    // 2nd order info as well, since it's so cheap to calculate
                    // and you can get a big improvement with it.
                    // h(eta) = c eta^3 - a eta^2 + b eta + g = 0
                    // dh/deta = 3 c eta^2 - 2 a eta + b
                    // d2h/deta2 = 6 c eta - 2 a
                    // Local quadratic is
                    // h~(de) = h(eta) + dh/deta de + 1/2 d2h/deta2 de^2
                    // Solve h~(de) = 0
                    // de = (-dh +- sqrt(dh^2 - 2 h d2h)) / d2h
                    // de = -2h / (dh -+ sqrt(dh^2 - 2 h d2h))
                    // Use:
                    //   de = -2h / dh(1 + sqrt(1 - 2 h d2h/dh^2))  (first)
                    //   de = -dh(1 + sqrt(1 - 2 h d2h/dh^2)) / d2h (second)
                    // The first choice is intended to be closer to the
                    // Newton solution when h is very small.
                    T dh = (T(3)*c*eta - T(2)*a)*eta + b;
                    T d2h = T(6)*c*eta - T(2)*a;
                    T x = T(1) - T(2)*(h/dh)*(d2h/dh);
                    dbgcout<<"dh,d2h,x = "<<dh<<", "<<d2h<<", "<<x<<endl;
                    x = (x <= T(0)) ? T(0) : TMV_SQRT(x);
                    T deta = -T(2)*h / (dh*(T(1) + x));
                    dbgcout<<"deta = "<<deta<<endl;
                    if ((eta2 > eta && (deta < 0 || deta > eta2-eta)) || 
                        (eta2 < eta && (deta > 0 || deta < eta2-eta)) ) {
                        dbgcout<<"deta is wrong direction or too big - "
                            "use other root.\n";
                        deta = -dh*(T(1)+x) / d2h;
                        dbgcout<<"deta = "<<deta<<endl;
                        if ((eta2 > eta && (deta < 0 || deta > eta2-eta)) || 
                            (eta2 < eta && (deta > 0 || deta < eta2-eta)) ) {
#if 1
                            // This seems to work a bit better when 
                            // adjacent poles differ by orders of magnitude.
                            // I'm not sure if it is always the better
                            // choice in general.
                            dbgcout<<"deta is wrong direction or too big - "
                                "use Newton step and break out of 3-poles.\n";
                            eta = (-f/df)/d2;
                            break;
#else
                            dbgcout<<"deta is wrong direction or too big - "
                                "use bisection.\n";
                            deta = (eta2-eta)/T(2);
                            dbgcout<<"deta -> "<<deta<<endl;
#endif
                        }
                    }
                    TMVAssert(eta + deta <= maxeta);
                    TMVAssert(eta + deta >= mineta);
                    T etanew = eta + deta;
                    T hnew = 
                        (TMV_ABS(deta) < RT(1.e-3)*TMV_ABS(eta)) ?
                        // If change in eta is small compared to eta, then
                        // do delta calculation for hnew
                        // c(e+de)^3 - a(e+de)^2 + b(e+de) + g
                        // h + 3ce^2 de + 3ce de^2 + cde^3 - 2ae de - a de^2 
                        //     + b de
                        // h + c de^3 + (3ce-a) de^2 + (3ce^2-2ae+b) de
                        h + ((c*deta + T(3)*c*eta-a)*deta + 
                             (T(3)*c*eta-T(2)*a)*eta+b )*deta :
                        // Else regular calculation is ok.
                        ((c*etanew - a)*etanew + b)*etanew + g;
                    dbgcout<<"etanew = "<<etanew<<endl;
                    dbgcout<<"hnew = "<<hnew<<endl;
                    if ( (h > T(0)) != (hnew > T(0)) ) { eta2 = eta; h2 = h; }
                    eta = etanew;
                    h = hnew;
                    if (TMV_ABS(deta) < eps*TMV_ABS(eta) || 
                        TMV_ABS(h) < eps*g) break;
                    if (iter3 == TMV_MAXITER-1) {
                        TMV_Warning(
                            "Warning - unable to converge for THREEPOLES "
                            "solution in FindDCSingularValue\n"
                            "No solution after 100 iterations.\n");
                    }
                }
                eta *= d2;
                dbgcout<<"rescaled eta => "<<eta<<endl;
            }
            dbgcout<<"eta = "<<eta<<endl;

            // Check for eta being the wrong sign.
            // If f < 0, then eta should be > 0.
            // If f > 0, then eta should be < 0.
            // If not, use a Newton step instead for which this is guaranteed.
            if (f*eta >= T(0)) {
                eta = -f/df;
                dbgcout<<"Newton eta = "<<eta<<endl;
            }

            // Also check that tau+eta is still within allowed bounds:
            if (tau+eta < lowerbound) {
                eta = (lowerbound-tau)/T(2);
                dbgcout<<"halfway to lowerbound eta = "<<eta<<endl;
            } else if (tau+eta > upperbound) {
                eta = (upperbound-tau)/T(2);
                dbgcout<<"halfway to upperbound eta = "<<eta<<endl;
            }

            tau += eta;
            TMVAssert(tau >= lowerbound && tau <= upperbound);
            // Update to t: 
            // d^2+tau = (d+t)^2 = d^2 + 2td + t^2
            // tau = 2td + t^2
            // tau + eta = 2(t+dt)d + (t+dt)^2
            // eta = 2 d dt + 2 t dt + dt^2
            // dt = eta / ((d+t) + sqrt((d+t)^2+eta) )
            eta /= s + TMV_SQRT(s*s+eta);
            t += eta;
            s += eta;
            for(ptrdiff_t j=0;j<N;j++) sum[j] = D[j]+s;
            if (TMV_ABS(s) < TMV_ABS(eta)) {
                // Then the iterative adjustment to diff isn't going to 
                // maintain the accuracy we need.
                dbgcout<<"Need to redo the diff calculations\n";
                for(ptrdiff_t j=0;j<N;j++) diff[j] = D[j]-s; 
            } else {
                for(ptrdiff_t j=0;j<N;j++) diff[j] -= eta; 
            }

            if (last) break;

            // Update psi, phi, etc.
            psi = T(0);
            dpsi = T(0);
            e = T(0);
            for(ptrdiff_t j=0;j<k1;j++) {
                T temp = (z[j] / sum[j]) / diff[j];
                psix = psi; dpsix = dpsi;
                psi += z[j] * temp;
                dpsi += temp * temp;
                e -= psi;
            }
            phi = T(0);
            dphi = T(0);
            for(ptrdiff_t j=N-1;j>k1;j--) {
                T temp = (z[j] / sum[j]) / diff[j];
                phix = phi; dphix = dphi;
                phi += z[j] * temp;
                dphi += temp * temp;
                e += phi;
            }
            fk = rho + psi + phi; 
            dbgcout<<"fk = "<<rho<<" + "<<psi<<" + "<<phi<<std::endl;
            dfk = dpsi+dphi;
            T temp = (z[k1] / sum[k1]) / diff[k1];
            k1term = temp*z[k1];
            dk1term = temp*temp;
            T fnew = fk + k1term;
            dbgcout<<"fnew = "<<fk<<" + "<<k1term<<" = "<<fnew<<std::endl;
            df = dfk + dk1term;
            e += T(2)*rho + T(6)*(phi-psi + TMV_ABS(k1term)) + 
                TMV_ABS(fnew) + TMV_ABS(tau)*df;

            // Use 3 poles when fk<0 for k1=k or fk>0 for k1=k+1
            bool newthreepoles = (((fk < T(0)) == (k1 == k)) && k1>0 && k1<N-1);
            if (newthreepoles != threepoles) {
                dbgcout<<"Change threepoles from "<<threepoles<<" to "<<
                    newthreepoles<<endl;
                dbgcout<<"fk = "<<fk<<"  fnew = "<<fnew<<endl;
                threepoles = newthreepoles;
            }

            // Update bounds:
            if (fnew < T(0)) lowerbound = tau;
            else upperbound = tau;
            dbgcout<<"New bounds = "<<lowerbound<<"  "<<upperbound<<endl;

            // Update scheme:
            if (fnew*f > T(0) && TMV_ABS(fnew) > T(0.1) * TMV_ABS(f) && k<N-1)
                fixedweight = !fixedweight;

            f = fnew;
            if (iter == TMV_MAXITER-1) {
                TMV_Warning(
                    "Warning - Unable to find solution in "
                    "FindDCSingularValue\n"
                    "No solution after 100 iterations.\n");
            }
        }
        dbgcout<<"Found Singularvalue S("<<k<<") = "<<s<<endl;
        dbgcout<<"eta = "<<-f/df<<", tau = "<<tau<<
            ", eta/tau = "<<-f/df/tau<<std::endl;
#ifdef XDEBUG
        // f(s) = rho + Sum_i=1..N z_i^2/(D_i^2-s^2)
        T ff = rho;
        dbgcout<<"f(s) = "<<rho<<" + ";
        for(ptrdiff_t j=0;j<N;j++) {
            dbgcout<<(z[j]/diff[j])*(z[j]/sum[j])<<" + ";
            ff += (z[j]/diff[j])*(z[j]/sum[j]);
        }
        dbgcout<<" = "<<ff<<std::endl;
        dbgcout<<"cf f = "<<f<<std::endl;
#endif
        return s;
    }

    template <class T> 
    void FindDCSingularValues(
        Vector<T>& S, T rho, const GenVector<T>& D, const GenVector<T>& z, 
        Matrix<T,ColMajor>& diffmat)
    {
        dbgcout<<"Start FindDCSV: "<<std::endl;
        dbgcout<<"D = "<<D<<endl;
        dbgcout<<"z = "<<z<<endl;
        dbgcout<<"rho = "<<rho<<endl;
        TMVAssert(rho > T(0));
        TMVAssert(S.size() == D.size());
        TMVAssert(S.size() == z.size());
        TMVAssert(S.size() == diffmat.colsize());
        TMVAssert(S.size() == diffmat.rowsize());

        const ptrdiff_t N = S.size();
        Vector<T> zsq(N);
        for(ptrdiff_t j=0;j<N;j++) zsq[j] = z[j]*z[j];
        const T normsqz = zsq.sumElements();

#ifdef _OPENMP
#pragma omp parallel 
        {
            Vector<T> diff(N);
            Vector<T> sum(N);
            T Sk;
            // For some reason pgCC requires an actual int for the omp
            // for loop, not just a signed integer type.  So ptrdiff_t
            // can't be used.
#ifdef __PGI
#define TMV_INT_OMP int
#else
#define TMV_INT_OMP ptrdiff_t
#endif
#pragma omp for
            for(TMV_INT_OMP k=0;k<N;k++) {
                Sk = FindDCSingularValue(
                    k,N,rho,D.cptr(),z.cptr(),
                    zsq.cptr(),normsqz,diff.ptr(),sum.ptr());
#pragma omp critical
                {
                    S[k] = Sk;
                    diffmat.col(k) = diff;
                }
            }
        }
#else
        Vector<T> sum(N);
        for(ptrdiff_t k=0;k<N;k++) {
            S[k] = FindDCSingularValue(
                k,N,rho,D.cptr(),z.cptr(),
                zsq.cptr(),normsqz,diffmat.col(k).ptr(),sum.ptr());
        }
#endif
        dbgcout<<"S => "<<S<<endl;
    }

    template <class T> 
    void FindDCSingularValues(
        Vector<T>& S, T rho, const GenVector<T>& D, const GenVector<T>& z)
    {  
        dbgcout<<"Start FindDCSV (No diffmat): "<<std::endl;
        dbgcout<<"D = "<<D<<endl;
        dbgcout<<"z = "<<z<<endl;
        TMVAssert(rho > T(0));
        TMVAssert(S.size() == D.size());
        TMVAssert(S.size() == z.size());

        const ptrdiff_t N = S.size();
        Vector<T> zsq(N);
        for(ptrdiff_t j=0;j<N;j++) zsq[j] = z[j]*z[j];
        T normsqz = zsq.sumElements();

#ifdef _OPENMP
#pragma omp parallel
        {
            Vector<T> diff(N);
            Vector<T> sum(N);
            T Sk;
#pragma omp for
            for(TMV_INT_OMP k=0;k<N;k++) {
                Sk = FindDCSingularValue(
                    k,N,rho,D.cptr(),z.cptr(),
                    zsq.cptr(),normsqz,diff.ptr(),sum.ptr());
#pragma omp critical
                {
                    S[k] = Sk;
                }
            }
#undef TMV_INT_OMP
        }
#else
        Vector<T> diff(N);
        Vector<T> sum(N);
        for(ptrdiff_t k=0;k<N;k++) {
            S[k] = FindDCSingularValue(
                k,N,rho,D.cptr(),z.cptr(),
                zsq.cptr(),normsqz,diff.ptr(),sum.ptr());
        }
#endif
        dbgcout<<"S => "<<S<<endl;
    }

    template <class T> 
    static void SmallProblem(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E, MatrixView<T> Vt)
    {
        dbgcout<<"Start SmallProblem: N = "<<D.size()<<std::endl;
        if (E.size() > 0)
            SV_DecomposeFromBidiagonal_QR(U,D,E,Vt);
        dbgcout<<"After QR"<<std::endl;
        // Make all of the singular values positive
        const ptrdiff_t N = D.size();
        RT* Di = D.ptr();
        for(ptrdiff_t i=0;i<N;++i,++Di) if (*Di < 0) {
            *Di = -(*Di);
            if (Vt.cptr()) Vt.row(i) = -Vt.row(i);
        }
        dbgcout<<"After make all Di positive"<<std::endl;
    }

    template <class T> 
    void SV_DecomposeFromBidiagonal_DC(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E, MatrixView<T> Vt,
        bool UisI, bool VisI)
    {
        // For these comments, I will use FortranStyle notation,
        // wherein the first element of D is D_1, and the last is D_N.
        // This notation just happens to be a bit more natural for
        // discussing this algorithm.
        //
        // Solve the SVD of unreduced Bidiagonal Matrix B (given by D,E).
        // Note: the input B must be unreduced - ie. all entries are non-zero.
        // This routine implements the divide and conquer approach, calling
        // itself for the recursion, and the QR algorithm when the
        // size gets too small for divide and conquer to be efficient.
        //
        // The basic idea of the divide and conquer algorithm is to split
        // the bidiagonal matrix B into two parts, B1 and B2 and a joining
        // row:
        //
        //     [    B1       0    ] K-1
        // B = [ D_K e_K  E_k e_1 ] 1
        //     [    0       B2    ] N-K
        //          K      N-K
        //
        // where e_k is a basis vector with all 0's, except for 1 at element k.
        //
        // The smaller bidiagonal matrices are first reduced recursively:
        // B1 = U1 (S1 0) V1t
        // B2 = U2 S2 V2t
        //
        //
        //     [ U1 0  0  ] [ (S1 0)  0  ] [ V1T  0 ]
        // B = [ 0  1  0  ] [   z1    z2 ] [ 0  V2T ]
        //     [ 0  0  U2 ] [   0     S2 ] 
        //
        // The vector z = [ z1 z2 ] is such that z VT = [ 0 D_K E_K 0 ]
        // which implies that:
        // z1 V1T = D_K e_K
        // z1 = D_K e_K V1
        //    = D_K V1.row(K)
        //    = D_K V1T.col(K)
        //
        // z2 V2T = E_K e_1
        // z2 = E_K e_1 V2
        //    = E_K V2.row(1)
        //    = E_K V2T.col(1)
        //
        // Now, all that remains is to find the SVD of the central matrix.
        // It turns out that this is not too difficult.
        // The first step is to permute the matrix into:
        //
        //     [ z1 z2 z3 ... zN ]
        //     [    d2           ]
        // M = [       d3        ] = D + e1 zT
        //     [          ...    ]
        //     [              dN ]
        //
        // where 0 < d2 < d3 < d4 < d5 < .... < dN,
        // We define d1 == 0 for simplicity in the following discussion.
        // We also define D = diag(d_i), z = {z_i}.
        //
        // The singular values of the above matrix are the solutions of 
        // the so-called secular equation:
        // f(s) = det(MTM - s^2I) = 0
        // 
        // MTM = (D + z e1T) (D + e1 zT) = D^2 + zzT
        // where we used the fact that d1 == 0.
        //
        // It turns out that if the d_i are all distinct, then
        // f(s) = 1 + zT (D^2 - s^2 I)^-1 z 
        //      = 1 + Sum_i zi^2 / (d_i^2-s^2)
        //
        // Furthermore, it can be shown that the solutions to f(s) = 0 are
        // such that:
        // 0 = d1 < s1 < d2 < s2 < d3 < s3 < ... < dN < sN < dN+|z|
        //
        // and that M = X S Y
        // 
        // with S = diag(s_i)
        // y_i = (D^2 - si^2 I)^-1 z
        // Y.row(i) = y_i / Norm(y_i);
        // x_i = ( 1 , d2 z2/(d2^2-si^2) , ... dN zN/(dN^2-si^2) )
        // X.col(i) = x_i / Norm(x_i);
        //
        // This then gives us the complete solution to our SVD.
        // However, there are still a couple of wrinkes to deal with.
        //
        // First, the solution above requires that the z_i are non-zero and
        // that the d_i are distinct.
        //
        // Zero-valued z_i's are easy to deal with.  If z_k is zero (or
        // in practice less than EPS*(d_max+|z|) ) then d_k is a singular value
        // itself, so you can permute the matrix to:
        //
        //     [ z1 z2 z3 ... zN-1 0 ]
        //     [    d2               ]
        // M = [       d3            ]
        //     [          ...        ]
        //     [              dN-1   ]
        //     [                  dN ]
        // 
        // So the now dN is a singular value and the corresponding singular 
        // vectors are just e_N.
        //
        // If there are two d_i which are very close to each other (difference
        // < EPS*(d_max+|z|), we need to use a Givens rotation to rotate
        // z_j into z_i.  Doing this on the rows as well preserves the two d's,
        // but now we have a 0 in the z vector, and we can "deflate" it as
        // described above.
        //
        // The second concern that these same authors address is that it is
        // numerically more stable to calculate the s values as above, 
        // but then to recalculate the z_i as:
        //
        // z_i' = sqrt( (sn^2-di^2) prod_1..i-1 (sk^2-di^2)/(dk^2-di^2) *
        //             prod_i..N-1 (sk^2-di^2)/(dk+1^2-di^2) ) * sign(z_i)
        //
        // and then using these z_i' values for the calculation of x_i, y_i 
        // as described above.  I haven't really tried to understand why this 
        // is the case, but I take their word for it and use this method.
        //

        TMVAssert(D.size()>0);
        TMVAssert(E.size()+1 == D.size());
        if (U.cptr()) TMVAssert(U.rowsize() == D.size()); 
        if (Vt.cptr()) TMVAssert(Vt.colsize() == D.size()); 
        TMVAssert(D.step()==1);
        TMVAssert(E.step()==1);

#ifdef XDEBUG
        dbgcout<<"Start SVD DC: \n";
        Vector<RT> D0 = D;
        dbgcout<<"D0 = "<<D0<<endl;
        Vector<RT> E0 = E;
        dbgcout<<"E0 = "<<E0<<endl;
        Matrix<T> U0(U.cptr()?U.colsize():D.size(),D.size());
        if (U.cptr()) U0 = U; else U0.setToIdentity();
        Matrix<T> Vt0(D.size(),Vt.cptr()?Vt.rowsize():D.size());
        if (Vt.cptr()) Vt0 = Vt; else Vt0.setToIdentity();
        Matrix<T> A0(U0.colsize(),Vt0.rowsize());
        Matrix<T> B0(D.size(),D.size(),RT(0));
        B0.diag() = D;
        if (E0.size() > 0) B0.diag(1) = E;
        A0 = U0 * B0 * Vt0;
        //dbgcout<<"A0 = "<<A0<<endl;
        double normA0 = Norm(A0);
        if (normA0 < 1.) normA0 = 1.;
        dbgcout<<"Norm(A0) = "<<Norm(A0)<<endl;
        if (U.cptr()) dbgcout<<"Norm(UtU -1) = "<<Norm(U0.adjoint()*U0-T(1))<<endl;
        else dbgcout<<"No input U\n";
        if (Vt.cptr()) dbgcout<<"Norm(VtVtt -1) = "<<Norm(Vt0*Vt0.adjoint()-T(1))<<endl;
        else dbgcout<<"No input Vt\n";
#endif

        ptrdiff_t N = D.size();
        // If N is too small, use the QR method
        if (N <= DC_LIMIT) SmallProblem(U,D,E,Vt);

        // It seems that the QR is always faster than DC for N < ~ 5000
        // And when N is much larger than this, then we start to worry about
        // memory issues of making the U and Vt matrices for the subproblems,
        // and QR is certainly never very much slower than doing 1 or 2 divides,
        // which seems to be optimal for large problems.
        // So I always use QR for the S only calculation.
        else if (!(U.cptr() || Vt.cptr())) SmallProblem(U,D,E,Vt);

        else {
            dbgcout<<"N > "<<DC_LIMIT<<endl;
            ptrdiff_t K = (N-1)/2;
            const RT DK = D(K);
            const RT EK = E(K);
            Vector<RT> z(N,RT(0));

            // Do the left sub-problem
            VectorView<RT> D1 = D.subVector(0,K);
            VectorView<RT> E1 = E.subVector(0,K);
            Matrix<RT,RowMajor> Vt1(K+1,Vt.cptr()?K+1:1);
            if (Vt.cptr()) Vt1.setToIdentity();
            else Vt1.col(0).makeBasis(K); // only need col(K)
            BidiagonalZeroLastCol<RT>(D1,E1,Vt1.view());
            if (U.cptr()) {
                Matrix<RT,ColMajor> U1(K,K);
                U1.setToIdentity();
                SV_DecomposeFromBidiagonal_DC<RT>(
                    U1.view(),D1,
                    E1.subVector(0,K-1),Vt1.rowRange(0,K),true,false);
                if (UisI) U.subMatrix(0,K,0,K) = U1;
                else U.colRange(0,K) *= U1;
            } else {
                MatrixView<RT> U1(0,0,0,1,1,NonConj);
                SV_DecomposeFromBidiagonal_DC<RT>(
                    U1,D1,E1.subVector(0,K-1),
                    Vt1.rowRange(0,K),false,false);
            }
            z.subVector(0,K+1) = DK * Vt1.col(Vt.cptr() ? K : 0);
            if (Vt.cptr()) {
                if (VisI) Vt.subMatrix(0,K+1,0,K+1) = Vt1;
                else Vt.rowRange(0,K+1) = Vt1 * Vt.rowRange(0,K+1);
            }

            // Do the right sub-problem
            VectorView<RT> D2 = D.subVector(K+1,N);
            VectorView<RT> E2 = E.subVector(K+1,N-1);
            Matrix<RT,RowMajor> Vt2(N-K-1,Vt.cptr()?N-K-1:1);
            if (Vt.cptr()) Vt2.setToIdentity();
            else Vt2.col(0).makeBasis(0); // only need col(0)
            if (U.cptr()) {
                Matrix<RT,ColMajor> U2(N-K-1,N-K-1);
                U2.setToIdentity();
                SV_DecomposeFromBidiagonal_DC<RT>(
                    U2.view(),D2,E2,Vt2.view(),true,Vt.cptr());
                if (UisI) U.subMatrix(K+1,N,K+1,N) = U2;
                else U.colRange(K+1,N) *= U2;
            } else {
                MatrixView<RT> U1(0,0,0,1,1,NonConj);
                SV_DecomposeFromBidiagonal_DC<RT>(
                    U1,D2,E2,Vt2.view(),false,Vt.cptr());
            }
            z.subVector(K+1,N) = EK * Vt2.col(0);
            D(K) = RT(0);
            if (Vt.cptr()) {
                if (VisI) Vt.subMatrix(K+1,N,K+1,N) = Vt2;
                else Vt.rowRange(K+1,N) = Vt2 * Vt.rowRange(K+1,N);
            }
#ifdef XDEBUG
            dbgcout<<"Done subproblems (N="<<N<<")\n";
            dbgcout<<"D = "<<D<<endl;
            dbgcout<<"z = "<<z<<endl;
            Matrix<RT> M(N,N,RT(0));
            M.diag() = D;
            M.row(K) = z;
            if (U.cptr() && Vt.cptr()) {
                //dbgcout<<"M = "<<M<<endl;
                //dbgcout<<"UMVt = "<<U * M * Vt<<endl;
                dbgcout<<"Norm(UMVt-A0) = "<<Norm(U * M * Vt -A0)<<endl;
                dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<endl;
                dbgcout<<"Norm(VtV-1) = "<<Norm(Vt*Vt.adjoint()-T(1))<<endl;
                if (Norm(U*M*Vt-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                if (Norm(Vt*Vt.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
            }
#endif

            // Check for deflation step 1:
            const RT eps = TMV_Epsilon<T>();
            RT norminf = D.subVector(0,N).normInf();
            RT tol = eps*norminf;
            dbgcout<<"tol = "<<tol<<endl;
            dbgcout<<"Initially D("<<K<<") = "<<D(K)<<std::endl;
            ptrdiff_t pivot = K;
            for(ptrdiff_t i=0;i<N-1;i++) {
                dbgcout<<"D("<<i<<") = "<<D(i)<<std::endl;
                dbgcout<<"z("<<i<<") = "<<z(i)<<std::endl;
                dbgcout<<"pivot = "<<pivot<<"  "<<D(pivot)<<std::endl;
                if (i != pivot && // Make sure not to deflate the pivot value.
                    (TMV_ABS(z(i)) < tol || z(i)*z(i)*eps == RT(0))) {
                    dbgcout<<"Deflate this one.\n";
                    z.swap(i,N-1);
                    D.swap(i,N-1);
                    if (U.cptr()) U.swapCols(i,N-1);
                    if (Vt.cptr()) Vt.swapRows(i,N-1);
                    if (pivot == N-1) pivot = i;
#ifdef XDEBUG
                    M.swapCols(i,N-1);
                    M.swapRows(i,N-1);
#endif
                    --i; --N;
                    RT absdn = TMV_ABS(D(N));
                    if (absdn == norminf) {
                        norminf = D.subVector(0,N).normInf();
                        tol = eps*norminf;
                    }
                }
            }
            if (TMV_ABS(z(N-1)) < tol) --N;
#ifdef XDEBUG
            dbgcout<<"After deflation\n";
            dbgcout<<"DN = "<<D.subVector(0,N)<<endl;
            dbgcout<<"zN = "<<z.subVector(0,N)<<endl;
            if (U.cptr() && Vt.cptr()) {
                //dbgcout<<"M = "<<M<<endl;
                //dbgcout<<"UMVt = "<<U * M * Vt<<endl;
                dbgcout<<"Norm(UMVt-A0) = "<<Norm(U * M * Vt -A0)<<endl;
                dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<endl;
                dbgcout<<"Norm(VtV-1) = "<<Norm(Vt*Vt.adjoint()-T(1))<<endl;
                if (Norm(U*M*Vt-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                if (Norm(Vt*Vt.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
            }
#endif

            if (N > 1) {
                // Sort the inner matrix so z is in the first row, D are increasing
                AlignedArray<ptrdiff_t> P(N);
                D.subVector(0,N).sort(P.get(),Ascend);
                z.subVector(0,N).permute(P.get());
                if (U.cptr()) U.colRange(0,N).permuteCols(P.get());
                if (Vt.cptr()) Vt.rowRange(0,N).permuteRows(P.get());
#ifdef XDEBUG
                dbgcout<<"After sort\n";
                dbgcout<<"DN = "<<D.subVector(0,N)<<endl;
                dbgcout<<"zN = "<<z.subVector(0,N)<<endl;
                if (U.cptr() && Vt.cptr()) {
                    M.colRange(0,N).permuteCols(P.get());
                    M.rowRange(0,N).permuteRows(P.get());
                    //dbgcout<<"M = "<<M<<endl;
                    //dbgcout<<"U = "<<U<<endl;
                    //dbgcout<<"Vt = "<<Vt<<endl;
                    //dbgcout<<"UMVt = "<<U * M * Vt<<endl;
                    dbgcout<<"Norm(UMVt-A0) = "<<Norm(U * M * Vt -A0)<<endl;
                    dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<endl;
                    dbgcout<<"Norm(VtV-1) = "<<Norm(Vt*Vt.adjoint()-T(1))<<endl;
                    if (Norm(U*M*Vt-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                    if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                    if (Norm(Vt*Vt.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
                }
#endif

                // Check for deflation step 2:
                ptrdiff_t i_firstswap = N;
                tol = eps*D(N-1);
                dbgcout<<"tol = "<<tol<<std::endl;
                for(ptrdiff_t i=N-1;i>0;--i) {
                    dbgcout<<"D("<<i-1<<") = "<<D(i-1)<<std::endl;
                    dbgcout<<"D("<<i<<") = "<<D(i)<<std::endl;
                    dbgcout<<"diff = "<<TMV_ABS(D(i) - D(i-1))<<std::endl;
                    if (TMV_ABS(D(i) - D(i-1)) < tol) {
                        dbgcout<<"Deflate this pair\n";
                        Givens<RT> G = GivensRotate(z(i-1),z(i));
                        TMVAssert(z(i) == RT(0));
                        if (Vt.cptr()) G.mult(Vt.rowPair(i-1,i));
                        if (i > 1) {
                            if (U.cptr()) G.mult(U.colPair(i-1,i).transpose());
#ifdef XDEBUG
                            G.mult(M.colPair(i-1,i).transpose());
                            G.mult(M.rowPair(i-1,i));
#endif
                        } else {
                            // i == 1 messes up a bit of the symmetry, because
                            // we would rotate into i==0, but that is the row
                            // in the matrix M where z is stored.
                            // Doing a Givens rotation from the left would 
                            // rotate M(0,0) into M(1,0), which would be bad.
                            // However, in this case, D(0) == 0, so this 
                            // deflation means that D(1) ~= 0 too.
                            // Therefore, the rotation from the right is 
                            // all you need.
#ifdef XDEBUG
                            G.mult(M.colPair(0,1).transpose());
#endif
                        }
                        if (i < N-1) {
                            z.swap(i,N-1);
                            D.swap(i,N-1);
                            if (U.cptr()) U.swapCols(i,N-1);
                            if (Vt.cptr()) Vt.swapRows(i,N-1);
#ifdef XDEBUG
                            M.swapCols(i,N-1);
                            M.swapRows(i,N-1);
#endif
                            i_firstswap = i;
                        } else {
                            tol = eps*D(N-2);
                        }
                        --N;
                    }
                }
#ifdef XDEBUG
                dbgcout<<"After second deflation\n";
                dbgcout<<"DN = "<<D.subVector(0,N)<<endl;
                dbgcout<<"zN = "<<z.subVector(0,N)<<endl;
                if (U.cptr() && Vt.cptr()) {
                    //dbgcout<<"M = "<<M<<endl;
                    //dbgcout<<"U = "<<U<<endl;
                    //dbgcout<<"Vt = "<<Vt<<endl;
                    //dbgcout<<"UMVt = "<<U * M * Vt<<endl;
                    dbgcout<<"Norm(UMVt-A0) = "<<Norm(U * M * Vt -A0)<<endl;
                    dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<endl;
                    dbgcout<<"Norm(VtV-1) = "<<Norm(Vt*Vt.adjoint()-T(1))<<endl;
                    if (Norm(U*M*Vt-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                    if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                    if (Norm(Vt*Vt.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
                }
#endif

                // Re-sort if we screwed up the sort in the second deflation.
                // The D(i) for i >= i_firstswap might be unsorted
                if (i_firstswap < N-1) {
                    D.subVector(i_firstswap,N).sort(P.get(),Ascend);
                    z.subVector(i_firstswap,N).permute(P.get());
                    if (U.cptr()) U.colRange(i_firstswap,N).permuteCols(P.get());
                    if (Vt.cptr()) Vt.rowRange(i_firstswap,N).permuteRows(P.get());
                }
#ifdef XDEBUG
                dbgcout<<"After second sort\n";
                dbgcout<<"DN = "<<D.subVector(0,N)<<endl;
                dbgcout<<"zN = "<<z.subVector(0,N)<<endl;
                if (U.cptr() && Vt.cptr()) {
                    if (i_firstswap < N-1) {
                        M.colRange(i_firstswap,N).permuteCols(P.get());
                        M.rowRange(i_firstswap,N).permuteRows(P.get());
                    }
                    //dbgcout<<"M = "<<M<<endl;
                    //dbgcout<<"UMVt = "<<U * M * Vt<<endl;
                    dbgcout<<"Norm(UMVt-A0) = "<<Norm(U * M * Vt -A0)<<endl;
                    dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<endl;
                    dbgcout<<"Norm(VtV-1) = "<<Norm(Vt*Vt.adjoint()-T(1))<<endl;
                    if (Norm(U*M*Vt-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                    if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                    if (Norm(Vt*Vt.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
                }
#endif
            }

            // After deflation, it is possible for N to be 1, in which case,
            // the SVD is done.
            if (N > 1) {
                VectorView<RT> DN = D.subVector(0,N);
                VectorView<RT> zN = z.subVector(0,N);

                // Find the new singular values. 
                Vector<RT> S(N);

                // First rescale DN,zN by their maximum value to 
                // avoid issues with underflow/overflow:
                RT scale = TMV_MAX(zN.maxAbs2Element(),TMV_ABS2(DN[N-1]));
                dbgcout<<"scale = "<<scale<<std::endl;
                DN /= scale;
                zN /= scale;
                dbgcout<<"After scale: DN = "<<DN<<std::endl;
                dbgcout<<"zN = "<<zN<<std::endl;

                if (U.cptr() || Vt.cptr()) {
                    Matrix<RT,ColMajor> W(N,N); // W(i,j) = D(i)-S(j)
                    FindDCSingularValues(S,RT(1),DN,zN,W);
                    dbgcout<<"S = "<<S<<endl;

                    // Update z's to numerically more stable values:
                    // z_i' = sqrt( (sn^2-di^2) 
                    //             prod_1..i-1 (sk^2-di^2)/(dk^2-di^2) *
                    //             prod_i..N-1 (sk^2-di^2)/(dk+1^2-di^2) ) 
                    //        * sign(z_i)
                    for(ptrdiff_t i=0;i<N;i++) {
                        RT di = D(i);
                        RT prod = -W(i,N-1)*(S(N-1) + di);
                        for(ptrdiff_t k=0;k<i;k++) {
                            prod *= -W(i,k)/(D(k)-di);
                            prod *= (S(k)+di)/(D(k)+di);
                        }
                        for(ptrdiff_t k=i+1;k<N;k++) {
                            prod *= -W(i,k-1)/(D(k)-di);
                            prod *= (S(k-1)+di)/(D(k)+di);
                        }
                        prod = TMV_SQRT(prod);
                        if (z(i) > 0) z(i) = prod;
                        else z(i) = -prod;
                    }

                    // Make X, Y matrices
                    // x_j = ( -1 , d2 z2/(d2^2-sj^2) , ... dN zN/(dN^2-sj^2) )
                    // X.col(j) = x_j / Norm(x_j);
                    // y_j = ( z1/(d1^2-sj^2) , z2/(d2^2-sj^2) , 
                    //         ... zN/(dN^2-sj^2) )
                    // Y.row(j) = y_j / Norm(y_j);

                    // Currently W(i,j) = di-sj
                    // First convert it to Y.transpose()
                    for(ptrdiff_t j=0;j<N;j++) {
                        VectorView<RT> yj = W.col(j);
                        Vector<RT> diff_j = yj; // copy current values
                        yj(0) = -z(0)/(S(j)*S(j));
                        for(ptrdiff_t i=1;i<N;i++) 
                            yj(i) = (z(i) / diff_j(i)) / (D(i)+S(j));
                    }
                    if (Vt.cptr()) {
                        Vector<RT> normyj(N);
                        for(ptrdiff_t j=0;j<N;j++) {
                            W.col(j) /= (normyj(j) = Norm(W.col(j)));
                        }
                        // Vt = Y * Vt
                        Vt.rowRange(0,N) = W.transpose() * Vt.rowRange(0,N);
                        if (U.cptr()) for(ptrdiff_t j=0;j<N;j++) W.col(j) *= normyj(j);
                    }
                    if (U.cptr()) {
                        // Now convert W from Y.transpose() to X
                        W.row(0).setAllTo(RT(-1));
                        for(ptrdiff_t i=1;i<N;i++) W.row(i) *= D(i);
                        for(ptrdiff_t j=0;j<N;j++) W.col(j) /= Norm(W.col(j));
                        // U = U * X
                        U.colRange(0,N) *= W;
                    }
                } else {
                    FindDCSingularValues(S,RT(1),DN,zN);
                    dbgcout<<"S = "<<S<<endl;
                }

                // Done.  Copy S to D and undo scaling.
                DN = scale * S;

                dbgcout<<"Done D = "<<D<<endl;
#ifdef XDEBUG
                if (U.cptr() && Vt.cptr()) {
                    M.subMatrix(0,N,0,N) = DiagMatrixViewOf(S);
                    //dbgcout<<"M = "<<M<<endl;
                    //dbgcout<<"UMVt = "<<U * M * Vt<<endl;
                    dbgcout<<"Norm(UMVt-A0) = "<<Norm(U * M * Vt -A0)<<endl;
                    dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<endl;
                    dbgcout<<"Norm(VtV-1) = "<<Norm(Vt*Vt.adjoint()-T(1))<<endl;
                    if (Norm(U*M*Vt-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                    if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                    if (Norm(Vt*Vt.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
                }
#endif
            } else if (N==1) {
                if (D(0) == RT(0)) {
                    D(0) = TMV_ABS(z(0));
                    if (Vt.cptr() && z(0) < RT(0)) Vt.row(0) *= RT(-1);
                }
#ifdef XDEBUG
                if (U.cptr() && Vt.cptr()) {
                    M(0,0) = D(0);
                    //dbgcout<<"M = "<<M<<endl;
                    //dbgcout<<"UMVt = "<<U * M * Vt<<endl;
                    dbgcout<<"Norm(UMVt-A0) = "<<Norm(U * M * Vt -A0)<<endl;
                    dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<endl;
                    dbgcout<<"Norm(VtV-1) = "<<Norm(Vt*Vt.adjoint()-T(1))<<endl;
                    if (Norm(U*M*Vt-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                    if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                    if (Norm(Vt*Vt.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
                }
#endif
            }
#ifdef XDEBUG
            if (U.cptr() && Vt.cptr()) {
                M = DiagMatrixViewOf(D);
                //dbgcout<<"M = "<<M<<endl;
                //dbgcout<<"UMVt = "<<U * M * Vt<<endl;
                dbgcout<<"Norm(UMVt-A0) = "<<Norm(U * M * Vt -A0)<<endl;
                dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<endl;
                dbgcout<<"Norm(VtV-1) = "<<Norm(Vt*Vt.adjoint()-T(1))<<endl;
                if (Norm(U*M*Vt-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                if (Norm(Vt*Vt.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
            }
#endif
        }

#ifdef XDEBUG
        DiagMatrix<RT> S(D);
        Matrix<T> A2(A0.colsize(),A0.rowsize());
        dbgcout<<"S = "<<S.diag()<<endl;
        if (U.cptr() && Vt.cptr()) {
            A2 = U * S * Vt;
            dbgcout<<"Norm(A2-A0) = "<<Norm(A2-A0)<<endl;
            dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<endl;
            dbgcout<<"Norm(VtV-1) = "<<Norm(Vt*Vt.adjoint()-T(1))<<endl;
            dbgcout<<"THRESH * normA0 = "<<THRESH<<" * "<<normA0<<" = "<<THRESH*normA0<<std::endl;
            if (!(Norm(A2-A0) < THRESH*normA0) 
#ifdef TESTUV
                || ( Vt.colsize()>=Vt.rowsize() && 
                     !(Norm(Vt.adjoint()*Vt-T(1)) < THRESH*normA0) )
                || ( Vt.rowsize()>=Vt.colsize() && 
                     !(Norm(Vt*Vt.adjoint()-T(1)) < THRESH*normA0) )
                || ( U.colsize()>=U.rowsize() && 
                     !(Norm(U.adjoint()*U-T(1)) < THRESH*normA0) ) 
#endif
            ) {
                cerr<<"SV_Decompose_DC: \n";
                cerr<<"input D = "<<D0<<endl;
                cerr<<"input E = "<<E0<<endl;
                cerr<<"output S = "<<D<<endl;
                Matrix<T> U_QR = U0;
                Matrix<T> Vt_QR = Vt0;
                Vector<RT> D_QR = D0;
                Vector<RT> E_QR = E0;
                SmallProblem<T>(
                    U_QR.view(),D_QR.view(),E_QR.view(),Vt_QR.view());
                cerr<<"QR solution: S = "<<D_QR<<endl;
                cerr<<"U = "<<U<<endl;
                cerr<<"QR solution: U = "<<U_QR<<endl;
                cerr<<"U-U_QR = "<<U-U_QR<<endl;
                cerr<<"Norm(U-U_QR) = "<<Norm(U-U_QR)<<endl; 
                cerr<<"Vt = "<<Vt<<endl;
                cerr<<"QR solution: Vt = "<<Vt_QR<<endl;
                cerr<<"Vt-Vt_QR = "<<Vt-Vt_QR<<endl;
                cerr<<"Norm(Vt-Vt_QR) = "<<Norm(Vt-Vt_QR)<<endl; 
                cerr<<"UBVt = "<<A0<<endl;
                cerr<<"USVt = "<<A2<<endl;
                abort();
            }
        }
#endif
    }

#undef RT

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SVDecompose_DC.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


