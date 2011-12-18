

#ifndef TMV_SVDecompose_DC_H
#define TMV_SVDecompose_DC_H

#include "TMV_SVDecompose.h"

#ifdef XDEBUG_SVD
#define THRESH 1.e-4
//#define TESTUV // Should only use this for full USV decomposition.
#include "TMV_SVDecompose_QR.h"
#include "TMV_TriMatrix.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef PRINTALGO_SVD
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#define TMV_USE_WRITER
#endif

#define dbgcout SafeWriter()
#include "TMV_SafeWriter.h"

#define TMV_DC_LIMIT 32
#define TMV_MAXITER 100

namespace tmv {

    // Defined in TMV_SVDecompose_DC.cpp
    template <class Tu, class RT, class Tv>
    void InstSV_DecomposeFromBidiagonal_DC(
        MatrixView<Tu> U, VectorView<RT> D, VectorView<RT> E,
        MatrixView<Tv> V, bool UisI, bool VisI);

    // First the function that finds the singular values for the join
    // step of the divide and conquer algorithm.
    template <class T, class V1, class V2, class V3, class V4, class V5>
    static T FindDCSingularValue(
        const int k, const int N, const T rho,
        const BaseVector_Calc<V1>& D, const BaseVector_Calc<V2>& z,
        const BaseVector_Calc<V3>& zsq, T normsqz,
        BaseVector_Mutable<V4>& diff, BaseVector_Mutable<V5>& sum)
    {
        dbgcout<<"start FindDCV: k = "<<k<<std::endl;
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

        dbgcout<<"k = "<<k<<std::endl;
        if (k<N-1) {
            dbgcout<<"root bracketed by "<<D[k]<<"  "<<D[k+1]<<std::endl;
            dbgcout<<"   range = "<<D[k+1]-D[k]<<std::endl;
            dbgcout<<"s^2 bracketed by "<<D[k]*D[k]<<"  "<<D[k+1]*D[k+1]<<std::endl;
            dbgcout<<"   range = "<<D[k+1]*D[k+1]-D[k]*D[k]<<std::endl;
        } else {
            dbgcout<<"root bracketed by "<<D[k]<<"  "<<
                TMV_SQRT(D[k]*D[k]+normsqz)<<std::endl;
            dbgcout<<"   range = "<<TMV_SQRT(D[k]*D[k]+normsqz)-D[k]*D[k]<<std::endl;
            dbgcout<<"s^2 bracketed by "<<D[k]*D[k]<<"  "<<
                D[k]*D[k]+normsqz<<std::endl;
            dbgcout<<"   range = "<<normsqz<<std::endl;
        }

        // Most of the routine is the same for k=N-1 as for the normal
        // case, but there are some lines where the normal case uses
        // k, but when k==N-1, I need k-1.  So I use kk for these.
        int kk = k < N-1 ? k : k-1;

        // First we need an initial guess for s_k 
        dbgcout<<"D[kk+1]-D[kk] = "<<D[kk+1]-D[kk]<<std::endl;
        dbgcout<<"eps*D["<<kk+1<<"] = "<<eps*D[kk+1]<<std::endl;
        T delta = (D[kk+1] - D[kk]) * (D[kk+1] + D[kk]);
        // It is known that D_k < s_k < D_k+1
        // Let s^2 = (D_k^2 + D_k+1^2)/2
        // Define t,tau: s = D_k + t, and s^2 = D_k^2 + tau
        T tau = k < N-1 ? delta/T(2) : normsqz/(T(2)*rho);
        T t = tau / (D[k] + TMV_SQRT(D[k]*D[k]+tau));
        T s = D[k] + t;
        dbgcout<<"Initial s = "<<s<<" = "<<D[k]<<" + "<<t<<std::endl;

        for(int j=0;j<N;j++) sum[j] = D[j]+s;
        for(int j=0;j<N;j++) diff[j] = D[j]-D[k];
        for(int j=0;j<N;j++) diff[j] -= t;

        // Let c = f(s) - z_k^2/(D_k^2-s^2) - z_k+1^2/(D_k+1^2-s^2)
        // i.e. c is f(s) without the two terms for the poles that
        // bracket the root we are seraching for.
        // Since the largest terms are normally the ones for the poles nearer
        // to s, we perform the calculation from each end and proceed towards k.
        // Note: for k==N-1, the two poles do not bracket the root, but 
        // we proceed similarly anyway.
        T psi(0);
        for(int j=0;j<kk;j++) {
            psi += (zsq[j] / sum[j]) / diff[j];
            //dbgcout<<"psi += "<<zsq[j]<<" / ("<<diff[j]<<" * "<<sum[k]<<") = "<<psi<<std::endl;
        }
        T phi(0);
        for(int j=N-1;j>kk+1;j--) {
            phi += (zsq[j] / sum[j]) / diff[j];
            //dbgcout<<"phi += "<<zsq[j]<<" / ("<<diff[j]<<" * "<<sum[k]<<") = "<<phi<<std::endl;
        }
        T c = rho + psi + phi;
        dbgcout<<"c = "<<c<<std::endl;

        // Finish calculating f = f(s)
        T f = c + (zsq[kk] / sum[kk]) / diff[kk] +
            (zsq[kk+1] / sum[kk+1]) / diff[kk+1];
        dbgcout<<"f = c + "<<zsq[kk]<<" / ("<<diff[kk]<<" * "<<sum[kk]<<")\n";
        dbgcout<<"      + "<<zsq[kk+1]<<" / ("<<
            diff[kk+1]<<" * "<<sum[kk+1]<<")\n";
        dbgcout<<"f("<<tau<<") = "<<f<<std::endl;

        T lowerbound, upperbound; // the allowed range for solution tau
        int k1; // The index of the pole closest to the root
        if (k < N-1) {
            dbgcout<<"k < N-1\n";
            if (f >= T(0)) { // D_k < s_k < s
                dbgcout<<"f >= 0\n";
                k1 = k;
                lowerbound = T(0);
                upperbound = tau;
                dbgcout<<"bounds = 0,tau = "<<tau<<std::endl;
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
                dbgcout<<"bounds = -tau,0 = "<<-tau<<std::endl;
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
                    dbgcout<<"bounds = 0,tau = "<<tau<<std::endl;
                } else { // s < s_k < sqrt(D_k^2 + |z|^2)
                    dbgcout<<"f < 0\n";
                    lowerbound = tau;
                    upperbound = normsqz/rho;
                    dbgcout<<"bounds = tau,2tau = "<<tau<<"  "<<
                        normsqz/rho<<std::endl;
                    if (c <= T(0)) tau = upperbound; break;
                    // In this case there is no positive solution.  Best guess
                    // is just to start at upper bound.
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
                // But this time, if f <  0, the solution to the quadratic
                // can really be much > upperbound.  So the check is especially
                // important in this case.
            } while (false);
        }
        dbgcout<<"Initial tau = "<<tau<<std::endl;
        dbgcout<<"Allowed bounds for tau = "<<
            lowerbound<<"  "<<upperbound<<std::endl;
        TMVAssert(tau >= lowerbound && tau <= upperbound);
        t = tau / ( TMV_SQRT(D[k1]*D[k1] + tau) + D[k1] );

        s = D[k1] + t;
        dbgcout<<"First refinement: s = "<<s<<" = "<<D[k1]<<" + "<<t<<std::endl;
        // We should be able to write diff[j] = (D[j]-D[k1)-t,
        // but some compilers do optimizations that lose accuracy 
        // and can end up with diff = 0.
        // So need to do this in two steps.
        for(int j=0;j<N;j++) sum[j] = D[j]+s;
        for(int j=0;j<N;j++) diff[j] = D[j]-D[k1];
        for(int j=0;j<N;j++) diff[j] -= t;

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
        int iter = 0;
        do {
            // Calculate psi, phi
            psi = phi = e = dpsi = dphi = T(0);
            for(int j=0;j<k1;j++) {
                T temp = (z[j] / sum[j]) / diff[j];
                psix = psi; dpsix = dpsi;
                psi += z[j] * temp;
                dpsi += temp * temp;
                e -= psi;
            }
            for(int j=N-1;j>k1;j--) {
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
            dbgcout<<"rho = "<<rho<<std::endl;
            dbgcout<<"psi = "<<psi<<std::endl;
            dbgcout<<"phi = "<<phi<<std::endl;
            dbgcout<<"fk = "<<fk<<std::endl;
            dbgcout<<"k1 = "<<k1<<"  k = "<<k<<std::endl;
            dbgcout<<"sum = "<<sum[k1]<<std::endl;
            dbgcout<<"diff = "<<diff[k1]<<std::endl;
            dbgcout<<"z = "<<z[k1]<<std::endl;
            dbgcout<<"temp = "<<temp<<std::endl;
            dbgcout<<"k1term = "<<k1term<<std::endl;
            dbgcout<<"f = "<<f<<std::endl;
            dbgcout<<"df = "<<df<<std::endl;

            if (k == N-1) break;
            //if (TMV_ABS(f) <= eps*e) break; 
            // Do it this way instead so nan will break out too.
            if (!(TMV_ABS(f) > eps*e)) break;
            if (TMV_ABS(f) <= T(1.e3)*TMV_ABS(f1) && f*f1 < T(0)) break;

            // If we get to here, then the initial guess was not very good.  
            // This can mean that the function is not very well described
            // by the nearby poles.
            // So use bisection until we at least do better than the 
            // midpoint guess we started with.
            // Then we can move on to the fancy, quick-converging techniques 
            // below.
            // This is quite rare, but when the above calculation gives a very 
            // bad initial guess, this adjustment can be critical.
            dbgcout<<"f = "<<f<<", f1 = "<<f1<<std::endl;
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
            dbgcout<<"dt = "<<dt<<std::endl;
            tau += dt;
            dbgcout<<"tau -> "<<tau<<std::endl;
            dbgcout<<"bounds -> "<<lowerbound<<"  "<<upperbound<<std::endl;

            // change dt from d(tau) to d(t)
            dt /= s + TMV_SQRT(s*s+dt);
            t += dt;
            s += dt;
            for(int j=0;j<N;j++) sum[j] = D[j]+s;
            if (TMV_ABS(s) < TMV_ABS(dt)) {
                dbgcout<<"Need to redo the diff calculations\n";
                for(int j=0;j<N;j++) diff[j] = D[j]-s;
            } else  {
                for(int j=0;j<N;j++) diff[j] -= dt;
            }

            if (iter == TMV_MAXITER-1) {
                TMV_Warning(
                    "Warning - Unable to find appropriate initial "
                    "guess in FindDCSingularValue\n"
                    "Unsuccessful after 100 iterations.\n");
#ifdef TMVDEBUG
                exit(1);
#endif
            }
        } while (++iter < TMV_MAXITER);

        // Use 3 poles when fk<0 for k1=k or fk>0 for k1=k+1
        bool threepoles = (((fk < T(0)) == (k1 == k)) && k1>0 && k1<N-1);
        bool fixedweight = true;
        if (threepoles) dbgcout<<"USE THREE POLES\n";
        bool last = false;
        iter = 0;
        do {
            dbgcout<<"k = "<<k<<", loop 2: iter = "<<iter<<std::endl;
            dbgcout<<"Main iter = "<<iter<<", f = "<<f<<
                ", eps*e = "<<eps*e<<std::endl;
            dbgcout<<"f("<<tau<<") = "<<f<<"  s = "<<s<<
                ", s^2 = "<<D[k1]*D[k1]+tau<<std::endl;
            dbgcout<<"df = "<<df<<", eta_newt = -f/df = "<<-f/df<<std::endl;

            T eta; // = the change to be added to s

            if (!(TMV_ABS(f) > eps*e)) {
                last = true;
                eta = T(0); // This will trigger a Newton step below.
            } else if (!threepoles) {
                // The normal case: only 2 poles used for approx
                T d1 = diff[kk] * sum[kk];
                T d2 = diff[kk+1] * sum[kk+1];
                dbgcout<<"d1,d2 = "<<d1<<", "<<d2<<std::endl;
                if (fixedweight)
                    if (kk==k1) c = fk - d2*dfk;
                    else c = fk - d1*dfk;
                else c = fk - d1*dpsi - d2*dphi;
                // c eta^2 - a eta + b = 0
                T d1x = d1/d2;
                T a = (d1x+T(1))*f - d1*df;
                T b = d1x*f;
                dbgcout<<"c,a,b = "<<c<<", "<<a<<", "<<b<<std::endl;
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
                    dbgcout<<"c=0  eta = "<<eta<<std::endl;
                    if (k==N-1)
                        TMVAssert(eta > -tau);
                    else
                        TMVAssert(eta > d1 && eta < d2);
                } else {
                    T d = a*a-T(4)*b*c;
                    dbgcout<<"d = "<<d<<std::endl;
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
                        dbgcout<<"eta = "<<eta<<std::endl;
                        if (eta <= -tau) // Then rounding errors - use bisect
                            eta = -tau/T(2);
                        TMVAssert(eta > -tau);
                    } else {
                        eta = a <= 0 ? (a-d)/(T(2)*c) : T(2)*b/(a+d);
                        eta *= d2;
                        dbgcout<<"eta = "<<eta<<std::endl;
                        if (eta <= d1) // Then rounding errors - use bisect
                            eta = d1/T(2);
                        if (eta >= d2) // Then rounding errors - use bisect
                            eta = d2/T(2);
                        dbgcout<<"eta => "<<eta<<std::endl;
                        TMVAssert(eta > d1 && eta < d2);
                    }
                }
            } else {
                // More complicated case using 3 poles. 
                dbgcout<<"threepoles section\n";
                T d1 = diff[k1-1] * sum[k1-1];
                T d2 = diff[k1+1] * sum[k1+1];
                dbgcout<<"d1,d2 = "<<d1<<" "<<d2<<std::endl;
                if (fixedweight)
                    if (k1 == k) c = (rho+psix+phi) - d2*(dpsix+dphi);
                    else  c = (rho+psi+phix) - d1*(dpsi+dphix);
                else
                    if (k1 == k) c = (rho+psix+phi) - d1*dpsix - d2*dphi;
                    else c = (rho+psi+phix) - d1*dpsi - d2*dphix;
                dbgcout<<"c = "<<c<<std::endl;
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
                dbgcout<<"temp = "<<temp<<std::endl;
                dbgcout<<"a = "<<a<<std::endl;
                dbgcout<<"b = "<<b<<std::endl;
                dbgcout<<"g = "<<g<<std::endl;
                dbgcout<<"mineta, maxeta = "<<mineta<<", "<<maxeta<<std::endl;

                // Bounds on eta are:
                // lb < tau + eta < ub 
                // lb - tau < eta < ub - tau
                // h(0) = g > 0 regardless of the sign of f.
                // if f > 0, then deta needs to be negative
                // if f < 0, then deta needs to be positive
                // This dictates which bound is the one to use for eta2
                // In below iteration, we maintain h(eta2) * h(eta) < 0
                T eta2 = (f>0) ? mineta : maxeta;
                dbgcout<<"eta2 = "<<eta2<<std::endl;
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
                    dbgcout<<"h(0) = g = "<<h<<std::endl;
                    dbgcout<<"h(eta2) = h2 = "<<h2<<std::endl;
                    dbgcout<<"Same sign, so use Newton\n";
                    eta = (-f/df)/d2;
                } else for(int iter3 = 0; iter3 < TMV_MAXITER; iter3++) {
                    dbgcout<<"iter3 = "<<iter3<<", eta = "<<eta<<
                        "  h = "<<h<<std::endl;
                    dbgcout<<"eta2 = "<<eta2<<"  h2 = "<<h2<<std::endl;
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
                    dbgcout<<"dh,d2h,x = "<<dh<<", "<<d2h<<", "<<x<<std::endl;
                    x = (x <= T(0)) ? T(0) : TMV_SQRT(x);
                    T deta = -T(2)*h / (dh*(T(1) + x));
                    dbgcout<<"deta = "<<deta<<std::endl;
                    if ((eta2 > eta && (deta < 0 || deta > eta2-eta)) ||
                        (eta2 < eta && (deta > 0 || deta < eta2-eta)) ) {
                        dbgcout<<"deta is wrong direction or too big - "
                            "use other root.\n";
                        deta = -dh*(T(1)+x) / d2h;
                        dbgcout<<"deta = "<<deta<<std::endl;
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
                            dbgcout<<"deta -> "<<deta<<std::endl;
#endif
                        }
                    }
                    TMVAssert(eta + deta <= maxeta);
                    TMVAssert(eta + deta >= mineta);
                    T etanew = eta + deta;
                    T hnew =
                        (TMV_ABS(deta) < T(1.e-3)*TMV_ABS(eta)) ?
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
                    dbgcout<<"etanew = "<<etanew<<std::endl;
                    dbgcout<<"hnew = "<<hnew<<std::endl;
                    if ( (h > T(0)) != (hnew > T(0)) ) { eta2 = eta; h2 = h; }
                    eta = etanew;
                    h = hnew;
                    if (TMV_ABS(deta) < eps*TMV_ABS(eta) ||
                        TMV_ABS(h) < eps*g) break;
                    if (iter3 == TMV_MAXITER-1) {
                        TMV_Warning(
                            "Warning - unable to converge for THREEPOLES "
                            "solution\n"
                            "No solution after 100 iterations.\n");
#ifdef TMVDEBUG
                        exit(1);
#endif
                    }
                }
                eta *= d2;
                dbgcout<<"rescaled eta => "<<eta<<std::endl;
            }
            dbgcout<<"eta = "<<eta<<std::endl;

            // Check for eta being the wrong sign.
            // If f < 0, then eta should be > 0.
            // If f > 0, then eta should be < 0.
            // If not, use a Newton step instead for which this is guaranteed.
            if (f*eta >= T(0)) {
                eta = -f/df;
                dbgcout<<"Newton eta = "<<eta<<std::endl;
            }

            // Also check that tau+eta is still within allowed bounds:
            if (tau+eta < lowerbound) {
                eta = (lowerbound-tau)/T(2);
                dbgcout<<"halfway to lowerbound eta = "<<eta<<std::endl;
            } else if (tau+eta > upperbound) {
                eta = (upperbound-tau)/T(2);
                dbgcout<<"halfway to upperbound eta = "<<eta<<std::endl;
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
            for(int j=0;j<N;j++) sum[j] = D[j]+s;
            if (TMV_ABS(s) < TMV_ABS(eta)) {
                // Then the iterative adjustment to diff isn't going to 
                // maintain the accuracy we need.
                dbgcout<<"Need to redo the diff calculations\n";
                for(int j=0;j<N;j++) diff[j] = D[j]-s;
            } else {
                for(int j=0;j<N;j++) diff[j] -= eta;
            }
            if (last) break;

            // Update psi, phi, etc.
            psi = T(0);
            dpsi = T(0);
            e = T(0);
            for(int j=0;j<k1;j++) {
                T temp = (z[j] / sum[j]) / diff[j];
                psix = psi; dpsix = dpsi;
                psi += z[j] * temp;
                dpsi += temp * temp;
                e -= psi;
            }
            phi = T(0);
            dphi = T(0);
            for(int j=N-1;j>k1;j--) {
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
                    newthreepoles<<std::endl;
                dbgcout<<"fk = "<<fk<<"  fnew = "<<fnew<<std::endl;
                threepoles = newthreepoles;
            }

            // Update bounds:
            if (fnew < T(0)) lowerbound = tau;
            else upperbound = tau;
            dbgcout<<"New bounds = "<<lowerbound<<"  "<<upperbound<<std::endl;

            // Update scheme:
            if (fnew*f > T(0) && TMV_ABS(fnew) > T(0.1) * TMV_ABS(f) && k<N-1)
                fixedweight = !fixedweight;

            f = fnew;
            if (iter == TMV_MAXITER-1) {
                TMV_Warning(
                    "Warning - Unable to find solution in "
                    "FindDCSingularValue\n"
                    "No solution after 100 iterations.\n");
#ifdef TMVDEBUG
                exit(1);
#endif
            }
        } while (++iter < TMV_MAXITER);
        dbgcout<<"Found Singularvalue S("<<k<<") = "<<s<<std::endl;
        dbgcout<<"eta = "<<-f/df<<", tau = "<<tau<<
            ", eta/tau = "<<-f/df/tau<<std::endl;
#ifdef XDEBUG_SVD
        // f(s) = rho + Sum_i=1..N z_i^2/(D_i^2-s^2)
        T ff = rho;
        dbgcout<<"f(s) = "<<rho<<" + ";
        for(int j=0;j<N;j++) {
            dbgcout<<(z[j]/diff[j])*(z[j]/sum[j])<<" + ";
            ff += (z[j]/diff[j])*(z[j]/sum[j]);
        }
        dbgcout<<" = "<<ff<<std::endl;
        dbgcout<<"cf f = "<<f<<std::endl;
#endif
        return s;
    }

    template <class V1, class T, class V2, class V3, class M1>
    void FindDCSingularValues(
        BaseVector_Mutable<V1>& S, T rho,
        const BaseVector_Calc<V2>& D, const BaseVector_Calc<V3>& z,
        BaseMatrix_Rec_Mutable<M1>& diffmat)
    {
        dbgcout<<"Start FindDCSV: "<<std::endl;
        dbgcout<<"D = "<<D<<std::endl;
        dbgcout<<"z = "<<z<<std::endl;
        dbgcout<<"rho = "<<rho<<std::endl;
        TMVAssert(rho > T(0));
        TMVAssert(S.size() == D.size());
        TMVAssert(S.size() == z.size());
        TMVAssert(S.size() == diffmat.colsize());
        TMVAssert(S.size() == diffmat.rowsize());

        const int N = S.size();
        Vector<T> zsq(N);
        for(int j=0;j<N;j++) zsq[j] = z[j]*z[j];
        const T normsqz = zsq.sumElements();

#ifdef _OPENMP
#pragma omp parallel 
        {
            Vector<T> diff(N);
            Vector<T> sum(N);
            T Sk;
#pragma omp for
            for(int k=0;k<N;k++) {
                Sk = FindDCSingularValue(k,N,rho,D,z,zsq,normsqz,diff,sum);
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    dbgcout<<"Before assign answer for  k = "<<k<<std::endl;
                    S[k] = Sk;
                    diffmat.col(k) = diff;
                    dbgcout<<"Done assign answer for  k = "<<k<<std::endl;
                }
            }
            dbgcout<<"After omp for loop\n";
        }
        dbgcout<<"After omp parallel region\n";
#else
        Vector<T> sum(N);
        for(int k=0;k<N;k++) {
            typename M1::col_type diff = diffmat.col(k);
            S[k] = FindDCSingularValue(k,N,rho,D,z,zsq,normsqz,diff,sum);
        }
#endif
        dbgcout<<"S => "<<S<<std::endl;
    }

    template <class V1, class T, class V2, class V3>
    void FindDCSingularValues(
        BaseVector_Mutable<V1>& S, T rho,
        const BaseVector_Calc<V2>& D, const BaseVector_Calc<V3>& z)
    {
        dbgcout<<"Start FindDCSV (No diffmat): "<<std::endl;
        dbgcout<<"D = "<<D<<std::endl;
        dbgcout<<"z = "<<z<<std::endl;
        TMVAssert(rho > T(0));
        TMVAssert(S.size() == D.size());
        TMVAssert(S.size() == z.size());

        const int N = S.size();
        Vector<T> zsq(N);
        for(int j=0;j<N;j++) zsq[j] = z[j]*z[j];
        T normsqz = zsq.sumElements();

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            Vector<T> diff(N);
            Vector<T> sum(N);
            T Sk;
#ifdef _OPENMP
#pragma omp for
#endif
            for(int k=0;k<N;k++) {
                Sk = FindDCSingularValue(k,N,rho,D,z,zsq,normsqz,diff,sum);
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    S[k] = Sk;
                }
            }
        }
        dbgcout<<"S => "<<S<<std::endl;
    }

    template <int algo, int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_DC_Helper;

    // algo 0: Trivial, nothing to do (M == 0, or N == 0)
    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_DC_Helper<0,cs,rs,Mu,Vd,Ve,Mv>
    { static TMV_INLINE void call(Mu& , Vd& , Ve& , Mv& , bool , bool ) {} };

    // algo 1: Small problem: use the QR algorithm
    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_DC_Helper<1,cs,rs,Mu,Vd,Ve,Mv>
    {
        static void call(Mu& U, Vd& D, Ve& E, Mv& V, bool UisI, bool VisI)
        {
            const int N = rs==TMV_UNKNOWN ? int(D.size()) : rs;
#ifdef PRINTALGO_SVD
            const int M = cs==TMV_UNKNOWN ? int(U.colsize()) : cs;
            std::cout<<"SVDecomposeFromBidiagonal algo 1: M,N,cs,rs = "<<
                M<<','<<N<<','<<cs<<','<<rs<<std::endl;
#endif
            dbgcout<<"Start SmallProblem: N = "<<D.size()<<std::endl;
            SV_DecomposeFromBidiagonal_QR(U,D,E,V,UisI,VisI);
            dbgcout<<"After QR"<<std::endl;
            // Make all of the singular values positive
            typename Vd::iterator Di = D.begin();
            for(int i=0;i<N;++i,++Di) if (*Di < 0) {
                *Di = -(*Di);
                if (V.cptr()) V.row(i) = -V.row(i);
            }
            dbgcout<<"After make all Di positive"<<std::endl;
        }
    };

    // algo 11: Normal divide and conquer algorithm
    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_DC_Helper<11,cs,rs,Mu,Vd,Ve,Mv>
    {
        static void call(Mu& U, Vd& D, Ve& E, Mv& V, bool UisI, bool VisI)
        {
            typedef typename Vd::value_type RT;
            TMVStaticAssert(Traits<RT>::isreal);

            int N = rs==TMV_UNKNOWN ? int(D.size()) : rs;
#ifdef PRINTALGO_SVD
            const int M1 = cs==TMV_UNKNOWN ? int(U.colsize()) : cs;
            std::cout<<"SVDecomposeFromBidiagonal algo 11: M,N,cs,rs = "<<
                M1<<','<<N<<','<<cs<<','<<rs<<std::endl;
#endif
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
            // B1 = U1 (S1 0) V1
            // B2 = U2 S2 V2
            //
            //
            //     [ U1 0  0  ] [ (S1 0)  0  ] [ V1  0 ]
            // B = [ 0  1  0  ] [   z1    z2 ] [ 0  V2 ]
            //     [ 0  0  U2 ] [   0     S2 ] 
            //
            // The vector z = [ z1 z2 ] is such that z V = [ 0 D_K E_K 0 ]
            // which implies that:
            // z1 V1 = D_K e_K
            // z1 = D_K e_K V1T
            //    = D_K V1T.row(K)
            //    = D_K V1.col(K)
            //
            // z2 V2 = E_K e_1
            // z2 = E_K e_1 V2T
            //    = E_K V2T.row(1)
            //    = E_K V2.col(1)
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
            if (V.cptr()) TMVAssert(V.colsize() == D.size());
            TMVAssert(D.step()==1);
            TMVAssert(E.step()==1);

            dbgcout<<"Start SVD DC: \n";
            dbgcout<<"U: "<<U.cptr()<<"  "<<U.colsize()<<"  "<<U.rowsize()<<"  "<<U.stepi()<<"  "<<U.stepj()<<std::endl;
            dbgcout<<"V: "<<V.cptr()<<"  "<<V.colsize()<<"  "<<V.rowsize()<<"  "<<V.stepi()<<"  "<<V.stepj()<<std::endl;
#ifdef XDEBUG_SVD
            typedef typename Traits2<typename Mu::value_type, typename Mv::value_type>::type T;
            Vector<RT> D0 = D;
            dbgcout<<"D0 = "<<D0<<std::endl;
            Vector<RT> E0 = E;
            dbgcout<<"E0 = "<<E0<<std::endl;
            Matrix<T> U0(U.cptr()?U.colsize():D.size(),D.size());
            if (U.cptr()) U0 = U; else U0.setToIdentity();
            Matrix<T> V0(D.size(),V.cptr()?V.rowsize():D.size());
            if (V.cptr()) V0 = V; else V0.setToIdentity();
            Matrix<T> A0(U0.colsize(),V0.rowsize());
            Matrix<T> B0(D.size(),D.size(),RT(0));
            B0.diag() = D;
            if (E0.size() > 0) B0.diag(1) = E;
            A0 = U0 * B0 * V0;
            //dbgcout<<"A0 = "<<A0<<std::endl;
            double normA0 = Norm(A0);
            if (normA0 < 1.) normA0 = 1.;
            dbgcout<<"Norm(A0) = "<<Norm(A0)<<std::endl;
            if (U.cptr()) dbgcout<<"Norm(UtU -1) = "<<Norm(U0.adjoint()*U0-T(1))<<std::endl;
            else dbgcout<<"No input U\n";
            if (V.cptr()) dbgcout<<"Norm(VVt -1) = "<<Norm(V0*V0.adjoint()-T(1))<<std::endl;
            else dbgcout<<"No input V\n";
#endif
            typedef typename Vd::subvector_type Vds;
            typedef typename Ve::subvector_type Ves;
            typedef Matrix<RT,RowMajor|NoDivider> Mvc;
            typedef typename Mvc::view_type Mvcv;
            typedef typename Mvc::rowrange_type Mvcs;
            typedef Matrix<RT,ColMajor|NoDivider> Muc;
            typedef typename Muc::view_type Mucv;
            const int xx = TMV_UNKNOWN;

            // If N is too small, use the QR method

            // Also, when we aren't making the U and V matrices, it seems 
            // that the QR algorithm is faster than DC for N < ~ 5000
            // And when N is much larger than this, then we start to worry about
            // memory issues of making the U and V matrices for the subproblems,
            // and QR is certainly never very much slower than doing 
            // 1 or 2 divides, which seems to be optimal for large problems.
            // So I always use QR for the S-only calculation.
            if (N <= TMV_DC_LIMIT || (!U.cptr() && !V.cptr())) {
                dbgcout<<"Do small problem:\n";
                return SVDecomposeFromBidiagonal_DC_Helper<
                    1,cs,rs,Mu,Vd,Ve,Mv>::call(U,D,E,V,UisI,VisI);
            }
            dbgcout<<"N > "<<TMV_DC_LIMIT<<std::endl;
            int K = (N-1)/2;
            const RT DK = D(K);
            const RT EK = E(K);
            Vector<RT> z(N,RT(0));

            // Do the left sub-problem
            Vds D1 = D.subVector(0,K);
            Ves E1 = E.subVector(0,K);
            Ves E1x = E.subVector(0,K-1);
            Mvc V1(K+1,V.cptr()?K+1:1);
            dbgcout<<"V1: "<<V1.cptr()<<"  "<<V1.colsize()<<"  "<<V1.rowsize()<<"  "<<V1.stepi()<<"  "<<V1.stepj()<<std::endl;
            Mvcv V1v = V1.view();
            Mvcs V1x = V1.rowRange(0,K);
            if (V.cptr()) V1.setToIdentity();
            else V1.col(0).makeBasis(K); // only need col(K)
            BidiagonalZeroLastCol(D1,E1,V1v);
            if (U.cptr()) {
                Muc U1(K,K);
                dbgcout<<"U1: "<<U1.cptr()<<"  "<<U1.colsize()<<"  "<<U1.rowsize()<<"  "<<U1.stepi()<<"  "<<U1.stepj()<<std::endl;
                Mucv U1v = U1.view();
                U1.setToIdentity();
                SVDecomposeFromBidiagonal_DC_Helper<
                    11,xx,xx,Mucv,Vds,Ves,Mvcv>::call(
                        U1v,D1,E1x,V1x,true,false);
                if (UisI) U.subMatrix(0,K,0,K) = U1;
                else U.colRange(0,K) *= U1;
            } else {
                SVDecomposeFromBidiagonal_DC_Helper<
                    11,xx,xx,Mu,Vds,Ves,Mvcv>::call(
                        U,D1,E1x,V1x,false,false);
            }
            dbgcout<<"After left sub-problem\n";
            z.subVector(0,K+1) = DK * V1.col(V.cptr() ? K : 0);
            dbgcout<<"After z = DK * V1.col\n";
            if (V.cptr()) {
                dbgcout<<"VisI = "<<VisI<<std::endl;
                dbgcout<<"K = "<<K<<std::endl;
                dbgcout<<"V steps = "<<V.stepi()<<"  "<<V.stepj()<<std::endl;
                dbgcout<<"V.rowRange(0,K+1) = "<<V.rowRange(0,K+1)<<std::endl;
                dbgcout<<"steps = "<<V.rowRange(0,K+1).stepi()<<"  "<<V.rowRange(0,K+1).stepj()<<std::endl;
                dbgcout<<"V const steps = "<<V.constView().stepi()<<"  "<<V.constView().stepj()<<std::endl;
                dbgcout<<"V.const.rowRange(0,K+1) = "<<V.constView().rowRange(0,K+1)<<std::endl;
                dbgcout<<"steps = "<<V.constView().rowRange(0,K+1).stepi()<<"  "<<V.constView().rowRange(0,K+1).stepj()<<std::endl;
                if (VisI) V.subMatrix(0,K+1,0,K+1) = V1;
                else V.rowRange(0,K+1) = V1 * V.rowRange(0,K+1);
                dbgcout<<"After set V\n";
            }

            // Do the right sub-problem
            Vds D2 = D.subVector(K+1,N);
            Ves E2 = E.subVector(K+1,N-1);
            Mvc V2(N-K-1,V.cptr()?N-K-1:1);
            dbgcout<<"V2: "<<V2.cptr()<<"  "<<V2.colsize()<<"  "<<V2.rowsize()<<"  "<<V2.stepi()<<"  "<<V2.stepj()<<std::endl;
            Mvcv V2v = V2.view();
            if (V.cptr()) V2.setToIdentity();
            else V2.col(0).makeBasis(0); // only need col(0)
            if (U.cptr()) {
                Muc U2(N-K-1,N-K-1);
                dbgcout<<"U2: "<<U2.cptr()<<"  "<<U2.colsize()<<"  "<<U2.rowsize()<<"  "<<U2.stepi()<<"  "<<U2.stepj()<<std::endl;
                Mucv U2v = U2.view();
                U2.setToIdentity();
                SVDecomposeFromBidiagonal_DC_Helper<
                    11,xx,xx,Mucv,Vds,Ves,Mvcv>::call(
                        U2v,D2,E2,V2v,true,V.cptr());
                if (UisI) U.subMatrix(K+1,N,K+1,N) = U2;
                else U.colRange(K+1,N) *= U2;
            } else {
                SVDecomposeFromBidiagonal_DC_Helper<
                    11,xx,xx,Mu,Vds,Ves,Mvcv>::call(
                        U,D2,E2,V2v,false,V.cptr());
            }
            dbgcout<<"After right sub-problem\n";
            z.subVector(K+1,N) = EK * V2.col(0);
            D(K) = RT(0);
            if (V.cptr()) {
                if (VisI) V.subMatrix(K+1,N,K+1,N) = V2;
                else V.rowRange(K+1,N) = V2 * V.rowRange(K+1,N);
            }
            dbgcout<<"Done subproblems (N="<<N<<")\n";
#ifdef XDEBUG_SVD
            dbgcout<<"D = "<<D<<std::endl;
            dbgcout<<"z = "<<z<<std::endl;
            Matrix<T> M(N,N,T(0));
            M.diag() = D;
            M.row(K) = z;
            if (U.cptr() && V.cptr()) {
                //dbgcout<<"M = "<<M<<std::endl;
                //dbgcout<<"UMV = "<<U*M*V<<std::endl;
                dbgcout<<"Norm(UMV-A0) = "<<Norm(U*M*V -A0)<<std::endl;
                dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
                dbgcout<<"Norm(VVt-1) = "<<Norm(V*V.adjoint()-T(1))<<std::endl;
                if (Norm(U*M*V-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                if (Norm(V*V.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
            }
#endif

            // Check for deflation step 1:
            const RT eps = TMV_Epsilon<RT>();
            RT norminf = D.subVector(0,N).normInf();
            RT tol = eps*norminf;
            dbgcout<<"tol = "<<tol<<std::endl;
            dbgcout<<"Initially D("<<K<<") = "<<D(K)<<std::endl;
            int pivot = K;
            for(int i=0;i<N-1;i++) {
                dbgcout<<"D("<<i<<") = "<<D(i)<<std::endl;
                dbgcout<<"z("<<i<<") = "<<z(i)<<std::endl;
                dbgcout<<"pivot = "<<pivot<<"  "<<D(pivot)<<std::endl;
                if (i != pivot && // Make sure not to deflate the pivot value.
                    (TMV_ABS(z(i)) < tol || z(i)*z(i)*eps == RT(0))) {
                    dbgcout<<"Deflate this one.\n";
                    z.swap(i,N-1);
                    D.swap(i,N-1);
                    if (U.cptr()) U.swapCols(i,N-1);
                    if (V.cptr()) V.swapRows(i,N-1);
                    if (pivot == N-1) pivot = i;
#ifdef XDEBUG_SVD
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
            dbgcout<<"After deflation\n";
#ifdef XDEBUG_SVD
            dbgcout<<"DN = "<<D.subVector(0,N)<<std::endl;
            dbgcout<<"zN = "<<z.subVector(0,N)<<std::endl;
            if (U.cptr() && V.cptr()) {
                //dbgcout<<"M = "<<M<<std::endl;
                //dbgcout<<"UMV = "<<U*M*V<<std::endl;
                dbgcout<<"Norm(UMV-A0) = "<<Norm(U*M*V -A0)<<std::endl;
                dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
                dbgcout<<"Norm(VVt-1) = "<<Norm(V*V.adjoint()-T(1))<<std::endl;
                if (Norm(U*M*V-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                if (Norm(V*V.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
            }
#endif

            if (N > 1) {
                // Sort the inner matrix so z is in the first row, D are increasing
                Permutation P(N);
                D.subVector(0,N).sort(P,Ascend);
                z.subVector(0,N) = P * z.subVector(0,N);
                if (U.cptr()) U.colRange(0,N) = U.colRange(0,N) * P.transpose();
                if (V.cptr()) V.rowRange(0,N) = P * V.rowRange(0,N);
                dbgcout<<"After sort\n";
#ifdef XDEBUG_SVD
                dbgcout<<"DN = "<<D.subVector(0,N)<<std::endl;
                dbgcout<<"zN = "<<z.subVector(0,N)<<std::endl;
                if (U.cptr() && V.cptr()) {
                    M.colRange(0,N) = M.colRange(0,N) * P.transpose();
                    M.rowRange(0,N) = P * M.rowRange(0,N);
                    //dbgcout<<"M = "<<M<<std::endl;
                    //dbgcout<<"U = "<<U<<std::endl;
                    //dbgcout<<"V = "<<V<<std::endl;
                    //dbgcout<<"UMV = "<<U*M*V<<std::endl;
                    dbgcout<<"Norm(UMV-A0) = "<<Norm(U*M*V -A0)<<std::endl;
                    dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
                    dbgcout<<"Norm(VVt-1) = "<<Norm(V*V.adjoint()-T(1))<<std::endl;
                    if (Norm(U*M*V-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                    if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                    if (Norm(V*V.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
                }
#endif

                // Check for deflation step 2:
                int i_firstswap = N;
                tol = eps*D(N-1);
                dbgcout<<"tol = "<<tol<<std::endl;
                for(int i=N-1;i>0;--i) {
                    dbgcout<<"D("<<i-1<<") = "<<D(i-1)<<std::endl;
                    dbgcout<<"D("<<i<<") = "<<D(i)<<std::endl;
                    dbgcout<<"diff = "<<TMV_ABS(D(i) - D(i-1))<<std::endl;
                    if (TMV_ABS(D(i) - D(i-1)) < tol) {
                        dbgcout<<"Deflate this pair\n";
                        Givens<RT> G = GivensRotate(z(i-1),z(i));
                        TMVAssert(z(i) == RT(0));
                        if (V.cptr()) {
                            typename Mv::rowpair_type Vrp = V.rowPair(i-1,i);
                            G.mult(Vrp);
                        }
                        if (i > 1) {
                            if (U.cptr()) {
                                typename Mu::colpair_type::transpose_type Ucpt = 
                                    U.colPair(i-1,i).transpose();
                                G.mult(Ucpt);
                            }
#ifdef XDEBUG_SVD
                            MatrixView<T> Mcpt = M.colPair(i-1,i).transpose();
                            MatrixView<T> Mrp = M.rowPair(i-1,i);
                            G.mult(Mcpt);
                            G.mult(Mrp);
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
#ifdef XDEBUG_SVD
                            MatrixView<T> Mcpt = M.colPair(0,1).transpose();
                            G.mult(Mcpt);
#endif
                        }
                        if (i < N-1) {
                            z.swap(i,N-1);
                            D.swap(i,N-1);
                            if (U.cptr()) U.swapCols(i,N-1);
                            if (V.cptr()) V.swapRows(i,N-1);
#ifdef XDEBUG_SVD
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
                dbgcout<<"After second deflation\n";
#ifdef XDEBUG_SVD
                dbgcout<<"DN = "<<D.subVector(0,N)<<std::endl;
                dbgcout<<"zN = "<<z.subVector(0,N)<<std::endl;
                if (U.cptr() && V.cptr()) {
                    //dbgcout<<"M = "<<M<<std::endl;
                    //dbgcout<<"U = "<<U<<std::endl;
                    //dbgcout<<"V = "<<V<<std::endl;
                    //dbgcout<<"UMV = "<<U*M*V<<std::endl;
                    dbgcout<<"Norm(UMV-A0) = "<<Norm(U*M*V -A0)<<std::endl;
                    dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
                    dbgcout<<"Norm(VVt-1) = "<<Norm(V*V.adjoint()-T(1))<<std::endl;
                    if (Norm(U*M*V-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                    if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                    if (Norm(V*V.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
                }
#endif

                // Re-sort if we screwed up the sort in the second deflation.
                // The D(i) for i >= i_firstswap might be unsorted
                dbgcout<<"i_firstswap = "<<i_firstswap<<std::endl;
                dbgcout<<"N-i_firstswap = "<<N-i_firstswap<<std::endl;
                if (i_firstswap < N-1) {
                    P.resize(N-i_firstswap);
                    D.subVector(i_firstswap,N).sort(P,Ascend);
                    z.subVector(i_firstswap,N) = 
                        P * z.subVector(i_firstswap,N);
                    if (U.cptr()) U.colRange(i_firstswap,N) =
                        U.colRange(i_firstswap,N) * P.transpose();
                    if (V.cptr()) V.rowRange(i_firstswap,N) =
                        P * V.rowRange(i_firstswap,N);
                }
                dbgcout<<"After second sort\n";
#ifdef XDEBUG_SVD
                dbgcout<<"DN = "<<D.subVector(0,N)<<std::endl;
                dbgcout<<"zN = "<<z.subVector(0,N)<<std::endl;
                if (U.cptr() && V.cptr()) {
                    if (i_firstswap < N-1) {
                        M.colRange(i_firstswap,N) =
                            M.colRange(i_firstswap,N) * P.transpose();
                        M.rowRange(i_firstswap,N) =
                            P * M.rowRange(i_firstswap,N);
                    }
                    //dbgcout<<"M = "<<M<<std::endl;
                    //dbgcout<<"UMV = "<<U*M*V<<std::endl;
                    dbgcout<<"Norm(UMV-A0) = "<<Norm(U*M*V -A0)<<std::endl;
                    dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
                    dbgcout<<"Norm(VVt-1) = "<<Norm(V*V.adjoint()-T(1))<<std::endl;
                    if (Norm(U*M*V-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                    if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                    if (Norm(V*V.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
                }
#endif
            }

            // After deflation, it is possible for N to be 1, in which case,
            // the SVD is done.
            if (N > 1) {
                typename Vd::subvector_type DN = D.subVector(0,N);
                VectorView<RT,Unit> zN = z.subVector(0,N);

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

                if (U.cptr() || V.cptr()) {
                    Matrix<RT,ColMajor> W(N,N); // W(i,j) = D(i)-S(j)
                    FindDCSingularValues(S,RT(1),DN,zN,W);
                    dbgcout<<"S = "<<S<<std::endl;

                    // Update z's to numerically more stable values:
                    // z_i' = sqrt( (sn^2-di^2) 
                    //             prod_1..i-1 (sk^2-di^2)/(dk^2-di^2) *
                    //             prod_i..N-1 (sk^2-di^2)/(dk+1^2-di^2) ) 
                    //        * sign(z_i)
                    for(int i=0;i<N;i++) {
                        RT di = D(i);
                        RT prod = -W(i,N-1)*(S(N-1) + di);
                        for(int k=0;k<i;k++) {
                            prod *= -W(i,k)/(D(k)-di);
                            prod *= (S(k)+di)/(D(k)+di);
                        }
                        for(int k=i+1;k<N;k++) {
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
                    for(int j=0;j<N;j++) {
                        VectorView<RT,Unit> yj = W.col(j);
                        Vector<RT> diff_j = yj; // copy current values
                        yj(0) = -z(0)/(S(j)*S(j));
                        for(int i=1;i<N;i++)
                            yj(i) = (z(i) / diff_j(i)) / (D(i)+S(j));
                    }
                    if (V.cptr()) {
                        Vector<RT> normyj(N);
                        for(int j=0;j<N;j++) {
                            W.col(j) /= (normyj(j) = Norm(W.col(j)));
                        }
                        // V = Y * V
                        V.rowRange(0,N) = W.transpose() * V.rowRange(0,N);
                        if (U.cptr()) for(int j=0;j<N;j++) W.col(j) *= normyj(j);
                    }
                    if (U.cptr()) {
                        // Now convert W from Y.transpose() to X
                        W.row(0).setAllTo(RT(-1));
                        for(int i=1;i<N;i++) W.row(i) *= D(i);
                        for(int j=0;j<N;j++) W.col(j) /= Norm(W.col(j));
                        // U = U*X
                        U.colRange(0,N) *= W;
                    }
                } else {
                    FindDCSingularValues(S,RT(1),DN,zN);
                    dbgcout<<"S = "<<S<<std::endl;
                }

                // Done.  Copy S to D and undo scaling.
                DN = scale * S;

                dbgcout<<"Done D = "<<D<<std::endl;
#ifdef XDEBUG_SVD
                if (U.cptr() && V.cptr()) {
                    M.subMatrix(0,N,0,N) = DiagMatrixViewOf(DN);
                    //dbgcout<<"M = "<<M<<std::endl;
                    //dbgcout<<"UMV = "<<U*M*V<<std::endl;
                    dbgcout<<"Norm(UMV-A0) = "<<Norm(U*M*V -A0)<<std::endl;
                    dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
                    dbgcout<<"Norm(VVt-1) = "<<Norm(V*V.adjoint()-T(1))<<std::endl;
                    if (Norm(U*M*V-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                    if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                    if (Norm(V*V.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
                }
#endif
            } else if (N==1) {
                if (D(0) == RT(0)) {
                    D(0) = TMV_ABS(z(0));
                    if (V.cptr() && z(0) < RT(0)) V.row(0) = -V.row(0);
                }
#ifdef XDEBUG_SVD
                if (U.cptr() && V.cptr()) {
                    M(0,0) = D(0);
                    //dbgcout<<"M = "<<M<<std::endl;
                    //dbgcout<<"UMV = "<<U*M*V<<std::endl;
                    dbgcout<<"Norm(UMV-A0) = "<<Norm(U*M*V -A0)<<std::endl;
                    dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
                    dbgcout<<"Norm(VVt-1) = "<<Norm(V*V.adjoint()-T(1))<<std::endl;
                    if (Norm(U*M*V-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                    if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                    if (Norm(V*V.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
                }
#endif
            }
            dbgcout<<"Done SVDecompose algo 11: D => "<<D<<std::endl;
#ifdef XDEBUG_SVD
            if (U.cptr() && V.cptr()) {
                M = DiagMatrixViewOf(D);
                //dbgcout<<"M = "<<M<<std::endl;
                //dbgcout<<"UMV = "<<U*M*V<<std::endl;
                dbgcout<<"Norm(UMV-A0) = "<<Norm(U*M*V -A0)<<std::endl;
                dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
                dbgcout<<"Norm(VVt-1) = "<<Norm(V*V.adjoint()-T(1))<<std::endl;
                if (Norm(U*M*V-A0) > THRESH*normA0) abort();
#ifdef TESTUV
                if (Norm(U.adjoint()*U-T(1)) > THRESH*normA0) abort();
                if (Norm(V*V.adjoint()-T(1)) > THRESH*normA0) abort();
#endif
            }
#endif

#ifdef XDEBUG_SVD
            DiagMatrix<RT> S(D);
            Matrix<T> A2(A0.colsize(),A0.rowsize());
            dbgcout<<"S = "<<S.diag()<<std::endl;
            if (U.cptr() && V.cptr()) {
                A2 = U*S*V;
                dbgcout<<"Norm(A2-A0) = "<<Norm(A2-A0)<<std::endl;
                dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
                dbgcout<<"Norm(VVt-1) = "<<Norm(V*V.adjoint()-T(1))<<std::endl;
                dbgcout<<"THRESH * normA0 = "<<THRESH<<" * "<<normA0<<" = "<<THRESH*normA0<<std::endl;
                if (!(Norm(A2-A0) <= THRESH*normA0)
#ifdef TESTUV
                    || ( V.colsize()>=V.rowsize() &&
                         !(Norm(V.adjoint()*V-T(1)) <= THRESH*normA0) )
                    || ( V.rowsize()>=V.colsize() &&
                         !(Norm(V*V.adjoint()-T(1)) <= THRESH*normA0) )
                    || ( U.colsize()>=U.rowsize() &&
                         !(Norm(U.adjoint()*U-T(1)) <= THRESH*normA0) )
#endif
                ) {
                    std::cerr<<"SV_Decompose_DC: \n";
                    std::cerr<<"input D = "<<D0<<std::endl;
                    std::cerr<<"input E = "<<E0<<std::endl;
                    std::cerr<<"output S = "<<D<<std::endl;
                    Matrix<T> U_QR = U0;
                    Matrix<T> V_QR = V0;
                    Vector<RT> D_QR = D0;
                    Vector<RT> E_QR = E0;
                    SV_DecomposeFromBidiagonal_QR(U_QR,D,E,V,UisI,VisI);
                    std::cerr<<"QR solution: S = "<<D_QR<<std::endl;
                    std::cerr<<"U = "<<U<<std::endl;
                    std::cerr<<"QR solution: U = "<<U_QR<<std::endl;
                    std::cerr<<"U-U_QR = "<<U-U_QR<<std::endl;
                    std::cerr<<"Norm(U-U_QR) = "<<Norm(U-U_QR)<<std::endl;
                    std::cerr<<"V = "<<V<<std::endl;
                    std::cerr<<"QR solution: V = "<<V_QR<<std::endl;
                    std::cerr<<"V-V_QR = "<<V-V_QR<<std::endl;
                    std::cerr<<"Norm(V-V_QR) = "<<Norm(V-V_QR)<<std::endl;
                    std::cerr<<"UBV = "<<A0<<std::endl;
                    std::cerr<<"USV = "<<A2<<std::endl;
                    std::cerr<<"Norm(diff) = "<<Norm(A2-A0)<<std::endl;
                    std::cerr<<"cf. Norm(A) = "<<normA0<<std::endl;
                    std::cerr<<"THRESH * Norm(A) = "<<THRESH*normA0<<std::endl;
                    abort();
                }
            }
#endif
        }
    };

    // TODO: Make OpenMP algorithm.
    // A note about OpenMP here:
    // Both the divide step (the recursion in algo 11)
    // and the conquer step (the loop in FindDCSingularValues) are
    // pontentially parallelizable.  However, recursive functions do not
    // lend themselves to being parallelized by OpenMP very easily.
    // So I think I need to rewrite the divide step in a way that does
    // not rely on recursion before I apply OPENMP pragmas to it.
    // (Possibly use recursion to make a list of problems to do, and
    // then have OpenMP parallelize the list of tasks at each level.)
    // Anyway, for now we only have recursion in the FindDCSingularValues 
    // function.
    //
    // Update: The new OpenMP 3.0 task feature may make this a lot easier.
    // Try that first!

    // algo 90: call InstSV_DecomposeFromBidiagonal_DC
    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_DC_Helper<90,cs,rs,Mu,Vd,Ve,Mv>
    {
        static TMV_INLINE void call(
            Mu& U, Vd& D, Ve& E, Mv& V, bool UisI, bool VisI)
        { 
            InstSV_DecomposeFromBidiagonal_DC(
                U.xView(),D.xView(),E.xView(),V.xView(),UisI,VisI); 
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_DC_Helper<-3,cs,rs,Mu,Vd,Ve,Mv>
    {
        static TMV_INLINE void call(
            Mu& U, Vd& D, Ve& E, Mv& V, bool UisI, bool VisI)
        {
            const int algo = 11;
#ifdef PRINTALGO_SVD
            std::cout<<"Inline SVDecomposeFromBidiagonal: \n";
            std::cout<<"U = "<<TMV_Text(U)<<std::endl;
            std::cout<<"D = "<<TMV_Text(D)<<std::endl;
            std::cout<<"E = "<<TMV_Text(E)<<std::endl;
            std::cout<<"V = "<<TMV_Text(V)<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<std::endl;
            std::cout<<"sizes = "<<U.colsize()<<"  "<<U.rowsize()<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
#ifdef XDEBUG_SVD
            typedef typename Traits2<typename Mu::value_type, typename Mv::value_type>::type T;
            typedef typename Mu::real_type RT;
            TMVAssert(D.minAbsElement() > RT(0));
            TMVAssert(E.minAbsElement() > RT(0));

            const int N = rs==TMV_UNKNOWN ? int(D.size()) : rs;
            dbgcout<<"Start Decompose from Bidiag:\n";
            if (U.cptr()) dbgcout<<"U = "<<TMV_Text(U)<<std::endl;
            if (V.cptr()) dbgcout<<"V = "<<TMV_Text(V)<<std::endl;
            dbgcout<<"D = "<<TMV_Text(D)<<"  step "<<D.step()<<"  "<<D<<std::endl;
            dbgcout<<"E = "<<TMV_Text(E)<<"  step "<<E.step()<<"  "<<E<<std::endl;
            //if (U.cptr()) dbgcout<<"U = "<<U<<std::endl;
            //if (V.cptr()) dbgcout<<"V = "<<V<<std::endl;

            dbgcout<<"UisI, VisI = "<<UisI<<"  "<<VisI<<std::endl;
            Matrix<RT> B(N,N,RT(0));
            B.diag() = D;
            B.diag(1) = E;
            const int M1 = U.cptr() && V.cptr() ? int(U.colsize()) : 0;
            const int N1 = U.cptr() && V.cptr() ? int(V.rowsize()) : 0;
            Matrix<T> A0(M1,N1);
            if (U.cptr() && V.cptr()) A0 = U*B*V;
            //dbgcout<<"A0 = "<<A0<<std::endl;
#endif
            SVDecomposeFromBidiagonal_DC_Helper<algo,cs,rs,Mu,Vd,Ve,Mv>::call( 
                U,D,E,V,UisI,VisI);
            //std::cout<<"U = "<<U<<std::endl;
            //std::cout<<"S = "<<S<<std::endl;
            //std::cout<<"V = "<<V<<std::endl;
#ifdef XDEBUG_SVD
            if (U.cptr() && V.cptr()) {
                Matrix<T> AA = U * DiagMatrixViewOf(D) * V;
                dbgcout<<"Done DC Norm(A0-AA) = "<<Norm(A0-AA)<<std::endl;
                dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
                dbgcout<<"Norm(VtV-1) = "<<Norm(V.adjoint()*V-T(1))<<std::endl;
                dbgcout<<"Norm(VVt-1) = "<<Norm(V*V.adjoint()-T(1))<<std::endl;
                dbgcout<<"U = "<<U<<std::endl;
                dbgcout<<"S = "<<D<<std::endl;
                dbgcout<<"V = "<<V<<std::endl;
                if (!(Norm(A0-AA) <= THRESH*Norm(A0))) {
                    std::cerr<<"SV_DecomposeFromBidiagonal DC: \n";
                    std::cerr<<"UBV = "<<A0<<std::endl;
                    std::cerr<<"USV = "<<AA<<std::endl;
                    //std::cerr<<"U = "<<U<<std::endl;
                    //std::cerr<<"V = "<<V<<std::endl;
                    std::cerr<<"input B = "<<B<<std::endl;
                    std::cerr<<"S = "<<D<<std::endl;
                    dbgcout<<"Norm(UBV-USV) = "<<Norm(A0-AA)<<std::endl;
                    dbgcout<<"cf "<<THRESH<<"*Norm(UBV) = "<<THRESH*Norm(A0)<<std::endl;
                    abort();
                }
            }
#endif
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_DC_Helper<-2,cs,rs,Mu,Vd,Ve,Mv>
    {
        static TMV_INLINE void call(
            Mu& U, Vd& D, Ve& E, Mv& V, bool UisI, bool VisI)
        {
            typedef typename Mu::value_type T;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                Traits<T>::isinst;
            const int algo = 
                inst ? 90 :
                -3;
            SVDecomposeFromBidiagonal_DC_Helper<algo,cs,rs,Mu,Vd,Ve,Mv>::call(
                U,D,E,V,UisI,VisI);
        }
    };

    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_DC_Helper<-1,cs,rs,Mu,Vd,Ve,Mv>
    {
        static TMV_INLINE void call(
            Mu& U, Vd& D, Ve& E, Mv& V, bool UisI, bool VisI)
        {
            SVDecomposeFromBidiagonal_DC_Helper<-2,cs,rs,Mu,Vd,Ve,Mv>::call(
                U,D,E,V,UisI,VisI); 
        }
    };

    template <class Mu, class Vd, class Ve, class Mv>
    inline void InlineSV_DecomposeFromBidiagonal_DC(
        BaseMatrix_Rec_Mutable<Mu>& U,
        BaseVector_Mutable<Vd>& D, BaseVector_Mutable<Ve>& E, 
        BaseMatrix_Rec_Mutable<Mv>& V, bool UisI, bool VisI)
    {
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Vd::value_type>::samebase));
        TMVStaticAssert((Traits2<
                         typename Vd::value_type,
                         typename Ve::value_type>::sametype));
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Mv::value_type>::samebase));
        TMVStaticAssert(Traits<typename Vd::value_type>::isreal);
        TMVStaticAssert(!Mu::_conj);
        TMVStaticAssert(!Mv::_conj);
        TMVStaticAssert((Sizes<Mu::_rowsize,Vd::_size>::same));
        TMVStaticAssert((Sizes<Vd::_size,IntTraits<Ve::_size>::Sp1>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_colsize>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_rowsize>::same));
        TMVAssert(D.size() == E.size()+1);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        if (V.cptr()) {
            TMVAssert(V.colsize() == D.size());
            TMVAssert(V.rowsize() >= V.colsize());
        }
        const int cs = Mu::_colsize;
        const int rs1 = Sizes<Mu::_rowsize,Vd::_size>::size;
        const int rs2 = Sizes<Mv::_rowsize,Mv::_colsize>::size;
        const int rs = Sizes<rs1,rs2>::size;
        typedef typename Mu::cview_type Muv;
        typedef typename Vd::cview_type Vdv;
        typedef typename Ve::cview_type Vev;
        typedef typename Mv::cview_type Mvv;
        TMV_MAYBE_REF(Mu,Muv) Uv = U.cView();
        TMV_MAYBE_REF(Vd,Vdv) Dv = D.cView();
        TMV_MAYBE_REF(Ve,Vev) Ev = E.cView();
        TMV_MAYBE_REF(Mv,Mvv) Vv = V.cView();
        SVDecomposeFromBidiagonal_DC_Helper<-3,cs,rs,Muv,Vdv,Vev,Mvv>::call(
            Uv,Dv,Ev,Vv,UisI,VisI);
    }

    template <class Mu, class Vd, class Ve, class Mv>
    inline void SV_DecomposeFromBidiagonal_DC(
        BaseMatrix_Rec_Mutable<Mu>& U,
        BaseVector_Mutable<Vd>& D, BaseVector_Mutable<Ve>& E, 
        BaseMatrix_Rec_Mutable<Mv>& V, bool UisI, bool VisI)
    {
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Vd::value_type>::samebase));
        TMVStaticAssert((Traits2<
                         typename Vd::value_type,
                         typename Ve::value_type>::sametype));
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Mv::value_type>::samebase));
        TMVStaticAssert(Traits<typename Vd::value_type>::isreal);
        TMVStaticAssert(!Mu::_conj);
        TMVStaticAssert(!Mv::_conj);
        TMVStaticAssert((Sizes<Mu::_rowsize,Vd::_size>::same));
        TMVStaticAssert((Sizes<Vd::_size,IntTraits<Ve::_size>::Sp1>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_colsize>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_rowsize>::same));
        TMVAssert(D.size() == E.size()+1);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        if (V.cptr()) {
            TMVAssert(V.colsize() == D.size());
            TMVAssert(V.rowsize() >= V.colsize());
        }
        const int cs = Mu::_colsize;
        const int rs1 = Sizes<Mu::_rowsize,Vd::_size>::size;
        const int rs2 = Sizes<Mv::_rowsize,Mv::_colsize>::size;
        const int rs = Sizes<rs1,rs2>::size;
        typedef typename Mu::cview_type Muv;
        typedef typename Vd::cview_type Vdv;
        typedef typename Ve::cview_type Vev;
        typedef typename Mv::cview_type Mvv;
        TMV_MAYBE_REF(Mu,Muv) Uv = U.cView();
        TMV_MAYBE_REF(Vd,Vdv) Dv = D.cView();
        TMV_MAYBE_REF(Ve,Vev) Ev = E.cView();
        TMV_MAYBE_REF(Mv,Mvv) Vv = V.cView();
        SVDecomposeFromBidiagonal_DC_Helper<-2,cs,rs,Muv,Vdv,Vev,Mvv>::call(
            Uv,Dv,Ev,Vv,UisI,VisI);
    }

} // namespace tmv

#undef TMV_DC_LIMIT

#endif

