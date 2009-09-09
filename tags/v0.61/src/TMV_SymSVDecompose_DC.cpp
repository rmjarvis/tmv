///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_Blas.h"
#include "TMV_SymSVDiv.h"
#include "TMV_SymMatrix.h"
#include "TMV_Givens.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"
#include <iostream>
using std::endl;

//#define XDEBUG
//#define TIME

#ifdef XDEBUG
#include "TMV_DiagMatrixArith.h"
//#define dbgcout std::cout
#define dbgcout if (false) std::cout
using std::cerr;
#else
#define dbgcout if (false) std::cout
#endif

namespace tmv {

#define RT RealType(T)

#ifdef TMV_BLOCKSIZE
#define DC_LIMIT TMV_BLOCKSIZE/2
#else
#define DC_LIMIT 2
#endif

  template <class T> static T FindDCEigenValue(
      const int k, const int N, const T rho, const T* D, const T* z, 
      const T* zsq, const T normsqz, T* diff)
  {
    // Find solutions to the equation:
    // f(s) = rho + Sum_i=1..N z_i^2/(D_i-s)
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
    const T eps = Epsilon<T>();
    const int MAXITER = 20;

    dbgcout<<"k = "<<k<<endl;
    if (k<N-1) {
      dbgcout<<"root bracketed by "<<D[k]<<"  "<<D[k+1]<<endl;
      dbgcout<<"   range = "<<D[k+1]-D[k]<<endl;
    } else {
      dbgcout<<"root bracketed by "<<D[k]<<"  "<<D[k]+normsqz<<endl;;
      dbgcout<<"   range = "<<normsqz<<endl;
    }

    // Most of the routine is the same for k=N-1 as for the normal
    // case, but there are some lines where the normal case uses
    // k, but when k==N-1, I need k-1.  So I use kk for these.
    int kk = k < N-1 ? k : k-1;

    // First we need an initial guess for s_k 
    dbgcout<<"D[kk+1]-D[kk] = "<<D[kk+1]-D[kk]<<endl;
    dbgcout<<"eps*D["<<kk+1<<"] = "<<eps*D[kk+1]<<endl;
    T delta = D[kk+1] - D[kk];
    // It is known that D_k < s_k < D_k+1
    // Let s = (D_k + D_k+1)/2
    // Define tau: s = D_k + tau
    T tau = k < N-1 ? delta/T(2) : normsqz/(T(2)*rho);
    T s = D[k] + tau;

    for(int j=0;j<N;j++) { diff[j] = (D[j]-D[k])-tau; }

    // Let c = f(s) - z_k^2/(D_k-s) - z_k+1^2/(D_k+1-s)
    // i.e. c is f(s) without the two terms for the poles that
    // bracket the root we are seraching for.
    // Since the largest terms are normally the ones for the poles nearer
    // to s, we perform the calculation from each end and proceed towards k.
    // Note: for k==N-1, the two poles do not bracket the root, but 
    // we proceed similarly anyway.
    T psi(0);
    for(int j=0;j<kk;j++) {
      psi += zsq[j] / diff[j];
      dbgcout<<"psi += "<<zsq[j]<<" / "<<diff[j]<<" = "<<psi<<endl;
    }
    T phi(0);
    for(int j=N-1;j>kk+1;j--) {
      phi += zsq[j] / diff[j];
      dbgcout<<"phi += "<<zsq[j]<<" / "<<diff[j]<<" = "<<phi<<endl;
    }
    T c = rho + psi + phi;
    dbgcout<<"c = "<<c<<endl;

    // Finish calculating f = f(s)
    T f = c + zsq[kk] / diff[kk] + zsq[kk+1] / diff[kk+1];
    dbgcout<<"f = c + "<<zsq[kk]<<" / "<<diff[kk]<<"\n";
    dbgcout<<"      + "<<zsq[kk+1]<<" / "<<diff[kk+1]<<"\n";
    dbgcout<<"f("<<tau<<") = "<<f<<endl;

    T lowerbound, upperbound; // the allowed range for solution tau
    int k1; // The index of the pole closest to the root
    if (k < N-1) {
      dbgcout<<"k < N-1\n";
      if (f >= T(0)) { // D_k < s_k < s
	dbgcout<<"f >= 0\n";
	k1 = k;
	lowerbound = T(0);
	upperbound = tau;
	dbgcout<<"bounds = 0,tau = "<<tau<<endl;
	T a = c*delta + zsq[k] + zsq[k+1];
	T b = zsq[k]*delta;
	T d = a*a-T(4)*b*c;
	if (d < T(0)) d = T(0); // Then rounding error.  Use d = 0
	else d = SQRT(d);
	// Our initial estimate for s_k^2 is D_k^2 + tau
	tau = (a > 0) ? T(2)*b/(a+d) : (a-d)/(T(2)*c);
	if (tau > upperbound) tau = upperbound;
      } else {
	dbgcout<<"f < 0\n";
	k1 = k+1;
	lowerbound = -tau;
	upperbound = T(0);
	dbgcout<<"bounds = -tau,0 = "<<-tau<<endl;
	T a = -c*delta + zsq[k] + zsq[k+1];
	T b = -zsq[k+1]*delta;
	T d = a*a-T(4)*b*c;
	if (d < T(0)) d = T(0); 
	else d = SQRT(d);
	// Our initial estimate for s_k^2 is D_k+1^2 + tau (where tau < 0)
	tau = (a > 0) ? T(2)*b/(a+d) : (a-d)/(T(2)*c);
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
	  dbgcout<<"bounds = tau,2tau = "<<tau<<"  "<<normsqz/rho<<endl;
	  if (c <= T(0)) tau = upperbound; break;
	  // In this case there is no positive solution.  Best guess
	  // is just to start at upper bound.
	}
	T a = -c*delta + zsq[k-1] + zsq[k];
	T b = -zsq[k]*delta;
	T d = a*a-T(4)*b*c;
	if (d < T(0)) d = T(0); 
	else d = SQRT(d);
	tau = (a < 0) ? T(2)*b/(a-d) : (a+d)/(T(2)*c);
	if (tau > upperbound) tau = upperbound;
	// The other times we check this are basically for rounding error
	// reasons.  But this time, if f <  0, the solution to the quadratic
	// can really be much > upperbound.  So the check is especially
	// important in this case.
      } while (false);
    }
    dbgcout<<"Initial tau = "<<tau<<endl;
    dbgcout<<"Allowed bounds for tau = "<<lowerbound<<"  "<<upperbound<<endl;
    TMVAssert(tau >= lowerbound && tau <= upperbound);

    s = D[k1] + tau;
    for(int j=0;j<N;j++) diff[j] = D[j]-D[k1]; 
    for(int j=0;j<N;j++) diff[j] -= tau; 

    // Define f(s) = rho + psi(s) + z_k1^2/(D_k1^2-s^2) + phi(s)
    // psi(s) = Sum_k=1..k1-1 z_k^2/(D_k^2-s^2)
    // phi(s) = Sum_k=k1+1..N z_k^2/(D_k^2-s^2)

    T e, dpsi, dphi, fk, dfk, k1term, dk1term, df;
    // e is used for an error bound on f, which we will accumulate
    // as we go.  The full formula is :
    // e = 2rho + Sum_j=1..k (k-j+6) z_j^2/abs(D_j-s) 
    //       + Sum_j=k+1..N (j-k+5) z_j^2/(D_j-s)
    //       + abs(f(s)) + abs(tau)*abs(f'(s))

    // psix, phix store the value of psi,phi without the last term added.
    // These are useful in the threepoles section below.
    T psix=0, phix=0, dpsix=0, dphix=0;

    T f1 = f;
    for(int iter = 0; iter < MAXITER; iter++) {
      // Calculate psi, phi
      psi = phi = e = dpsi = dphi = T(0);
      for(int j=0;j<k1;j++) {
	T temp = z[j] / diff[j];
	psix = psi; dpsix = dpsi;
	psi += z[j] * temp;
	dpsi += temp * temp;
	e -= psi;
      }
      for(int j=N-1;j>k1;j--) {
	T temp = z[j] / diff[j];
	phix = phi; dphix = dphi;
	phi += z[j] * temp;
	dphi += temp * temp;
	e += phi;
      }

      // Finish calculating f
      fk = rho + psi + phi; // == f(s) without z_k1^2/(D_k1-s) term
      dfk = dpsi + dphi;
      T temp = z[k1] / diff[k1];
      k1term = temp*z[k1];
      dk1term = temp*temp;
      f = fk + k1term;
      df = dfk + dk1term;
      e += T(2)*rho + T(6)*(phi-psi + ABS(k1term)) + ABS(f) + ABS(tau)*df;
      dbgcout<<"fk = "<<fk<<endl;
      dbgcout<<"k1 = "<<k1<<"  k = "<<k<<endl;
      dbgcout<<"f = "<<f<<endl;
      dbgcout<<"psi = "<<psi<<endl;
      dbgcout<<"k1term = "<<k1term<<endl;
      dbgcout<<"phi = "<<phi<<endl;

      if (k == N-1) break;
      if (ABS(f) <= eps*e) break; 
      if (ABS(f) <= RT(1.e3)*ABS(f1) && f*f1 < T(0)) break;

      // If we get to here, then the initial guess was not very good.  
      // This can mean that the function is not very well described
      // by the nearby poles.  So use bisection until we at least do better
      // than the midpoint guess we started with.
      // Then we can move on to the fancy, quick-converging techniques below.
      // This is quite rare, but when the above calculation gives a very bad
      // initial guess, this adjustment can be critical.
      dbgcout<<"f = "<<f<<", f1 = "<<f1<<endl;
      dbgcout<<"Use bisection to update initial tau:\n";

      T t1 = tau;

      if (f < T(0)) { lowerbound = tau; tau = (upperbound+tau)/T(2); }
      else { upperbound = tau; tau = (lowerbound+tau)/T(2); }
      dbgcout<<"tau -> "<<tau<<endl;

      s = D[k1] + tau;
      T dt = tau-t1; 
      for(int j=0;j<N;j++) diff[j] -= dt;

#ifndef NWARN
      if (iter == MAXITER-1) {
	std::cout<<"Warning - Unable to find appropriate initial guess\n";
	std::cout<<"in FindDCEigenValue\n";
	std::cout<<"Unseccessful after "<<MAXITER<<" iterations.\n";
      }
#endif
    }

    // Use 3 poles when fk<0 for k1=k or fk>0 for k1=k+1
    bool threepoles = (((fk < T(0)) == (k1 == k)) && k1>0 && k1<N-1);
    bool fixedweight = true;
    if (threepoles) dbgcout<<"USE THREE POLES\n";

    for(int iter = 0; iter < MAXITER; iter++) {
      dbgcout<<"Main iter = "<<iter<<", f = "<<f<<", eps*e = "<<eps*e<<endl;
      dbgcout<<"f("<<tau<<") = "<<f<<"  s = "<<s<<endl;

      if (ABS(f) <= eps*e) break;

      T eta; // = the change to be added to s
      if (!threepoles) { // The normal case - only 2 poles used for approx.
	T d1 = diff[kk];
	T d2 = diff[kk+1];
	dbgcout<<"d1,d2 = "<<d1<<", "<<d2<<endl;
	if (fixedweight) 
	  if (kk==k1) c = fk - d2*dfk;
	  else c = fk - d1*dfk;
	else  c = fk - d1*dpsi - d2*dphi;
	// c eta^2 - a eta + b = 0
	T a = (d1+d2)*f - d1*d2*df;
	T b = d1*d2*f;
	dbgcout<<"c,a,b = "<<c<<", "<<a<<", "<<b<<endl;
	if (c == T(0)) {
	  if (a == T(0)) {
	    // Then the above has rounding error.
	    // This version can't have a == 0:
	    if (fixedweight)
	      // For k1 == k case:
	      // a = (d1+d2)*f - d1*d2*df
	      // with f = d2*df - delta*dk1term  (i.e. c==0)
	      // a = (d1+d2)*d2*df - (d1+d2)*(d2-d1)*dk1term - d1*d2*df
	      //   = d2^2*df + (d1+d2)*(d1-d2)*dk1term
	      //   = d2^2*(df-dk1term) + d1^2 zsq(k1)/d1^2
	      // For k1 == k+1 case, swap d1,d2
	      a = zsq[k1] + (k1==kk ? d2*d2 : d1*d1)*dfk;
	    else
	      // a = (d1+d2)*f - d1*d2*df
	      // with f = d1*dpsi + d2*dphi + k1term
	      // a = (d1+d2)*d1*dpsi + (d1+d2)*d2*dphi - d1*d2*(dpsi+dphi) 
	      //           + (d1+d2) k1term - d1*d2*dk1term
	      //   = d1*d1*dpsi + d2*d2*dphi + d1 k1term
	      a = d1*d1*dpsi + d2*d2*dphi + (k1==kk ? d1 : d2)*k1term;
	  }
	  eta = b/a;
	  dbgcout<<"c=0  eta = "<<eta<<endl;
	  if (k==N-1)
	    TMVAssert(eta > -tau);
	  else
	    TMVAssert(eta > d1 && eta < d2);
	} else {
	  T d = a*a-T(4)*b*c;
	  dbgcout<<"d = "<<d<<endl;
	  if (d < T(0)) d = T(0); 
	  else d = SQRT(d);
	  if (k==N-1) { // Then we want the solutionn with other sign for d
	    eta = a >= 0 ? (a+d)/(T(2)*c) : T(2)*b/(a-d);
	    dbgcout<<"eta = "<<eta<<endl;
	    if (eta <= -tau) // Then rounding errors - use bisect
	      eta = -tau/RT(2);
	    TMVAssert(eta > -tau);
	  }
	  else {
	    eta = a <= 0 ? (a-d)/(T(2)*c) : T(2)*b/(a+d);
	    dbgcout<<"eta = "<<eta<<endl;
	    if (eta <= d1) // Then rounding errors - use bisect
	      eta = d1/RT(2);
	    if (eta >= d2) // Then rounding errors - use bisect
	      eta = d2/RT(2);
	    dbgcout<<"eta => "<<eta<<endl;
	    TMVAssert(eta > d1 && eta < d2);
	  }
	}
      } else { // More complicated case using 3 poles. 
	dbgcout<<"threepoles section\n";
	T d1 = diff[k1-1];
	T d2 = diff[k1+1];
	dbgcout<<"d1,d2 = "<<d1<<" "<<d2<<endl;
	if (fixedweight) 
	  if (k1 == k) c = (rho+psix+phi) - d2*(dpsix+dphi);
	  else c = (rho+psi+phix) - d1*(dpsi+dphix);
	else 
	  if (k1 == k) c = (rho+psix+phi) - d1*dpsix - d2*dphi;
	  else c = (rho+psi+phix) - d1*dpsi - d2*dphix;
	dbgcout<<"c = "<<c<<endl;
	// Cubic equation is h(eta) = c eta^3 - a eta^2 + b eta + g = 0
	// with a = {fk(d1+d2) - dfk d1 d2} + z^2[k1] - c tau
	//      b = d1 d2 (fk + tau dfk) - f (d1+d2) tau 
	//      b = d1 d2 fk - tau ( {fk(d1+d2) - dfk d1 d2} + (d1+d2) k1term)
	//      g = d1 d2 tau f
	// We need to be careful with the terms in {} above (in a and b),
	// since the terms can nearly cancel.  Do better with phix,psix:
	// For k1 == k:
	//      {} = (rho + psix + zsq[k1-1]/d1 + phi)(d1+d2) 
	//           - (dpsix + zsq[k1-1]/d1^2 + dphi)(d1 d2)
	//         = (rho+psix+phi)(d1+d2) - (dpsix+dphi)(d1 d2) + zsq[k1-1]
	T temp = (k1==k ?
	    (d1+d2)*(rho+psix+phi) - d1*d2*(dpsix+dphi) + zsq[k1-1] :
	    (d1+d2)*(rho+psi+phix) - d1*d2*(dpsi+dphix) + zsq[k1+1]);
	dbgcout<<"temp = "<<temp<<endl;
	T a = temp + zsq[k1] - c*tau;
	dbgcout<<"a = "<<a<<endl;
	T b = d1*d2*fk - tau*(temp + (d1+d2)*k1term);
	dbgcout<<"b = "<<b<<endl;
	T g = d1*d2*tau*f;
	dbgcout<<"g = "<<g<<endl;
	T mineta = lowerbound - tau;
	T maxeta = upperbound - tau;
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
	  // The best guess then is to just go to eta2.
	  // But instead we use bisection to be safe.
	  eta = eta2/2;
	} else for(int iter3 = 0; iter3 < MAXITER; iter3++) {
	  dbgcout<<"iter3 = "<<iter3<<", eta = "<<eta<<"  h = "<<h<<endl;
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
	  T x = T(1) - T(2) *(h/dh)*(d2h/dh);
	  dbgcout<<"dh,d2h,x = "<<dh<<", "<<d2h<<", "<<x<<endl;
	  x = (x <= T(0)) ? T(0) : SQRT(x);
	  T deta = -T(2)*h / (dh*(T(1) + x));
	  dbgcout<<"deta = "<<deta<<endl;
	  if ((eta2 > eta && (deta < 0 || deta > eta2-eta)) || 
	      (eta2 < eta && (deta > 0 || deta < eta2-eta)) ) {
	    dbgcout<<"deta is wrong direction or too big - use other root.\n";
	    deta = -dh*(T(1) + x) / d2h;
	    dbgcout<<"deta = "<<deta<<endl;
	    if ((eta2 > eta && (deta < 0 || deta > eta2-eta)) || 
		(eta2 < eta && (deta > 0 || deta < eta2-eta)) ) {
	      dbgcout<<"deta is wrong direction or too big - use bisection.\n";
	      deta = (eta2-eta)/T(2);
	      dbgcout<<"deta -> "<<deta<<endl;
	    }
	  }
	  TMVAssert(eta + deta <= maxeta);
	  TMVAssert(eta + deta >= mineta);
	  T etanew = eta + deta;
	  dbgcout<<"etanew = "<<etanew<<endl;
	  T hnew = (ABS(deta) < RT(1.e-3)*ABS(eta)) ?
	    // Do delta calculation for hnew
	    // c(e+de)^3 - a(e+de)^2 + b(e+de) + g
	    // h + 3ce^2 de + 3ce de^2 + cde^3 - 2ae de - a de^2 + b de
	    // h + c de^3 + (3ce-a) de^2 + (3ce^2-2ae+b) de
	    h + ((c*deta + T(3)*c*eta-a)*deta + 
		(T(3)*c*eta-T(2)*a)*eta+b )*deta
	    :
	    // Regular calculation is ok.
	    ((c*etanew - a)*etanew + b)*etanew + g;
	  dbgcout<<"hnew = "<<hnew<<endl;
	  if ( (h > T(0)) != (hnew > T(0)) ) { eta2 = eta; h2 = h; }
	  eta = etanew;
	  h = hnew;
	  if (ABS(deta) < eps*ABS(eta) || ABS(h) < eps*g) break;
#ifndef NWARN
	  if (iter3 == MAXITER-1) {
	    std::cout<<"Warning - unable to converge for THREEPOLES solution\n";
	    std::cout<<"No solution after "<<MAXITER<<" iterations.\n";
	  }
#endif
	}
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
      s += eta;
      for(int j=0;j<N;j++) diff[j] -= eta; 

      // Update psi, phi, etc.
      psi = T(0);
      dpsi = T(0);
      e = T(0);
      for(int j=0;j<k1;j++) {
	T temp = z[j] / diff[j];
	psix = psi; dpsix = dpsi;
	psi += z[j] * temp;
	dpsi += temp * temp;
	e -= psi;
      }
      phi = T(0);
      dphi = T(0);
      for(int j=N-1;j>k1;j--) {
	T temp = z[j] / diff[j];
	phix = phi; dphix = dphi;
	phi += z[j] * temp;
	dphi += temp * temp;
	e += phi;
      }
      fk = rho + psi + phi; 
      dfk = dpsi+dphi;
      T temp = z[k1] / diff[k1];
      k1term = temp*z[k1];
      dk1term = temp*temp;
      T fnew = fk + k1term;
      df = dfk + dk1term;
      e += T(2)*rho + T(6)*(phi-psi + ABS(k1term)) + ABS(fnew) + ABS(tau)*df;

      // Use 3 poles when fk<0 for k1=k or fk>0 for k1=k+1
      bool newthreepoles = (((fk < T(0)) == (k1 == k)) && k1>0 && k1<N-1);
      if (newthreepoles != threepoles) {
	dbgcout<<"Change threepoles from "<<threepoles<<" to "<<newthreepoles<<endl;
	dbgcout<<"fk = "<<fk<<"  fnew = "<<fnew<<endl;
	threepoles = newthreepoles;
      }

      // Update bounds:
      if (fnew < T(0)) lowerbound = tau;
      else upperbound = tau;
      dbgcout<<"New bounds = "<<lowerbound<<"  "<<upperbound<<endl;

      // Update scheme:
      if (fnew*f > T(0) && ABS(fnew) > T(0.1) * ABS(f) && k<N-1)
	fixedweight = !fixedweight;

      f = fnew;
#ifndef NWARN
      if (iter == MAXITER-1) {
	std::cout<<"Warning - Unable to find solution in FindDCEigenValue\n";
	std::cout<<"No solution after "<<MAXITER<<" iterations.\n";
      }
#endif
    }
    dbgcout<<"Found Eigenvalue S("<<k<<") = "<<s<<endl;

    return s;
  }

  template <class T> void FindDCEigenValues(
      Vector<T>& S, const T rho, const GenVector<T>& D, const GenVector<T>& z, 
      Matrix<T,ColMajor>& diffmat)
  {
    dbgcout<<"Start FindDCEV\n";
    dbgcout<<"D = "<<D<<endl;
    dbgcout<<"z = "<<z<<endl;
    TMVAssert(rho > T(0));
    TMVAssert(S.size() == D.size());
    TMVAssert(S.size() == z.size());
    TMVAssert(S.size() == diffmat.colsize());
    TMVAssert(S.size() == diffmat.rowsize());

    const int N = S.size();
    Vector<T> zsq(N);
    for(int j=0;j<N;j++) zsq[j] = z[j]*z[j];
    const T normsqz = zsq.SumElements();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int k=0;k<N;k++) {
      S[k] = FindDCEigenValue(k,N,rho,D.cptr(),z.cptr(),
	  zsq.cptr(),normsqz,diffmat.col(k).ptr());
    }
    dbgcout<<"S => "<<S<<endl;
  }

  template <class T> void FindDCEigenValues(
      Vector<T>& S, const T rho, const GenVector<T>& D, const GenVector<T>& z)
  {
    dbgcout<<"Start FindDCEV (No diffmat)\n";
    dbgcout<<"D = "<<D<<endl;
    dbgcout<<"z = "<<z<<endl;
    TMVAssert(rho > T(0));
    TMVAssert(S.size() == D.size());
    TMVAssert(S.size() == z.size());

    const int N = S.size();
    Vector<T> zsq(N);
    for(int j=0;j<N;j++) zsq[j] = z[j]*z[j];
    const T normsqz = zsq.SumElements();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      Vector<T> diff(N);
#ifdef _OPENMP
#pragma omp for
#endif
      for(int k=0;k<N;k++) {
	S[k] = FindDCEigenValue(k,N,rho,D.cptr(),z.cptr(),
	    zsq.cptr(),normsqz,diff.ptr());
      }
    }
    dbgcout<<"S => "<<S<<endl;
  }

  template <class T> static void SmallProblem(
      MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E)
  {
    if (E.size() > 0) Eigen_From_Tridiagonal_QR(U,D,E);
  }

  template <class T> void Eigen_From_Tridiagonal_DC(
      MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E, bool UisI)
  {
    // Solve the SVD of unreduced Tridiagonal Matrix T (given by D,E).
    // Note: the input T must be unreduced - ie. all off-diagonal entries 
    // are non-zero.
    // This routine implements the divide and conquer approach, calling
    // itself for the recursion, and the QR algorithm when the
    // size gets too small for divide and conquer to be efficient.
    //
    // The basic idea of the divide and conquer algorithm is to split
    // the tridiagonal matrix T into two parts, T1 and T2, plus a simple
    // rank-1 modification:
    // 
    // T = [ T1  0 ] + rho v ^ v
    //     [ 0  T2 ]
    //
    // v = [ 0 0 ... 0 a b 0 0 ... 0 ]
    // where rho * a * b = E_K,
    // T1 = T.SubMatrix(1,K,1,K) - rho*a*a e_K^e_K
    // T2 = T.SubMatrix(K,N,K,N) - rho*b*b e_1^e_1
    // 
    // The simplest version is rho = 1, a = 1, b = E_K.
    //
    // The smaller tridiagonal matrices are first reduced recursively:
    // T1 = U1 S1 U1t
    // T2 = U2 S2 U2t
    //
    //
    // T = [ U1  0 ] ( [ S1  0 ] + rho z ^ z ) [ U1  0 ]t
    //     [ 0  U2 ] ( [ 0  S2 ]             ) [ 0  U2 ]
    //
    // The vector z = [ z1 z2 ] is such that z Ut = [ 0 1 E_K 0 ]
    // which implies that:
    // z1 = U1.row(K)
    // z2 = E_K U2.row(1)
    //
    // Now, we are left with finding the eigen decomposition of a 
    // rank-1 modification of a diagonal matrix.
    // It turns out that this is not too difficult.
    //
    // Let the eigenvalues be Si, and the eigenvector for Si be
    // vi = (D - Si I)^-1 z
    //
    // Proof: (From Golub and Van Loan, 8.5.3)
    //
    // The definition of an eigenvector gives us:
    // (D + rho z zt) vi = Si vi
    // (D + rho z zt) vi - Si vi = 0
    // D vi + rho z (zt vi) - Si vi = 0
    // (D - Si I) vi + rho (zt vi) z = 0
    // vi = -rho (zt vi) (D - Si I)^-1 z
    // 
    // This assumes that (D - Si I) is non-singular, which is 
    // true so long as the elements of D are all different, and the 
    // elements of z are all non-zero. (Proof omitted.)
    //
    // So either vi = C (D - Si I)^-1 z (with C an arbitrary constant)
    // or (zt vi) = 0.
    //
    // zt vi leads to vi = 0, which is not a valid eigenvalue, so 
    // vi is a multiple of (D - Si I)^-1 z.  C is arbitrary, so set it = 1.
    //
    // (D - Si I) vi + rho (zt vi) z = 0
    // vit (D - Si I) vi + rho (zt vi) (vit z) = 0
    // zt vi + rho (zt vi) (vit z)  = 0
    // (zt vi) (1 + rho zt (D - Si I)^-1 z) = 0
    //
    // 1 + rho zt (D - Si I)^-1 z = 0
    //
    // The solutions to this equation are the eigenvalues Si.
    // And the corresponding eigenvectors are (D - Si I)^-1 z.
    //
    // The only thing left to worry about are how to deal with 
    // cases where the D(i) values are not distinct or where
    // some z_i = 0.
    //
    // Zero-valued z_i's are easy to deal with.  If z_k is zero (or
    // in practice less than EPS*(d_max+|z|) ) then d_k is an eigenvalue
    // itself, so you can permute it to the end and reduce the order
    // of D and z by one.
    //
    // If there are two d_i which are very close to each other (difference
    // < EPS*(d_max+|z|), we need to use a Givens rotation to rotate
    // z_j into z_i.  Doing this on the rows as well preserves the two d's,
    // but now we have a 0 in the z vector, and we can "deflate" it as
    // described above.
    //
    // Another concern, addressed by Gu and Eisenstat, is that it is 
    // numerically more stable to calculate the s values as above, 
    // but then to recalculate the z vector as:
    //
    // z_i' = sqrt( (sn-di) prod_1..i-1 (sk-di)/(dk-di) *
    //             prod_i..N-1 (sk-di)/(dk+1-di) ) * sign(z_i)
    //
    // Then v_i = (D - Si I)^-1 z'
    //
    // In exact arithmetic, it can be shown that z = z', but with 
    // rounding errors, the above value apparently gives better results.
    // Specifically, the resulting v_i are more nearly orthogonal.
    //

    TMVAssert(D.size()>0);
    TMVAssert(E.size()+1 == D.size());
    if (U) TMVAssert(U->rowsize() == D.size()); 
    TMVAssert(D.step()==1);
    TMVAssert(E.step()==1);

#ifdef XDEBUG
    dbgcout<<"Start SVD DC: \n";
    Vector<RT> D0 = D;
    dbgcout<<"D0 = "<<D0<<endl;
    Vector<RT> E0 = E;
    dbgcout<<"E0 = "<<E0<<endl;
    Matrix<T> U0(U?U->colsize():D.size(),D.size());
    if (U) U0 = *U; else U0.SetToIdentity();
    Matrix<T> A0(U0.colsize(),U0.colsize());
    Matrix<T> T0(D.size(),D.size(),RT(0));
    T0.diag() = D;
    if (E0.size() > 0) T0.diag(1) = T0.diag(-1) = E;
    A0 = U0 * T0 * U0.Adjoint();
    dbgcout<<"A0 = "<<A0<<endl;
    Matrix<T> U_QR = U0;
    Vector<RT> D_QR = D0;
    Vector<RT> E_QR = E0;
    SmallProblem<T>(U_QR.View(),D_QR.View(),E_QR.View());
#endif

    int N = D.size();
    // If N is too small, use the QR method
    if (N <= DC_LIMIT) SmallProblem(U,D,E);
    else {
      dbgcout<<"N > "<<DC_LIMIT<<endl;
      int K = (N-1)/2;
      Vector<RT> z(N,RT(0));

#ifdef _OPENMP
#pragma omp parallel sections 
#endif
      {
#ifdef _OPENMP
#pragma omp section
#endif
	{
	  // Do the left sub-problem
	  VectorView<RT> D1 = D.SubVector(0,K);
	  VectorView<RT> E1 = E.SubVector(0,K-1);
	  D1(K-1) -= RT(1);
	  Matrix<RT,ColMajor> U1(U?K:1,K);
	  if (U) {
	    U1.SetToIdentity();
	    Eigen_From_Tridiagonal_DC<RT>(U1.View(),D1,E1,true);
	    if (UisI) U->SubMatrix(0,K,0,K) = U1;
	    else U->Cols(0,K) *= U1;
	  } else {
	    U1.row(0).MakeBasis(K-1); // only need row(K-1)
	    Eigen_From_Tridiagonal_DC<RT>(U1.View(),D1,E1,false);
	  }
	  z.SubVector(0,K) = U1.row(U ? K-1 : 0);
	} // end left sub-problem

#ifdef _OPENMP
#pragma omp section
#endif
	{
	  // Do the right sub-problem
	  VectorView<RT> D2 = D.SubVector(K,N);
	  VectorView<RT> E2 = E.SubVector(K,N-1);
	  const RT EK = E(K-1);
	  D2(0) -= EK*EK;
	  Matrix<RT,ColMajor> U2(U?N-K:1,N-K);
	  if (U) {
	    U2.SetToIdentity();
	    Eigen_From_Tridiagonal_DC<RT>(U2.View(),D2,E2,true);
	    if (UisI) U->SubMatrix(K,N,K,N) = U2;
	    else U->Cols(K,N) *= U2;
	  } else {
	    U2.row(0).MakeBasis(0); // only need col(0)
	    Eigen_From_Tridiagonal_DC<RT>(U2.View(),D2,E2,false);
	  }
	  z.SubVector(K,N) = EK * U2.row(0);
	} // end right sub-problem
      } // end parallel sections
      dbgcout<<"Done subproblems\n";
      dbgcout<<"D = "<<D<<endl;
      dbgcout<<"z = "<<z<<endl;

#ifdef XDEBUG
      Matrix<RT> M = DiagMatrixViewOf(D) + (z^z);
      if (U) {
	dbgcout<<"M = "<<M<<endl;
	dbgcout<<"U M Ut = "<<*U*M*U->Adjoint()<<endl;
	dbgcout<<"Norm(UMUt-A0) = "<<Norm(*U*M*U->Adjoint()-A0)<<endl;
	if (Norm(*U*M*U->Adjoint()-A0) > 0.01*Norm(A0)) abort();
      }
#endif

      // Check for deflation step 1:
      const RT eps = Epsilon<T>();
      RT norminf = D.SubVector(0,N).NormInf();
      RT tol = eps*norminf;
      dbgcout<<"tol = "<<tol<<endl;
      for(int i=0;i<N-1;i++) {
	if (ABS(z(i)) < tol || z(i)*z(i)*eps == RT(0)) {
	  z.Swap(i,N-1);
	  D.Swap(i,N-1);
	  if (U) U->SwapCols(i,N-1);
#ifdef XDEBUG
	  M.SwapCols(i,N-1);
	  M.SwapRows(i,N-1);
#endif
	  --i; --N;
	  RT adn = ABS(D(N));
	  if (adn == norminf) {
	    norminf = D.SubVector(0,N).NormInf();
	    tol = eps*norminf;
	  }
	}
      }
      if (ABS(z(N-1)) < tol) --N;
      dbgcout<<"After deflation\n";
      dbgcout<<"DN = "<<D.SubVector(0,N)<<endl;
      dbgcout<<"zN = "<<z.SubVector(0,N)<<endl;
#ifdef XDEBUG
      if (U) {
	dbgcout<<"M = "<<M<<endl;
	dbgcout<<"U M Ut = "<<*U*M*U->Adjoint()<<endl;
	dbgcout<<"Norm(UMUt-A0) = "<<Norm(*U*M*U->Adjoint()-A0)<<endl;
	if (Norm(*U*M*U->Adjoint()-A0) > 0.01*Norm(A0)) abort();
      }
#endif

      if (N > 1) {
	// Sort the inner matrix so z is in the first row, D are increasing
	auto_array<int> P(new int[N]);
	D.SubVector(0,N).Sort(P.get(),ASCEND);
	z.SubVector(0,N).Permute(P.get());
	if (U) U->Cols(0,N).PermuteCols(P.get());
	dbgcout<<"After sort\n";
	dbgcout<<"DN = "<<D.SubVector(0,N)<<endl;
	dbgcout<<"zN = "<<z.SubVector(0,N)<<endl;
#ifdef XDEBUG
	if (U) {
	  M.Cols(0,N).PermuteCols(P.get());
	  M.Rows(0,N).PermuteRows(P.get());
	  dbgcout<<"M = "<<M<<endl;
	  dbgcout<<"U M Ut = "<<*U*M*U->Adjoint()<<endl;
	  dbgcout<<"Norm(UMUt-A0) = "<<Norm(*U*M*U->Adjoint()-A0)<<endl;
	  if (Norm(*U*M*U->Adjoint()-A0) > 0.01*Norm(A0)) abort();
	}
#endif

	// Check for deflation step 2:
	int i_firstswap = N;
	tol = eps*D(N-1);
	for(int i=N-1;i>0;--i) {
	  if (ABS(D(i) - D(i-1)) < tol) {
	    Givens<RT> G = Givens_Rotate(z(i-1),z(i));
	    TMVAssert(z(i) == RT(0));
	    if (U) G.Mult(U->ColPair(i-1,i).Transpose());
#ifdef XDEBUG
	    G.Mult(M.ColPair(i-1,i).Transpose());
	    G.Mult(M.RowPair(i-1,i));
#endif
	    if (i < N-1) {
	      z.Swap(i,N-1);
	      D.Swap(i,N-1);
	      if (U) U->SwapCols(i,N-1);
#ifdef XDEBUG
	      M.SwapCols(i,N-1);
	      M.SwapRows(i,N-1);
#endif
	      i_firstswap = i;
	    }
	    --N;
	  }
	}
	dbgcout<<"After second deflation\n";
	dbgcout<<"DN = "<<D.SubVector(0,N)<<endl;
	dbgcout<<"zN = "<<z.SubVector(0,N)<<endl;
#ifdef XDEBUG
	if (U) {
	  dbgcout<<"M = "<<M<<endl;
	  dbgcout<<"U M Ut = "<<*U*M*U->Adjoint()<<endl;
	  dbgcout<<"Norm(UMUt-A0) = "<<Norm(*U*M*U->Adjoint()-A0)<<endl;
	  if (Norm(*U*M*U->Adjoint()-A0) > 0.01*Norm(A0)) abort();
	}
#endif

	// Re-sort if we screwed up the sort in the second deflation.
	// The D(i) from i_firstswap on might be unsorted
	if (i_firstswap < N-1) {
	  D.SubVector(i_firstswap,N).Sort(P.get(),ASCEND);
	  z.SubVector(i_firstswap,N).Permute(P.get());
	  if (U) U->Cols(i_firstswap,N).PermuteCols(P.get());
	}
	dbgcout<<"After second sort\n";
	dbgcout<<"DN = "<<D.SubVector(0,N)<<endl;
	dbgcout<<"zN = "<<z.SubVector(0,N)<<endl;
#ifdef XDEBUG
	if (U) {
	  if (i_firstswap < N-1) {
	    M.Cols(i_firstswap,N).PermuteCols(P.get());
	    M.Rows(i_firstswap,N).PermuteRows(P.get());
	  }
	  dbgcout<<"M = "<<M<<endl;
	  dbgcout<<"U M Ut = "<<*U*M*U->Adjoint()<<endl;
	  dbgcout<<"Norm(UMUt-A0) = "<<Norm(*U*M*U->Adjoint()-A0)<<endl;
	  if (Norm(*U*M*U->Adjoint()-A0) > 0.01*Norm(A0)) abort();
	}
#endif
      }

      // After deflation, it is possible for N to be 1, in which case,
      // the SVD is done.
      if (N > 1) {
	VectorView<RT> DN = D.SubVector(0,N);
	VectorView<RT> zN = z.SubVector(0,N);

	// Find the new singular values. 
	Vector<RT> S(N);
	if (U) {
	  Matrix<RT,ColMajor> W(N,N); // W(i,j) = D(i)-S(j)
	  FindDCEigenValues(S,RT(1),DN,zN,W);
	  dbgcout<<"S = "<<S<<endl;
	  //dbgcout<<"D_QR = "<<D_QR<<endl;

	  // Update z's to numerically more stable values:
	  // z_i' = sqrt( (sn-di) prod_1..i-1 (sk-di)/(dk-di) *
	  //             prod_i..N-1 (sk-di)/(dk+1-di) ) * sign(z_i)
	  dbgcout<<"z = "<<z<<endl;
	  for(int i=0;i<N;i++) {
	    RT di = D(i);
	    RT prod = -W(i,N-1);
	    for(int k=0;k<i;k++) {
	      prod *= -W(i,k)/(D(k)-di);
	    }
	    for(int k=i+1;k<N;k++) {
	      prod *= -W(i,k-1)/(D(k)-di);
	    }
	    prod = SQRT(prod);
	    if (z(i) > 0) z(i) = prod;
	    else z(i) = -prod;
	  }
	  dbgcout<<"z => "<<z<<endl;

	  // Make X matrix
	  // x_j = ( z1/(d1-sj) , z2/(d2-sj) , ... zN/(dN-sj) )
	  // X.col(j) = x_j / Norm(x_j);

	  // Currently W(i,j) = di-sj
	  // Convert it to X
	  dbgcout<<"W = "<<W<<endl;
	  for(int j=0;j<N;j++) {
	    VectorView<RT> xj = W.col(j);
	    for(int i=0;i<N;i++) xj(i) = z(i) / xj(i);
	    xj /= Norm(xj);
	  }
	  dbgcout<<"W => (X) "<<W<<endl;
	  U->Cols(0,N) *= W;
#ifdef XDEBUG
	  if (U) {
	    M.SubMatrix(0,N,0,N) = W.Adjoint() * M.SubMatrix(0,N,0,N) * W;
	    dbgcout<<"M = "<<M<<endl;
	    dbgcout<<"U M Ut = "<<*U*M*U->Adjoint()<<endl;
	    dbgcout<<"Norm(UMUt-A0) = "<<Norm(*U*M*U->Adjoint()-A0)<<endl;
	    if (Norm(*U*M*U->Adjoint()-A0) > 0.01*Norm(A0)) abort();
	  }
#endif
	} else {
	  FindDCEigenValues(S,RT(1),DN,zN);
	  dbgcout<<"S = "<<S<<endl;
	}

	// Done.  Copy S to D
	D.SubVector(0,N) = S;
	dbgcout<<"Done D = "<<D<<endl;

#ifdef XDEBUG
	if (U) {
	  M.SubMatrix(0,N,0,N) = DiagMatrixViewOf(S);
	  dbgcout<<"M = "<<M<<endl;
	  dbgcout<<"U = "<<*U<<endl;
	  dbgcout<<"U M Ut = "<<*U*M*U->Adjoint()<<endl;
	  dbgcout<<"Norm(UMUt-A0) = "<<Norm(*U*M*U->Adjoint()-A0)<<endl;
	  if (Norm(*U*M*U->Adjoint()-A0) > 0.01*Norm(A0)) abort();
	}
#endif
      } else if (N == 1) {
	D(0) += z(0)*z(0);
#ifdef XDEBUG
	if (U) {
	  M(0,0) = D(0);
	  dbgcout<<"M = "<<M<<endl;
	  dbgcout<<"U M Ut = "<<*U*M*U->Adjoint()<<endl;
	  dbgcout<<"Norm(UMUt-A0) = "<<Norm(*U*M*U->Adjoint()-A0)<<endl;
	  if (Norm(*U*M*U->Adjoint()-A0) > 0.01*Norm(A0)) abort();
	}
#endif
      }
#ifdef XDEBUG
      if (U) {
	M = DiagMatrixViewOf(D);
	dbgcout<<"M = "<<M<<endl;
	dbgcout<<"U M Ut = "<<*U*M*U->Adjoint()<<endl;
	dbgcout<<"Norm(UMUt-A0) = "<<Norm(*U*M*U->Adjoint()-A0)<<endl;
	if (Norm(*U*M*U->Adjoint()-A0) > 0.01*Norm(A0)) abort();
      }
#endif
    }

#ifdef XDEBUG
    DiagMatrix<RT> S(D);
    Matrix<T> A2(A0.colsize(),A0.rowsize());
    if (U) A2 = *U * S * U->Adjoint();
    dbgcout<<"S = "<<S.diag()<<endl;
    dbgcout<<"QR solution: S = "<<D_QR<<endl;
    Vector<RT> D1 = D;
    D1.Sort(DESCEND,ABS_COMP);
    Vector<RT> D2 = D_QR;
    D2.Sort(DESCEND,ABS_COMP);
    if (U) dbgcout<<"Norm(A2-A0) = "<<Norm(A2-A0)<<endl;
    if (Norm(D1-D2) > 0.001*Norm(A0) ||
	(U && Norm(A2-A0) > 0.001*Norm(A0)) ) {
      cerr<<"Eigen_DC: \n";
      cerr<<"input D = "<<D0<<endl;
      cerr<<"input E = "<<E0<<endl;
      cerr<<"output S = "<<D<<endl;
      cerr<<"QR solution: S = "<<D_QR<<endl;
      if (U) {
	cerr<<"UBUt = "<<A0<<endl;
	cerr<<"USUt = "<<A2<<endl;
	cerr<<"U = "<<*U<<endl;
	cerr<<"S = "<<S<<endl;
      }
      abort();
    }
#endif
  }
   
#undef RT

#define InstFile "TMV_SymSVDecompose_DC.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv

