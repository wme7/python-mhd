
/*------------------------------------------------------------------------------
 * FILE: hllc_rmhd.c
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * PURPOSE: Obtain the intercell flux between two constant states using a 3-wave
 *   approximation.
 *
 * REFERENCES: Mignone & Bodo (2006) An HLLC Solver for Relativistic Flows
 *
 *------------------------------------------------------------------------------
 */

#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include <math.h>

#define SMALL_BX 1e-10

// Enums used for indexing through primitive/conserved
// data. The state of a given cell is described by 8
// contiguous doubles.

enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive

int rmhd_flux_and_eval(const double *U, const double *P, double *F, double *ap, double *am);
int prim_to_cons_point(const double *P, double *U);
int cons_to_prim_point(const double *U, double *P);

int hllc_flux(const double *pl, const double *pr, double *U, double *F, double s)
{
  int i;
  double epl, epr, eml, emr;
  double Ul[8], Ur[8];
  double Pl[8], Pr[8];
  double Fl[8], Fr[8];

  memcpy(Pl,pl,8*sizeof(double));
  memcpy(Pr,pr,8*sizeof(double));

  prim_to_cons_point(Pl,Ul);
  prim_to_cons_point(Pr,Ur);

  rmhd_flux_and_eval(Ul, Pl, Fl, &epl, &eml);
  rmhd_flux_and_eval(Ur, Pr, Fr, &epr, &emr);

  double ap = (epl>epr) ? epl : epr;
  double am = (eml<emr) ? eml : emr;

  Ul[tau] += Ul[ddd];  Fl[tau] += Fl[ddd]; // Change in convention of total energy
  Ur[tau] += Ur[ddd];  Fr[tau] += Fr[ddd];

  Ul[Bx] = Ur[Bx] = 0.5*(Ul[Bx] + Ur[Bx]); // Must have no normal jump in B
  Fl[Bx] = Fr[Bx] = 0.5*(Fl[Bx] + Fr[Bx]);

  double F_hll[8], U_hll[8];
  for (i=0; i<8; ++i)
    {
      U_hll[i] = (ap*Ur[i] - am*Ul[i] +       (Fl[i] - Fr[i])) / (ap - am);
      F_hll[i] = (ap*Fl[i] - am*Fr[i] + ap*am*(Ur[i] - Ul[i])) / (ap - am);
    }

  double Ul_[8], Ur_[8], P_[8], lc; // The star states
  if (fabs(U_hll[Bx]) > SMALL_BX)
    {
      const double B_dot_FB = U_hll[By]*F_hll[By] + U_hll[Bz]*F_hll[Bz];
      const double a =  F_hll[tau] - B_dot_FB; // eqns (42)
      const double b = -F_hll[Sx ] - U_hll[tau] +
	(U_hll[By]*U_hll[By] + U_hll[Bz]*U_hll[Bz]) +
	(F_hll[By]*F_hll[By] + F_hll[Bz]*F_hll[Bz]);
      const double c =  U_hll[Sx ] - B_dot_FB;

      P_[vx] = lc = (-b - sqrt(b*b - 4*a*c)) / (2*a); // take root with the minus sign
      P_[vy] = (U_hll[By]*P_[vx] - F_hll[By]) / U_hll[Bx]; // eqn (38)
      P_[vz] = (U_hll[Bz]*P_[vx] - F_hll[Bz]) / U_hll[Bx];

      P_[Bx] = U_hll[Bx];
      P_[By] = U_hll[By];
      P_[Bz] = U_hll[Bz];

      const double v_dotB_ = P_[vx]*P_[Bx] + P_[vy]*P_[By] + P_[vz]*P_[Bz];
      const double v_dotv_ = P_[vx]*P_[vx] + P_[vy]*P_[vy] + P_[vz]*P_[vz];
      const double gm2_ = 1.0 / (1.0 - v_dotv_);

      P_[pre] = F_hll[Sx] + P_[Bx]*P_[Bx]/gm2_ - (F_hll[tau] - P_[Bx]*v_dotB_)*P_[vx];

      Ul_[ddd] = (am - Pl[vx]) / (am - P_[vx]) * Ul[ddd];
      Ur_[ddd] = (ap - Pr[vx]) / (ap - P_[vx]) * Ur[ddd];

      Ul_[tau] = (am*Ul[tau] - Ul[Sx] + P_[pre]*P_[vx] - v_dotB_*P_[Bx]) / (am - P_[vx]);
      Ur_[tau] = (ap*Ur[tau] - Ur[Sx] + P_[pre]*P_[vx] - v_dotB_*P_[Bx]) / (ap - P_[vx]);

      Ul_[Sy ] = (-P_[Bx]*(P_[By]/gm2_ + v_dotB_*P_[vy]) + am*Ul[Sy] - Fl[Sy]) / (am - P_[vx]);
      Ur_[Sy ] = (-P_[Bx]*(P_[By]/gm2_ + v_dotB_*P_[vy]) + ap*Ur[Sy] - Fr[Sy]) / (ap - P_[vx]);

      Ul_[Sz ] = (-P_[Bx]*(P_[Bz]/gm2_ + v_dotB_*P_[vz]) + am*Ul[Sz] - Fl[Sz]) / (am - P_[vx]);
      Ur_[Sz ] = (-P_[Bx]*(P_[Bz]/gm2_ + v_dotB_*P_[vz]) + ap*Ur[Sz] - Fr[Sz]) / (ap - P_[vx]);

      Ul_[Sx ] = (Ul_[tau] + P_[pre])*P_[vx] - v_dotB_*P_[Bx];
      Ur_[Sx ] = (Ur_[tau] + P_[pre])*P_[vx] - v_dotB_*P_[Bx];

      Ul_[Bx] = Ur_[Bx] = P_[Bx];
      Ul_[By] = Ur_[By] = P_[By];
      Ul_[Bz] = Ur_[Bz] = P_[Bz];
    }
  else
    {
      const double a =  F_hll[tau];
      const double b = -F_hll[Sx ] - U_hll[tau];
      const double c =  U_hll[Sx ];

      const double vx_ = lc = (-b - sqrt(b*b - 4*a*c)) / (2*a);
      const double p_  = -F_hll[tau]*vx_ + F_hll[Sx];

      Ul_[ddd] = (am - Pl[vx]) / (am - vx_) * Ul[ddd];
      Ur_[ddd] = (ap - Pr[vx]) / (ap - vx_) * Ur[ddd];

      Ul_[tau] = (am*Ul[tau] - Ul[Sx] + p_*vx_) / (am - vx_);
      Ur_[tau] = (ap*Ur[tau] - Ur[Sx] + p_*vx_) / (ap - vx_);

      Ul_[Sx ] = (Ul_[tau] + p_)*vx_;
      Ur_[Sx ] = (Ur_[tau] + p_)*vx_;

      Ul_[Sy ] = (am - Pl[vx]) / (am - vx_) * Ul[Sy];
      Ur_[Sy ] = (ap - Pr[vx]) / (ap - vx_) * Ur[Sy];

      Ul_[Sz ] = (am - Pl[vx]) / (am - vx_) * Ul[Sz];
      Ur_[Sz ] = (ap - Pr[vx]) / (ap - vx_) * Ur[Sz];

      Ul_[Bx ] = U_hll[Bx];
      Ur_[Bx ] = U_hll[Bx];

      Ul_[By ] = (am - Pl[vx]) / (am - vx_) * Ul[By];
      Ur_[By ] = (ap - Pr[vx]) / (ap - vx_) * Ur[By];

      Ul_[Bz ] = (am - Pl[vx]) / (am - vx_) * Ul[Bz];
      Ur_[Bz ] = (ap - Pr[vx]) / (ap - vx_) * Ur[Bz];
    }

  if      (         s<=am ) for (i=0; i<8; ++i) U[i] = Ul [i];
  else if ( am<s && s<=lc ) for (i=0; i<8; ++i) U[i] = Ul_[i];
  else if ( lc<s && s<=ap ) for (i=0; i<8; ++i) U[i] = Ur_[i];
  else if ( ap<s          ) for (i=0; i<8; ++i) U[i] = Ur [i];

  if      (         s<=am ) for (i=0; i<8; ++i) F[i] = Fl[i];
  else if ( am<s && s<=lc ) for (i=0; i<8; ++i) F[i] = Fl[i] + am*(Ul_[i]-Ul[i]);
  else if ( lc<s && s<=ap ) for (i=0; i<8; ++i) F[i] = Fr[i] + ap*(Ur_[i]-Ur[i]);
  else if ( ap<s          ) for (i=0; i<8; ++i) F[i] = Fr[i];

  U[tau] -= U[ddd]; // Change in convention of total energy
  F[tau] -= F[ddd];

  return 0;
}
