
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

#define SMALL_BX 1e-12
#define SMALL_A  1e-10

// Enums used for indexing through primitive/conserved
// data. The state of a given cell is described by 8
// contiguous doubles.

enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive

static int dimension = 1;

int rmhd_flux_and_eval(const double *U, const double *P, double *F, double *ap, double *am);
int prim_to_cons_point(const double *P, double *U);
int cons_to_prim_point(const double *U, double *P);

int hllc_set_dimension(int d) { dimension = d; }
int hllc_flux(const double *pl, const double *pr, double *U, double *F, double s)
{
  int i;
  double epl, epr, eml, emr;
  double Ul[8], Ur[8];
  double Pl[8], Pr[8];
  double Fl[8], Fr[8];

  int B1,B2,B3;
  int S1,S2,S3;
  int v1,v2,v3;

  switch (dimension)
    {
    case 1:
      B1=Bx; B2=By; B3=Bz;
      S1=Sx; S2=Sy; S3=Sz;
      v1=vx; v2=vy; v3=vz;
      break;

    case 2:
      B1=By; B2=Bz; B3=Bx;
      S1=Sy; S2=Sz; S3=Sx;
      v1=vy; v2=vz; v3=vx;
      break;

    case 3:
      B1=Bz; B2=Bx; B3=By;
      S1=Sz; S2=Sx; S3=Sy;
      v1=vz; v2=vx; v3=vy;
      break;
    }

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

  Ul[B1] = Ur[B1] = 0.5*(Ul[B1] + Ur[B1]); // Must have no normal jump in B
  Fl[B1] = Fr[B1] = 0.5*(Fl[B1] + Fr[B1]);

  double F_hll[8], U_hll[8];
  for (i=0; i<8; ++i)
    {
      U_hll[i] = (ap*Ur[i] - am*Ul[i] +       (Fl[i] - Fr[i])) / (ap - am);
      F_hll[i] = (ap*Fl[i] - am*Fr[i] + ap*am*(Ur[i] - Ul[i])) / (ap - am);
    }

  double Ul_[8], Ur_[8], P_[8], lc; // The star states
  if (fabs(U_hll[B1]) > SMALL_BX)
    {
      const double B_dot_FB = U_hll[B2]*F_hll[B2] + U_hll[B3]*F_hll[B3];
      const double a =  F_hll[tau] - B_dot_FB; // eqns (42)
      const double b = -F_hll[S1 ] - U_hll[tau] +
	(U_hll[B2]*U_hll[B2] + U_hll[B3]*U_hll[B3]) +
	(F_hll[B2]*F_hll[B2] + F_hll[B3]*F_hll[B3]);
      const double c =  U_hll[S1 ] - B_dot_FB;

      P_[v1] = lc = (fabs(a) < SMALL_A) ? -c/b : (-b - sqrt(b*b - 4*a*c)) / (2*a);
      P_[v2] = (U_hll[B2]*P_[v1] - F_hll[B2]) / U_hll[B1]; // eqn (38)
      P_[v3] = (U_hll[B3]*P_[v1] - F_hll[B3]) / U_hll[B1];

      P_[B1] = U_hll[B1];
      P_[B2] = U_hll[B2];
      P_[B3] = U_hll[B3];

      const double v_dotB_ = P_[v1]*P_[B1] + P_[v2]*P_[B2] + P_[v3]*P_[B3];
      const double v_dotv_ = P_[v1]*P_[v1] + P_[v2]*P_[v2] + P_[v3]*P_[v3];
      const double gm2_ = 1.0 / (1.0 - v_dotv_);

      P_[pre] = F_hll[S1] + P_[B1]*P_[B1]/gm2_ - (F_hll[tau] - P_[B1]*v_dotB_)*P_[v1];

      Ul_[ddd] = (am - Pl[v1]) / (am - P_[v1]) * Ul[ddd];
      Ur_[ddd] = (ap - Pr[v1]) / (ap - P_[v1]) * Ur[ddd];

      Ul_[tau] = (am*Ul[tau] - Ul[S1] + P_[pre]*P_[v1] - v_dotB_*P_[B1]) / (am - P_[v1]);
      Ur_[tau] = (ap*Ur[tau] - Ur[S1] + P_[pre]*P_[v1] - v_dotB_*P_[B1]) / (ap - P_[v1]);

      Ul_[S1 ] = (Ul_[tau] + P_[pre])*P_[v1] - v_dotB_*P_[B1];
      Ur_[S1 ] = (Ur_[tau] + P_[pre])*P_[v1] - v_dotB_*P_[B1];

      Ul_[S2 ] = (-P_[B1]*(P_[B2]/gm2_ + v_dotB_*P_[v2]) + am*Ul[S2] - Fl[S2]) / (am - P_[v1]);
      Ur_[S2 ] = (-P_[B1]*(P_[B2]/gm2_ + v_dotB_*P_[v2]) + ap*Ur[S2] - Fr[S2]) / (ap - P_[v1]);

      Ul_[S3 ] = (-P_[B1]*(P_[B3]/gm2_ + v_dotB_*P_[v3]) + am*Ul[S3] - Fl[S3]) / (am - P_[v1]);
      Ur_[S3 ] = (-P_[B1]*(P_[B3]/gm2_ + v_dotB_*P_[v3]) + ap*Ur[S3] - Fr[S3]) / (ap - P_[v1]);

      Ul_[B1] = Ur_[B1] = P_[B1];
      Ul_[B2] = Ur_[B2] = P_[B2];
      Ul_[B3] = Ur_[B3] = P_[B3];
    }
  else
    {
      const double a =  F_hll[tau];
      const double b = -F_hll[S1 ] - U_hll[tau];
      const double c =  U_hll[S1 ];

      const double v1_ = lc = (fabs(a) < SMALL_A) ? -c/b : (-b - sqrt(b*b - 4*a*c)) / (2*a);
      const double p_  = -F_hll[tau]*v1_ + F_hll[S1];

      Ul_[ddd] = (am - Pl[v1]) / (am - v1_) * Ul[ddd];
      Ur_[ddd] = (ap - Pr[v1]) / (ap - v1_) * Ur[ddd];

      Ul_[tau] = (am*Ul[tau] - Ul[S1] + p_*v1_) / (am - v1_);
      Ur_[tau] = (ap*Ur[tau] - Ur[S1] + p_*v1_) / (ap - v1_);

      Ul_[S1 ] = (Ul_[tau] + p_)*v1_;
      Ur_[S1 ] = (Ur_[tau] + p_)*v1_;

      Ul_[S2 ] = (am - Pl[v1]) / (am - v1_) * Ul[S2];
      Ur_[S2 ] = (ap - Pr[v1]) / (ap - v1_) * Ur[S2];

      Ul_[S3 ] = (am - Pl[v1]) / (am - v1_) * Ul[S3];
      Ur_[S3 ] = (ap - Pr[v1]) / (ap - v1_) * Ur[S3];

      Ul_[B1 ] = U_hll[B1];
      Ur_[B1 ] = U_hll[B1];

      Ul_[B2 ] = (am - Pl[v1]) / (am - v1_) * Ul[B2];
      Ur_[B2 ] = (ap - Pr[v1]) / (ap - v1_) * Ur[B2];

      Ul_[B3 ] = (am - Pl[v1]) / (am - v1_) * Ul[B3];
      Ur_[B3 ] = (ap - Pr[v1]) / (ap - v1_) * Ur[B3];
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
