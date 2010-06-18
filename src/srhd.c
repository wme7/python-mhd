
/*------------------------------------------------------------------------------
 * FILE: srhd.c
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 *------------------------------------------------------------------------------
 */

#include <stdlib.h>
#include <memory.h>
#include <math.h>

/*------------------------------------------------------------------------------
 *
 * Public Interface
 *
 */
void set_adiabatic_gamma(double g);
int flux_and_eval(const double *U, const double *P, double *F,
		  double *ap, double *am, int dim);
int prim_to_cons_point(const double *P, double *U);
int cons_to_prim_point(const double *U, double *P);

int constrained_transport_2d(double *Fx, double *Fy, int stride[4]);
int constrained_transport_3d(double *Fx, double *Fy, double *Fz, int stride[4]);

/*------------------------------------------------------------------------------
 *
 * Private Data
 *
 */
enum { ddd, tau, Sx, Sy, Sz }; // Conserved
enum { rho, pre, vx, vy, vz }; // Primitive

static double adiabatic_gamma=1.4;
static double FailedPrimState[8];
static double FailedConsState[8];

/*------------------------------------------------------------------------------
 *
 * Inline functions
 *
 */
static inline double eos_pre(double Rho, double Sie) /* Adiabatic equation of state */
{
  return Sie * (Rho * (adiabatic_gamma - 1.0));
}
static inline double eos_sie(double Rho, double Pre)
{
  return Pre / (Rho * (adiabatic_gamma - 1.0));
}
static inline double eos_cs2(double Rho, double Pre)
{
  double e = eos_sie(Rho, Pre);
  return adiabatic_gamma * Pre / (Pre + Rho + Rho*e);
}

/*------------------------------------------------------------------------------
 *
 * Public Functions Definitions
 *
 */
void set_adiabatic_gamma(double g)
{
  adiabatic_gamma = g;
}

int flux_and_eval(const double *U, const double *P, double *F,
		  double *ap, double *am, int dim)
{
  switch (dim)
    {
    case 1:
      F[ddd] = U[ddd] * P[vx ];
      F[tau] = U[Sx ] - U[ddd] * P[vx ];
      F[Sx ] = U[Sx ] * P[vx ] + P[pre];
      F[Sy ] = U[Sy ] * P[vx ];
      F[Sz ] = U[Sz ] * P[vx ];
      break;
    case 2:
      F[ddd] = U[ddd] * P[vy ];
      F[tau] = U[Sy ] - U[ddd] * P[vy ];
      F[Sx ] = U[Sx ] * P[vy ];
      F[Sy ] = U[Sy ] * P[vy ] + P[pre];
      F[Sz ] = U[Sz ] * P[vy ];
      break;
    case 3:
      F[ddd] = U[ddd] * P[vz ];
      F[tau] = U[Sz ] - U[ddd] * P[vz ];
      F[Sx ] = U[Sx ] * P[vz ];
      F[Sy ] = U[Sy ] * P[vz ];
      F[Sz ] = U[Sz ] * P[vz ] + P[pre];
      break;
    }

  if (ap == 0 || am == 0) return 0; // User may skip eigenvalue calculation

  const double vx2 = P[vx]*P[vx];
  const double vy2 = P[vy]*P[vy];
  const double vz2 = P[vz]*P[vz];
  const double cs2 = eos_cs2(P[rho], P[pre]);
  const double v2 = vx2 + vy2 + vz2;

  switch (dim)
    {
    case 1:
      *ap = (P[vx]*(1-cs2) + sqrt(cs2*(1-v2)*(1-v2*cs2-vx2*(1-cs2))))/(1-v2*cs2);
      *am = (P[vx]*(1-cs2) - sqrt(cs2*(1-v2)*(1-v2*cs2-vx2*(1-cs2))))/(1-v2*cs2);
      break;
    case 2:
      *ap = (P[vy]*(1-cs2) + sqrt(cs2*(1-v2)*(1-v2*cs2-vy2*(1-cs2))))/(1-v2*cs2);
      *am = (P[vy]*(1-cs2) - sqrt(cs2*(1-v2)*(1-v2*cs2-vy2*(1-cs2))))/(1-v2*cs2);
      break;
    case 3:
      *ap = (P[vz]*(1-cs2) + sqrt(cs2*(1-v2)*(1-v2*cs2-vz2*(1-cs2))))/(1-v2*cs2);
      *am = (P[vz]*(1-cs2) - sqrt(cs2*(1-v2)*(1-v2*cs2-vz2*(1-cs2))))/(1-v2*cs2);
      break;
    }

  return 0;
}

int cons_to_prim_point(const double *U, double *P)
{
  static const double ERROR_TOLR = 1e-6;
  static const int NEWTON_MAX_ITER = 25;

  const double D     = U[ddd];
  const double Tau   = U[tau];
  const double S2    = U[Sx]*U[Sx] + U[Sy]*U[Sy] + U[Sz]*U[Sz];

  int soln_found     = 0;
  int n_iter         = 0;

  double f,g,W_soln=1.0,p=P[pre];
  while (!soln_found)
    {
      const double v2  = S2 / pow(Tau + D + p, 2);
      const double W2  = 1.0 / (1.0 - v2);
      const double W   = sqrt(W2);
      const double e   = (Tau + D*(1.0 - W) + p*(1.0 - W2)) / (D*W);
      const double Rho = D / W;
      const double h   = 1.0 + e + p/Rho;
      const double cs2 = (p/Rho + e)*(adiabatic_gamma - 1.0) / h;

      f = eos_pre(Rho, e) - p;
      g = v2*cs2 - 1.0;

      p -= f/g;

      if (fabs(f) < ERROR_TOLR)
	{
	  W_soln = W;
	  soln_found = 1;
	}
      if (n_iter++ == NEWTON_MAX_ITER)
	{
	  memcpy(FailedConsState, U, 5*sizeof(double));
	  memcpy(FailedPrimState, P, 5*sizeof(double));
	  return 1;
	}
    }

  P[rho] = D/W_soln;
  P[pre] = p;
  P[vx]  = U[Sx] / (Tau + D + p);
  P[vy]  = U[Sy] / (Tau + D + p);
  P[vz]  = U[Sz] / (Tau + D + p);

  return 0;
}
int prim_to_cons_point(const double *P, double *U)
{
  const double v2   =   P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
  const double W2   =   1.0 / (1.0 - v2);
  const double W    =   sqrt(W2);
  const double e    =   eos_sie(P[rho], P[pre]);
  const double h    =   1.0 + e + P[pre]/P[rho];

  U[ddd] = P[rho]*W;
  U[tau] = P[rho]*h*W2 - P[pre] - U[ddd];
  U[Sx]  = P[rho]*h*W2*P[vx];
  U[Sy]  = P[rho]*h*W2*P[vy];
  U[Sz]  = P[rho]*h*W2*P[vz];

  return 0;
}

int constrained_transport_2d(double *Fx, double *Fy, int stride[4])
{
  return 0;
}
int constrained_transport_3d(double *Fx, double *Fy, double *Fz, int stride[4])
{
  return 0;
}
