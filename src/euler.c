
/*------------------------------------------------------------------------------
 * FILE: euler.c
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 *------------------------------------------------------------------------------
 */

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
int prim_to_cons_array(const double *P, double *U, int N);
int cons_to_prim_array(const double *U, double *P, int N);

/*------------------------------------------------------------------------------
 *
 * Private Data
 *
 */
enum { rho, nrg, px, py, pz }; // Conserved
enum { RHO, pre, vx, vy, vz }; // Primitive

static double adiabatic_gamma=1.4;


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
      F[rho]  =  U[rho] * P[vx];
      F[nrg]  = (U[nrg] + P[pre])*P[vx];
      F[px]   =  U[px]  * P[vx] + P[pre];
      F[py]   =  U[py]  * P[vx];
      F[pz]   =  U[pz]  * P[vx];
      break;
    case 2:
      F[rho]  =  U[rho] * P[vy];
      F[nrg]  = (U[nrg] + P[pre])*P[vy];
      F[px]   =  U[px]  * P[vy];
      F[py]   =  U[py]  * P[vy] + P[pre];
      F[pz]   =  U[pz]  * P[vy];
      break;
    case 3:
      F[rho]  =  U[rho] * P[vz];
      F[nrg]  = (U[nrg] + P[pre])*P[vz];
      F[px]   =  U[px]  * P[vz];
      F[py]   =  U[py]  * P[vz];
      F[pz]   =  U[pz]  * P[vz] + P[pre];
      break;
    }

  switch (dim)
    {
    case 1:
      *ap = P[vx] + sqrt(adiabatic_gamma*P[pre] / U[rho]);
      *am = P[vx] - sqrt(adiabatic_gamma*P[pre] / U[rho]);
      break;
    case 2:
      *ap = P[vy] + sqrt(adiabatic_gamma*P[pre] / U[rho]);
      *am = P[vy] - sqrt(adiabatic_gamma*P[pre] / U[rho]);
      break;
    case 3:
      *ap = P[vz] + sqrt(adiabatic_gamma*P[pre] / U[rho]);
      *am = P[vz] - sqrt(adiabatic_gamma*P[pre] / U[rho]);
      break;
    }

  return 0;
}


int cons_to_prim_point(const double *U, double *P)
{
  const double gm1 = adiabatic_gamma-1.0;
  P[rho] = U[rho];
  P[pre] =(U[nrg] - 0.5*(U[px]*U[px] + U[py]*U[py] + U[pz]*U[pz])/U[rho])*gm1;
  P[vx ] = U[px ] / U[rho];
  P[vy ] = U[py ] / U[rho];
  P[vz ] = U[pz ] / U[rho];
  return 0;
}
int prim_to_cons_point(const double *P, double *U)
{
  const double gm1 = adiabatic_gamma-1.0;
  U[rho] = P[rho];
  U[px]  = P[rho] * P[vx];
  U[py]  = P[rho] * P[vy];
  U[pz]  = P[rho] * P[vz];
  U[nrg] = P[rho] * 0.5*(P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz]) + P[pre]/gm1;
  return 0;
}

int cons_to_prim_array(const double *U, double *P, int N)
{
  int i;
  for (i=0; i<N*5; i+=5) cons_to_prim_point(&U[i],&P[i]);
  return 0;
}
int prim_to_cons_array(const double *P, double *U, int N)
{
  int i;
  for (i=0; i<N*5; i+=5) prim_to_cons_point(&P[i],&U[i]);
  return 0;
}


int constrained_transport_2d(double *Fx, double *Fy)
{
  return 0;
}
int constrained_transport_3d(double *Fx, double *Fy, double *Fz)
{
  return 0;
}
