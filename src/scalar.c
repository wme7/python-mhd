
/*------------------------------------------------------------------------------
 * FILE: scalar.c
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
void set_flux_vector(double *A);
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
static double Ax = 1.0; // Wavespeeds
static double Ay = 0.0;
static double Az = 0.0;

/*------------------------------------------------------------------------------
 *
 * Public Functions Definitions
 *
 */
void set_flux_vector(double *A)
{
  Ax = A[0];
  Ay = A[1];
  Az = A[2];
}
int flux_and_eval(const double *U, const double *P, double *F,
		  double *ap, double *am, int dim)
{
  switch (dim)
    {
    case 1: F[0] = Ax*U[0]; break;
    case 2: F[0] = Ay*U[0]; break;
    case 3: F[0] = Az*U[0]; break;
    }

  if (ap == 0 || am == 0) return 0; // User may skip eigenvalue calculation

  switch (dim)
    {
    case 1:
      *ap = (Ax>0) ? Ax : 0;
      *am = (Ax<0) ? Ax : 0;
      break;

    case 2:
      *ap = (Ay>0) ? Ay : 0;
      *am = (Ay<0) ? Ay : 0;
      break;

    case 3:
      *ap = (Ay>0) ? Az : 0;
      *am = (Ay<0) ? Az : 0;
      break;
    }

  return 0;
}

int cons_to_prim_point(const double *U, double *P)
{
  P[0] = U[0];
  return 0;
}
int prim_to_cons_point(const double *P, double *U)
{
  U[0] = P[0];
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
