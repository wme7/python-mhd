
/*------------------------------------------------------------------------------
 * FILE: rmhd.c
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
enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive

static double adiabatic_gamma=1.4;
static double FailedPrimState[8];
static double FailedConsState[8];

/*------------------------------------------------------------------------------
 *
 * External Dependecies
 *
 */
int quartic_equation_init(double d4, double d3, double d2, double d1, double d0);
int quartic_equation_solve_apprx(double *x);
int quartic_equation_solve_exact(double *r1, double *r2, double *r3, double *r4,
				 int *nr12, int *nr34);

/*------------------------------------------------------------------------------
 *
 * Inline functions
 *
 */
static inline void invert_2by2_matrix(const double A[2][2], double B[2][2])
{
  double det = A[0][0]*A[1][1] - A[1][0]*A[0][1];

  B[0][0] =  A[1][1] / det;
  B[1][1] =  A[0][0] / det;
  B[0][1] = -A[0][1] / det;
  B[1][0] = -A[1][0] / det;
}
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
  const double v2   =   P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
  const double B2   =   P[Bx]*P[Bx] + P[By]*P[By] + P[Bz]*P[Bz];
  const double Bv   =   P[Bx]*P[vx] + P[By]*P[vy] + P[Bz]*P[vz];
  const double W2   =   1.0 / (1.0 - v2);
  const double W    =   sqrt(W2);
  const double b0   =   W*Bv;
  const double b2   =   (B2 + b0*b0) / W2;
  const double bx   =   (P[Bx] + b0 * W*P[vx]) / W;
  const double by   =   (P[By] + b0 * W*P[vy]) / W;
  const double bz   =   (P[Bz] + b0 * W*P[vz]) / W;
  const double e    =   eos_sie(P[rho], P[pre]);
  const double p    =   P[pre];
  const double h    =   1.0 + e + P[pre]/P[rho];
  const double p_   =   p + 0.5 * b2;

  switch (dim)
    {
    case 1:
      F[ddd] = U[ddd] * P[vx];
      F[tau] = U[tau] * P[vx] - b0*P[Bx] / W + p_*P[vx];
      F[Sx ] = U[Sx]  * P[vx] - bx*P[Bx] / W + p_;
      F[Sy ] = U[Sy]  * P[vx] - by*P[Bx] / W;
      F[Sz ] = U[Sz]  * P[vx] - bz*P[Bx] / W;
      F[Bx ] = 0.0;
      F[By ] = P[vx]*P[By] - P[vy]*P[Bx];
      F[Bz ] = P[vx]*P[Bz] - P[vz]*P[Bx];
      break;
    case 2:
      F[ddd] = U[ddd] * P[vy];
      F[tau] = U[tau] * P[vy] - b0*P[By] / W + p_*P[vy];
      F[Sx ] = U[Sx]  * P[vy] - bx*P[By] / W;
      F[Sy ] = U[Sy]  * P[vy] - by*P[By] / W + p_;
      F[Sz ] = U[Sz]  * P[vy] - bz*P[By] / W;
      F[Bx ] = P[vy]*P[Bx] - P[vx]*P[By];
      F[By ] = 0.0;
      F[Bz ] = P[vy]*P[Bz] - P[vz]*P[By];
      break;
    case 3:
      F[ddd] = U[ddd] * P[vz];
      F[tau] = U[tau] * P[vz] - b0*P[Bz] / W + p_*P[vz];
      F[Sx ] = U[Sx]  * P[vz] - bx*P[Bz] / W;
      F[Sy ] = U[Sy]  * P[vz] - by*P[Bz] / W;
      F[Sz ] = U[Sz]  * P[vz] - bz*P[Bz] / W + p_;
      F[Bx ] = P[vz]*P[Bx] - P[vx]*P[Bz];
      F[By ] = P[vz]*P[By] - P[vy]*P[Bz];
      F[Bz ] = 0.0;
      break;
    }

  double vi=0, bi=0;
  switch (dim)
    {
    case 1:
      vi = P[vx]; bi = bx;
      break;
    case 2:
      vi = P[vy]; bi = by;
      break;
    case 3:
      vi = P[vz]; bi = bz;
      break;
    }

  if (ap == 0 || am == 0) return 0; // User may skip eigenvalue calculation

  const double W4   =  W2*W2;
  const double cs2  =  eos_cs2(P[rho],P[pre]);
  const double V2   =  vi*vi;
  const double V3   =  vi*V2;
  const double V4   =  vi*V3;

  const double K  =    P[rho]*h * (1.0/cs2-1.0) * W4;
  const double L  =  -(P[rho]*h +   b2/cs2)     * W2;

  const double A4 =    K    - L          -   b0*b0;
  const double A3 = -4*K*vi + L*vi*2     + 2*b0*bi;
  const double A2 =  6*K*V2 + L*(1.0-V2) +   b0*b0 - bi*bi;
  const double A1 = -4*K*V3 - L*vi*2     - 2*b0*bi;
  const double A0 =    K*V4 + L*V2       +   bi*bi;

  quartic_equation_init(A4,A3,A2,A1,A0);

  double r1, r2, r3, r4;
  int nr12, nr34;
  int nr = quartic_equation_solve_exact(&r1, &r2, &r3, &r4, &nr12, &nr34);

  double ap12 = (r1>r2) ? r1 : r2;
  double ap34 = (r3>r4) ? r3 : r4;

  double am12 = (r1<r2) ? r1 : r2;
  double am34 = (r3<r4) ? r3 : r4;

  *ap = (nr==2) ? ((nr12==2) ? ap12 : ap34) : ((ap12>ap34) ? ap12 : ap34);
  *am = (nr==2) ? ((nr12==2) ? am12 : am34) : ((am12<am34) ? am12 : am34);

  if (fabs(*ap)>1.0 || fabs(*am)>1.0)
    {
      *ap =  1.0;
      *am = -1.0;
    }

  return 0;
}

int cons_to_prim_point(const double *U, double *P)
{
  static const double PRES_FLOOR = 1e-10;
  static const double ERROR_TOLR = 1e-6;
  static const int NEWTON_MAX_ITER = 25;

  // Quantites known from conserved
  const double gamf  = (adiabatic_gamma - 1.0) / adiabatic_gamma;
  const double D     = U[ddd];
  const double Tau   = U[tau];
  const double S2    = U[Sx]*U[Sx] + U[Sy]*U[Sy] + U[Sz]*U[Sz];
  const double B2    = U[Bx]*U[Bx] + U[By]*U[By] + U[Bz]*U[Bz];
  const double BS    = U[Bx]*U[Sx] + U[By]*U[Sy] + U[Bz]*U[Sz];
  const double BS2   = BS*BS;

  int est            = 0;
  int use_pres_floor = 0;
  int soln_found     = 0;
  int n_iter         = 0;

  // Quantites guessed from existing primitives, or estimated from conserved
  double v2 = P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
  double h_guess = 1.0 + eos_sie(P[rho], P[pre]) + P[pre]/P[rho];
  double W_guess = 1.0 / sqrt(1.0 - v2);
  double W = (est) ? sqrt(S2/(D*D) + 1) : W_guess;
  double Z = (est) ? D*W                : P[rho]*h_guess*W_guess*W_guess;

  const double bigZ = 1e20;
  const double bigW = 1e12;
  const double smlZ = 0.0;
  const double smlW = 1.0;

  while (!soln_found)
    {
      const double Z2 = Z*Z;
      const double Z3 = Z*Z2;
      const double W2 = W*W;
      const double W3 = W*W2;

      const double Pre = (use_pres_floor) ? PRES_FLOOR : (D/W) * (Z/(D*W) - 1.0) * gamf;

      const double f1 = -S2  +  (Z+B2)*(Z+B2)*(W2-1)/W2 - (2*Z+B2)*BS2/Z2;      // eqn (84)
      const double f2 = -Tau +   Z+B2 - Pre - 0.5*B2/W2 -      0.5*BS2/Z2 - D;  // eqn (85)

      const double df1dZ = 2*(B2+Z)*(BS2*W2 + (W2-1)*Z3) / (W2*Z3);
      const double df1dW = 2*(B2+Z)*(B2+Z) / W3;
      const double df2dZ = 1.0 + BS2/Z3 - gamf/W2;
      const double df2dW = B2/W3 + (2*Z - D*W)/W3 * gamf;

      const double J[2][2] = { { df1dZ,  df1dW },   // The Jacobian matrix
                               { df2dZ,  df2dW } };
      double       G[2][2];                         // The inverse Jacobian matrix

      invert_2by2_matrix(J, G);

      const double dZ = G[0][0] * f1 + G[0][1] * f2;
      const double dW = G[1][0] * f1 + G[1][1] * f2;

      double Z_new = Z - dZ;
      double W_new = W - dW;

      Z_new = (Z_new > smlZ) ? Z_new : -Z_new;
      Z_new = (Z_new < bigZ) ? Z_new :  Z;

      W_new = (W_new > smlW) ? W_new : smlW;
      W_new = (W_new < bigW) ? W_new : bigW;

      Z = Z_new;
      W = W_new;

      if (fabs(dZ/Z) + fabs(dW/W) < ERROR_TOLR)
        {
          if (Pre>=PRES_FLOOR)
            soln_found = 1;
          else
            {
              n_iter = 0;
              use_pres_floor = 1;
              W = (est) ? sqrt(S2/(D*D) + 1.0) : W_guess;
              Z = (est) ? D*W                  : P[rho]*h_guess*W_guess*W_guess;
            }
        }
      if (n_iter++ == NEWTON_MAX_ITER)
        {
          if (Pre < PRES_FLOOR)
            {
              n_iter = 0;
              use_pres_floor = 1;
              W = (est) ? sqrt(S2/(D*D) + 1.0) : W_guess;
              Z = (est) ? D*W                  : P[rho]*h_guess*W_guess*W_guess;
            }
          else
            {
              memcpy(FailedConsState, U, 8*sizeof(double));
              memcpy(FailedPrimState, P, 8*sizeof(double));
              return 1;
            }
        }
    }

  double b0 = BS * W / Z;
  P[rho] =   D/W;
  P[pre] =  (use_pres_floor) ? PRES_FLOOR : (D/W) * (Z/(D*W) - 1.0) * gamf;
  P[vx ] =  (U[Sx] + b0*U[Bx]/W) / (Z+B2);
  P[vy ] =  (U[Sy] + b0*U[By]/W) / (Z+B2);
  P[vz ] =  (U[Sz] + b0*U[Bz]/W) / (Z+B2);
  P[Bx ] =   U[Bx];
  P[By ] =   U[By];
  P[Bz ] =   U[Bz];

  return 0;
}
int prim_to_cons_point(const double *P, double *U)
{
  const double v2   =   P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
  const double B2   =   P[Bx]*P[Bx] + P[By]*P[By] + P[Bz]*P[Bz];
  const double Bv   =   P[Bx]*P[vx] + P[By]*P[vy] + P[Bz]*P[vz];
  const double W2   =   1.0 / (1.0 - v2);
  const double W    =   sqrt(W2);
  const double b0   =   W * Bv;
  const double b2   =   (B2 + b0*b0) / W2;
  const double bx   =   (P[Bx] + b0 * W*P[vx]) / W;
  const double by   =   (P[By] + b0 * W*P[vy]) / W;
  const double bz   =   (P[Bz] + b0 * W*P[vz]) / W;
  const double e    =   eos_sie(P[rho], P[pre]);
  const double p    =   P[pre];
  const double e_   =   e + 0.5 * b2 / P[rho];
  const double p_   =   p + 0.5 * b2;
  const double h_   =   1.0 + e_ + p_ / P[rho];

  U[ddd] = P[rho] * W;
  U[tau] = P[rho] * h_ * W2 - p_    - b0*b0 - U[ddd];
  U[Sx ] = P[rho] * h_ * W2 * P[vx] - b0*bx;
  U[Sy ] = P[rho] * h_ * W2 * P[vy] - b0*by;
  U[Sz ] = P[rho] * h_ * W2 * P[vz] - b0*bz;
  U[Bx ] = P[Bx ];
  U[By ] = P[By ];
  U[Bz ] = P[Bz ];
  return 0;
}

int constrained_transport_2d(double *Fx, double *Fy, int stride[4])
{
  const int N = stride[0]/8;
  double *FxBy = (double*) malloc(N*sizeof(double));
  double *FyBx = (double*) malloc(N*sizeof(double));

  double *F, *G;
  int i;

  const int sx=stride[1],sy=stride[2];
  for (i=sx; i<stride[0]-sx; i+=8)
    {
      F = &Fx[By+i];
      G = &Fy[Bx+i];

      FxBy[i/8] = (2*F[0]+F[sy]+F[-sy]-G[0]-G[sx]-G[-sy]-G[ sx-sy])*0.125;
      FyBx[i/8] = (2*G[0]+G[sx]+G[-sx]-F[0]-F[sy]-F[-sx]-F[-sx+sy])*0.125;
    }
  for (i=0; i<stride[0]; i+=8)
    {
      Fx[i+Bx] = 0.0;       Fx[i+By] = FxBy[i/8];
      Fy[i+Bx] = FyBx[i/8]; Fy[i+By] = 0.0;
    }

  free(FxBy);
  free(FyBx);

  return 0;
}
int constrained_transport_3d(double *Fx, double *Fy, double *Fz, int stride[4])
{
  const int N = stride[0]/8;
  double *FxBy = (double*) malloc(N*sizeof(double));
  double *FxBz = (double*) malloc(N*sizeof(double));

  double *FyBz = (double*) malloc(N*sizeof(double));
  double *FyBx = (double*) malloc(N*sizeof(double));

  double *FzBx = (double*) malloc(N*sizeof(double));
  double *FzBy = (double*) malloc(N*sizeof(double));

  double *F, *G, *H;
  int i;

  const int sx=stride[1],sy=stride[2],sz=stride[3];
  for (i=sx; i<stride[0]-sx; i+=8)
    {
      F = &Fx[By+i];
      G = &Fy[Bx+i];

      FxBy[i/8] = (2*F[0]+F[sy]+F[-sy]-G[0]-G[sx]-G[-sy]-G[ sx-sy])*0.125;
      FyBx[i/8] = (2*G[0]+G[sx]+G[-sx]-F[0]-F[sy]-F[-sx]-F[-sx+sy])*0.125;

      G = &Fy[Bz+i];
      H = &Fz[By+i];

      FyBz[i/8] = (2*G[0]+G[sz]+G[-sz]-H[0]-H[sy]-H[-sz]-H[ sy-sz])*0.125;
      FzBy[i/8] = (2*H[0]+H[sy]+H[-sy]-G[0]-G[sz]-G[-sy]-G[-sy+sz])*0.125;

      H = &Fz[Bx+i];
      F = &Fx[Bz+i];

      FzBx[i/8] = (2*H[0]+H[sx]+H[-sx]-F[0]-F[sz]-F[-sx]-F[ sz-sx])*0.125;
      FxBz[i/8] = (2*F[0]+F[sz]+F[-sz]-H[0]-H[sx]-H[-sz]-H[-sz+sx])*0.125;
    }
  for (i=0; i<stride[0]; i+=8)
    {
      Fx[i+Bx] = 0.0;        Fx[i+By] = FxBy[i/8];  Fx[i+Bz] = FxBz[i/8];
      Fy[i+Bx] = FyBx[i/8];  Fy[i+By] = 0.0;        Fy[i+Bz] = FyBz[i/8];
      Fz[i+Bx] = FzBx[i/8];  Fz[i+By] = FzBy[i/8];  Fz[i+Bz] = 0.0;
    }

  free(FxBy);  free(FyBz);  free(FzBx);
  free(FxBz);  free(FyBx);  free(FzBy);

  return 0;
}
