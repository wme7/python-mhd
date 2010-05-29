
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>

enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz };
enum { rho, pre, vx, vy, vz };

enum ReconstructMode { Reconstruct_PiecewiseConstant,
		       Reconstruct_PLM3Velocity,
		       Reconstruct_PLM4Velocity };

enum QuarticSolverMode { QuarticSolver_Exact,
			 QuarticSolver_Approx1,
			 QuarticSolver_Approx2 };

int dimension;
int stride[4];
double dx;

struct LibraryState
{
  int cons_to_prim_iter;
  int cons_to_prim_use_estimate;
  int cons_to_prim_verbose;

  double max_lambda;
  double adiabatic_gamma;
  double plm_theta;

  int mode_reconstruct;
  int mode_quartic_solver;
} lib_state;

double *PrimitiveArray;
double *FluxInterArray;


int set_state(struct LibraryState state)
{
  lib_state = state;
  return 0;
}
struct LibraryState get_state()
{
  return lib_state;
}

int initialize(double *P, int Nx)
{
  stride[0] = Nx;
  stride[1] = 1;
  stride[2] = 0;
  stride[3] = 0;
  dimension = 1;

  dx = 1.0 / Nx;

  PrimitiveArray = (double*) malloc(stride[0]*8*sizeof(double));
  FluxInterArray = (double*) malloc(stride[0]*8*sizeof(double));
  memcpy(PrimitiveArray, P, stride[0]*8*sizeof(double));

  printf("Initialized python ZMHD extension, version 2.\n");
  return 0;
}
int finalize()
{
  free(PrimitiveArray);
  free(FluxInterArray);
  printf("Finalized python ZMHD extension, version 2.\n");
  return 0;
}



inline double sign(double x)
{
  return (x>0)-(x<0);
}
inline double max2(double a, double b)
{
  return (a>b)?a:b;
}
inline double max3(double a, double b, double c)
{
  double ab=(a>b)?a:b;
  return (ab>c)?ab:c;
}
inline double min3(double a, double b, double c)
{
  double ab=(a<b)?a:b;
  return (ab<c)?ab:c;
}

int new_QuarticEquation(double d4, double d3, double d2, double d1, double d0);
int reconstruct_use_3vel(const double *P0, double *Pl, double *Pr);
int reconstruct_use_4vel(const double *P0, double *Pl, double *Pr);
int Fiph              (const double *U, double *F);
int dUdt_1d           (const double *U, double *L);
int prim_to_cons_point(const double *P, double *U);
int cons_to_prim_point(const double *U, double *P);
int prim_to_cons_array(const double *P, double *U);
int cons_to_prim_array(const double *U, double *P);
int rmhd_flux_and_eval(const double *U, const double *P, double *F, double *ap, double *am);

int solve_quartic_equation(double *r1, double *r2, double *r3, double *r4,
                           int *nr12, int *nr34);

int report_cons_to_prim_failure(const double *U, const double *P);
int report_nonphysical_failure (const double *U, const double *P);

double eos_pre(double Rho, double Sie)
{
  return Sie * (Rho * (lib_state.adiabatic_gamma - 1.0));
}
double eos_sie(double Rho, double Pre)
{
  return Pre / (Rho * (lib_state.adiabatic_gamma - 1.0));
}
double eos_cs2(double Rho, double Pre)
{
  double e = eos_sie(Rho, Pre);
  return lib_state.adiabatic_gamma * Pre / (Pre + Rho + Rho*e);
}


inline void invert_2by2_matrix(const double A[2][2], double B[2][2])
{
  double det = A[0][0]*A[1][1] - A[1][0]*A[0][1];

  B[0][0] =  A[1][1] / det;
  B[1][1] =  A[0][0] / det;
  B[0][1] = -A[0][1] / det;
  B[1][0] = -A[1][0] / det;
}
inline double plm_minmod(double ul, double u0, double ur)
{
  double a = lib_state.plm_theta * (u0 - ul);
  double b =               0.5   * (ur - ul);
  double c = lib_state.plm_theta * (ur - u0);

  return 0.25*fabs(sign(a) + sign(b)) * (sign(a) + sign(c))*min3(fabs(a), fabs(b), fabs(c));
}


int hll_flux(double *Ul, double *Ur, double *Fl, double *Fr,
             double eml, double epl, double emr, double epr, double *fiph)
{
  double am = min3(eml, emr, 0.0);
  double ap = max3(epl, epr, 0.0);

  int q;
  for (q=0; q<8; ++q)
    {
      fiph[q] = (ap*Fl[q] - am*Fr[q] + ap*am*(Ur[q] - Ul[q])) / (ap - am);
    }
  return 0;
}

int reconstruct_use_3vel(const double *P0, double *Pl, double *Pr)
{
  const size_t S = stride[dimension]*8;
  const size_t T = 2*S;

  Pr[rho] = P0[S+rho] - 0.5*plm_minmod(P0[ 0+rho], P0[S+rho], P0[T+rho]);
  Pl[rho] = P0[0+rho] + 0.5*plm_minmod(P0[-S+rho], P0[0+rho], P0[S+rho]);

  Pr[pre] = P0[S+pre] - 0.5*plm_minmod(P0[ 0+pre], P0[S+pre], P0[T+pre]);
  Pl[pre] = P0[0+pre] + 0.5*plm_minmod(P0[-S+pre], P0[0+pre], P0[S+pre]);

  Pr[vx] = P0[S+vx] - 0.5*plm_minmod(P0[ 0+vx], P0[S+vx], P0[T+vx]);
  Pl[vx] = P0[0+vx] + 0.5*plm_minmod(P0[-S+vx], P0[0+vx], P0[S+vx]);

  Pr[vy] = P0[S+vy] - 0.5*plm_minmod(P0[ 0+vy], P0[S+vy], P0[T+vy]);
  Pl[vy] = P0[0+vy] + 0.5*plm_minmod(P0[-S+vy], P0[0+vy], P0[S+vy]);

  Pr[vz] = P0[S+vz] - 0.5*plm_minmod(P0[ 0+vz], P0[S+vz], P0[T+vz]);
  Pl[vz] = P0[0+vz] + 0.5*plm_minmod(P0[-S+vz], P0[0+vz], P0[S+vz]);

  Pr[Bx] = P0[S+Bx] - 0.5*plm_minmod(P0[ 0+Bx], P0[S+Bx], P0[T+Bx]);
  Pl[Bx] = P0[0+Bx] + 0.5*plm_minmod(P0[-S+Bx], P0[0+Bx], P0[S+Bx]);

  Pr[By] = P0[S+By] - 0.5*plm_minmod(P0[ 0+By], P0[S+By], P0[T+By]);
  Pl[By] = P0[0+By] + 0.5*plm_minmod(P0[-S+By], P0[0+By], P0[S+By]);

  Pr[Bz] = P0[S+Bz] - 0.5*plm_minmod(P0[ 0+Bz], P0[S+Bz], P0[T+Bz]);
  Pl[Bz] = P0[0+Bz] + 0.5*plm_minmod(P0[-S+Bz], P0[0+Bz], P0[S+Bz]);

  return 0;
}
int reconstruct_use_4vel(const double *P0, double *Pl, double *Pr)
{
  const size_t S = stride[dimension]*8;
  const size_t T = 2*S;

  const double v0[4] = { 1, P0[-S+vx], P0[-S+vy], P0[-S+vz] };
  const double v1[4] = { 1, P0[ 0+vx], P0[ 0+vy], P0[ 0+vz] };
  const double v2[4] = { 1, P0[ S+vx], P0[ S+vy], P0[ S+vz] };
  const double v3[4] = { 1, P0[ T+vx], P0[ T+vy], P0[ T+vz] };

  const double W0 = 1.0 / sqrt(1.0 - (v0[1]*v0[1]+v0[2]*v0[2]+v0[3]*v0[3]));
  const double W1 = 1.0 / sqrt(1.0 - (v1[1]*v1[1]+v1[2]*v1[2]+v1[3]*v1[3]));
  const double W2 = 1.0 / sqrt(1.0 - (v2[1]*v2[1]+v2[2]*v2[2]+v2[3]*v2[3]));
  const double W3 = 1.0 / sqrt(1.0 - (v3[1]*v3[1]+v3[2]*v3[2]+v3[3]*v3[3]));

  Pr[vx] = P0[S+vx]*W2 - 0.5*plm_minmod(P0[ 0+vx]*W1, P0[S+vx]*W2, P0[T+vx]*W3);
  Pl[vx] = P0[0+vx]*W1 + 0.5*plm_minmod(P0[-S+vx]*W0, P0[0+vx]*W1, P0[S+vx]*W2);

  Pr[vy] = P0[S+vy]*W2 - 0.5*plm_minmod(P0[ 0+vy]*W1, P0[S+vy]*W2, P0[T+vy]*W3);
  Pl[vy] = P0[0+vy]*W1 + 0.5*plm_minmod(P0[-S+vy]*W0, P0[0+vy]*W1, P0[S+vy]*W2);

  Pr[vz] = P0[S+vz]*W2 - 0.5*plm_minmod(P0[ 0+vz]*W1, P0[S+vz]*W2, P0[T+vz]*W3);
  Pl[vz] = P0[0+vz]*W1 + 0.5*plm_minmod(P0[-S+vz]*W0, P0[0+vz]*W1, P0[S+vz]*W2);

  const double Wr = 1.0 / (1.0 - (Pr[vx]*Pr[vx] + Pr[vy]*Pr[vy] + Pr[vz]*Pr[vz]));
  const double Wl = 1.0 / (1.0 - (Pl[vx]*Pl[vx] + Pl[vy]*Pl[vy] + Pl[vz]*Pl[vz]));

  Pr[vx] /= Wr;  Pr[vy] /= Wr;  Pr[vz] /= Wr;
  Pl[vx] /= Wl;  Pl[vy] /= Wl;  Pl[vz] /= Wl;

  Pr[rho] = P0[S+rho] - 0.5*plm_minmod(P0[ 0+rho], P0[S+rho], P0[T+rho]);
  Pl[rho] = P0[0+rho] + 0.5*plm_minmod(P0[-S+rho], P0[0+rho], P0[S+rho]);

  Pr[pre] = P0[S+pre] - 0.5*plm_minmod(P0[ 0+pre], P0[S+pre], P0[T+pre]);
  Pl[pre] = P0[0+pre] + 0.5*plm_minmod(P0[-S+pre], P0[0+pre], P0[S+pre]);

  Pr[Bx] = P0[S+Bx] - 0.5*plm_minmod(P0[ 0+Bx], P0[S+Bx], P0[T+Bx]);
  Pl[Bx] = P0[0+Bx] + 0.5*plm_minmod(P0[-S+Bx], P0[0+Bx], P0[S+Bx]);

  Pr[By] = P0[S+By] - 0.5*plm_minmod(P0[ 0+By], P0[S+By], P0[T+By]);
  Pl[By] = P0[0+By] + 0.5*plm_minmod(P0[-S+By], P0[0+By], P0[S+By]);

  Pr[Bz] = P0[S+Bz] - 0.5*plm_minmod(P0[ 0+Bz], P0[S+Bz], P0[T+Bz]);
  Pl[Bz] = P0[0+Bz] + 0.5*plm_minmod(P0[-S+Bz], P0[0+Bz], P0[S+Bz]);

  return 0;
}


int dUdt_1d(const double *U, double *L)
{
  double *P = PrimitiveArray;
  double *F = FluxInterArray;

  int failures = cons_to_prim_array(U,P);

  dimension = 1;
  Fiph(U,F);

  size_t S = stride[dimension]*8;
  int i;
  for (i=S; i<stride[0]*8; ++i)
    {
      L[i] = -(F[i] - F[i-S]) / dx;
    }
  return failures;
}
int Fiph(const double *U, double *F)
{
  const int S = stride[dimension];

  int i;
  for (i=0; i<S*8; ++i)
    {
      F[i] = 0;
    }
  for (i=S*8; i<stride[0]*8-S*16; i+=8)
    {
      double Fl[8], Fr[8];
      double Ul[8], Ur[8];
      double Pl[8], Pr[8];
      double epl, epr, eml, emr;
      double *P = &PrimitiveArray[i];

      switch (lib_state.mode_reconstruct)
        {

        case Reconstruct_PiecewiseConstant:
          memcpy(Pl, P  , 8*sizeof(double));
          memcpy(Pr, P+8, 8*sizeof(double));
          break;

        case Reconstruct_PLM3Velocity:
          reconstruct_use_3vel(P, Pl, Pr);
          break;

        case Reconstruct_PLM4Velocity:
          reconstruct_use_4vel(P, Pl, Pr);
          break;

	default:
	  reconstruct_use_3vel(P, Pl, Pr);
          break;
        }

      prim_to_cons_point(Pl,Ul);
      prim_to_cons_point(Pr,Ur);

      rmhd_flux_and_eval(Ul, Pl, Fl, &epl, &eml);
      rmhd_flux_and_eval(Ur, Pr, Fr, &epr, &emr);

      hll_flux(Ul, Ur, Fl, Fr, eml, epl, emr, epr, &F[i]);
    }
  for (i=stride[0]*8-S*16; i<stride[0]*8; ++i)
    {
      F[i] = 0;
    }
  return 0;
}

int rmhd_flux_and_eval(const double *U, const double *P, double *F, double *ap, double *am)
{
  const double v2   =   P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
  const double B2   =   P[Bx]*P[Bx] + P[By]*P[By] + P[Bz]*P[Bz];
  const double Bv   =   P[Bx]*P[vx] + P[By]*P[vy] + P[Bz]*P[vz];
  const double W    =   1.0 / sqrt(1.0 - v2);
  const double W2   =   W*W;
  const double b0   =   W * Bv;
  const double b2   =   (B2 + b0*b0) / (W*W);
  const double bx   =   (P[Bx] + b0 * W*P[vx]) / W;
  const double by   =   (P[By] + b0 * W*P[vy]) / W;
  const double bz   =   (P[Bz] + b0 * W*P[vz]) / W;
  const double e    =   eos_sie(P[rho], P[pre]);
  const double p    =   P[pre];
  const double h    =   1.0 + e + P[pre]/P[rho];
  const double e_   =   e + 0.5 * b2 / P[rho];
  const double p_   =   p + 0.5 * b2;
  const double h_   =   1.0 + e_ + p_ / P[rho];

  switch (dimension)
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
  switch (dimension)
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

  const double W4   =  W2*W2;
  const double cs2  =  eos_cs2(P[rho],P[pre]);
  const double V2   =  pow(vi,2.);
  const double V3   =  pow(vi,3.);
  const double V4   =  pow(vi,4.);

  const double K  =    P[rho]*h * (1./cs2-1.) * W4;
  const double L  =  -(P[rho]*h +  b2/cs2)    * W2;

  const double A4 =    K    - L        -   b0*b0;
  const double A3 = -4*K*vi + L*vi*2   + 2*b0*bi;
  const double A2 =  6*K*V2 + L*(1-V2) +   b0*b0 - bi*bi;
  const double A1 = -4*K*V3 - L*vi*2   - 2*b0*bi;
  const double A0 =    K*V4 + L*V2     +   bi*bi;

  double r1, r2, r3, r4;
  int nr12, nr34, nr;

  new_QuarticEquation(A4,A3,A2,A1,A0);
  switch (lib_state.mode_quartic_solver)
    {
    case (QuarticSolver_Exact):
      nr = solve_quartic_equation (&r1, &r2, &r3, &r4, &nr12, &nr34);
      break;
    case (QuarticSolver_Approx1):
      nr = solve_quartic_approx1  (&r1, &r2, &r3, &r4, &nr12, &nr34);
      break;
    case (QuarticSolver_Approx2):
      nr = solve_quartic_approx2  (&r1, &r2, &r3, &r4, &nr12, &nr34);
      break;
    }

  double ap12 = (r1>r2) ? r1 : r2;
  double ap34 = (r3>r4) ? r3 : r4;

  double am12 = (r1<r2) ? r1 : r2;
  double am34 = (r3<r4) ? r3 : r4;

  *ap = (nr==2) ? ((nr12==2) ? ap12 : ap34) : ((ap12>ap34) ? ap12 : ap34);
  *am = (nr==2) ? ((nr12==2) ? am12 : am34) : ((am12<am34) ? am12 : am34);
  return 0;
}

int cons_to_prim_point(const double *U, double *P)
{
  static const double PRES_FLOOR = 1e-10;
  static const double ERROR_TOLR = 1e-6;
  static const int NEWTON_MAX_ITER = 15;

  // Quantites known from conserved
  const double gamf  = (lib_state.adiabatic_gamma - 1.0) / lib_state.adiabatic_gamma;
  const double D     = U[ddd];
  const double Tau   = U[tau];
  const double S2    = U[Sx]*U[Sx] + U[Sy]*U[Sy] + U[Sz]*U[Sz];
  const double B2    = U[Bx]*U[Bx] + U[By]*U[By] + U[Bz]*U[Bz];
  const double BS    = U[Bx]*U[Sx] + U[By]*U[Sy] + U[Bz]*U[Sz];
  const double BS2   = BS*BS;

  int est            = lib_state.cons_to_prim_use_estimate;
  int verbose        = lib_state.cons_to_prim_verbose;
  int use_pres_floor = 0;
  int soln_found     = 0;
  int n_iter         = 0;

  // Quantites guessed from existing primitives, or estimated from conserved
  double h_guess = 1 + eos_sie(P[rho], P[pre]) + P[pre]/P[rho];
  double W_guess = 1.0 / sqrt(1 - (P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz]));
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

      const double Pre = (use_pres_floor) ? PRES_FLOOR : (D/W) * (Z/(D*W) - 1) * gamf;

      const double f1 = -S2  +  (Z+B2)*(Z+B2)*(W2-1)/W2 - (2*Z+B2)*BS2/Z2;      // eqn (84)
      const double f2 = -Tau +   Z+B2 - Pre - 0.5*B2/W2 -      0.5*BS2/Z2 - D;  // eqn (85)

      const double df1dZ = 2*(B2+Z)*(BS2*W2 + (W2-1)*Z3) / (W2*Z3);
      const double df1dW = 2*(B2+Z)*(B2+Z) / W3;
      const double df2dZ = 1 + BS2/Z3 - gamf/W2;
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

      if (verbose)
        {
          printf("iteration number: %d, Pre = %8.6e, W = %8.6e, Z = %8.6e, f1 = %8.6e, f2 = %8.6e\n",
                 n_iter, Pre, W, Z, f1, f2);
        }
      if (fabs(dZ/Z) + fabs(dW/W) < ERROR_TOLR)
        {
          if (Pre>=PRES_FLOOR)
            soln_found = 1;
          else
            {
              n_iter = 0;
              use_pres_floor = 1;
              W = (est) ? sqrt(S2/(D*D) + 1) : W_guess;
              Z = (est) ? D*W                : P[rho]*h_guess*W_guess*W_guess;
            }
        }
      if (n_iter++ == NEWTON_MAX_ITER)
        {
          if (Pre < PRES_FLOOR)
            {
              n_iter = 0;
              use_pres_floor = 1;
              W = (est) ? sqrt(S2/(D*D) + 1) : W_guess;
              Z = (est) ? D*W                : P[rho]*h_guess*W_guess*W_guess;
            }
          else
            {
              return 1;
            }
        }
    }
  lib_state.cons_to_prim_iter += n_iter;

  double b0 = BS * W / Z;
  P[rho] =   D/W;
  P[pre] =  (use_pres_floor) ? PRES_FLOOR : (D/W) * (Z/(D*W) - 1) * gamf;
  P[vx ] =  (U[Sx] + b0*U[Bx]/W) / (Z+B2);
  P[vy ] =  (U[Sy] + b0*U[By]/W) / (Z+B2);
  P[vz ] =  (U[Sz] + b0*U[Bz]/W) / (Z+B2);
  P[Bx ] =   U[Bx];
  P[By ] =   U[By];
  P[Bz ] =   U[Bz];

  return 0;
}
int cons_to_prim_array(const double *U, double *P)
{
  int failures = 0;
  int i;

  if (P != PrimitiveArray) // Replace the input array if not the same
    memcpy(P, PrimitiveArray, 8*stride[0]*sizeof(double));

  for (i=0; i<stride[0]*8; i+=8)
    {
      const double *Ui = &U[i];
      double       *Pi = &P[i];

      if (Ui[ddd] < 0.0 || Ui[tau] < 0.0 ||
          Pi[pre] < 0.0 || Pi[rho] < 0.0)
        {
	  failures += 1;
        }
      else if (cons_to_prim_point(Ui, Pi))
        {
	  failures += 1;
        }
    }
  return failures;
}
int prim_to_cons_point(const double *P, double *U)
{
  const double v2   =   P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
  const double B2   =   P[Bx]*P[Bx] + P[By]*P[By] + P[Bz]*P[Bz];
  const double Bv   =   P[Bx]*P[vx] + P[By]*P[vy] + P[Bz]*P[vz];
  const double W2   =   1.0 / (1.0 - v2);
  const double W    =   sqrt(W2);
  const double b0   =   W * Bv;
  const double b2   =   (B2 + b0*b0) / (W*W);
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
int prim_to_cons_array(const double *P, double *U)
{
  int i;
  for (i=0; i<stride[0]*8; i+=8)
    {
      prim_to_cons_point(&P[i], &U[i]);
    }
  return 0;
}


int report_nonphysical_failure(const double *U, const double *P)
{
  printf("\n\n ## FATAL ERROR! ZMHD encountered nonphysical quantities!\n\n");
  printf("D   = %f\n", U[ddd]);
  printf("tau = %f\n", U[tau]);
  printf("Sx  = %f\n", U[Sx ]);
  printf("Sy  = %f\n", U[Sy ]);
  printf("Sz  = %f\n", U[Sz ]);
  printf("Bx  = %f\n", U[Bx ]);
  printf("By  = %f\n", U[By ]);
  printf("Bz  = %f\n", U[Bz ]);

  printf("Pressure = %f\n", P[pre]);
  printf("Density  = %f\n", P[rho]);
  printf("vx       = %f\n", P[vx]);
  printf("vy       = %f\n", P[vy]);
  printf("vz       = %f\n", P[vz]);
  return 0;
}


int report_cons_to_prim_failure(const double *U, const double *P)
{
  printf("\n\n ## FATAL ERROR! ZMHD encountered primitive recovery failure!\n\n");
  printf("The conserved state which choked the solver was:\n");

  double P_tmp[8];
  report_nonphysical_failure(U,P);
  memcpy(P_tmp, P, 8*sizeof(double));

  printf("\n\nHere is the output of the solver as it is failing:\n\n");
  lib_state.cons_to_prim_verbose = 1;
  cons_to_prim_point(U,P_tmp);
  return 0;
}
