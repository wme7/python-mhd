
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz };
enum { rho, pre, vx, vy, vz };

typedef struct
{
  double rho, pre, vx, vy, vz, Bx, By, Bz;
  double ux, uy, uz;

  double v2, B2, Bv, W2, W;
  double b0, bx, by, bz, b2;
  double e, h, cs2;
  double e_, p_, h_;
} AuxQnt;


int dimension;
int stride[4];
double dx;

enum ReconstructMode { Reconstruct_PiecewiseConstant,
		       Reconstruct_PLM3Velocity,
		       Reconstruct_PLM4Velocity };

enum QuarticSolverMode { QuarticSolver_Exact,
			 QuarticSolver_Approx1,
			 QuarticSolver_Approx2 };

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

AuxQnt *AuxlQuantArray;
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

  AuxlQuantArray = (AuxQnt*) malloc(stride[0]*1*sizeof(AuxQnt));
  FluxInterArray = (double*) malloc(stride[0]*8*sizeof(double));

  printf("Initialized python ZMHD extension, version 1.\n");
  return 0;
}
int finalize()
{
  free(AuxlQuantArray);
  free(FluxInterArray);
  printf("Finalized python ZMHD extension, version 1.\n");
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

int solve_quartic_equation(double d4, double d3, double d2, double d1, double d0,
                           double *r1, double *r2, double *r3, double *r4,
                           int *nr12, int *nr34);

int report_cons_to_prim_failure(const double *U, const AuxQnt *A);
int report_nonphysical_failure (const double *U, const AuxQnt *A);
int cons_to_prim_array(const double *U, double *P);
int prim_to_cons_array(const double *P, double *U);
int reconstruct(const AuxQnt *A, AuxQnt *Al, AuxQnt *Ar, double *Ul, double *Ur);
int Fiph(const double *U, double *fiph);
int rmhd_flux(const double *U, const AuxQnt *A, double *F);
int cons_to_aux(const double *U, AuxQnt *A, int use_guess, int verbose);
int recover_aux_array(const double *U);
int eigenvalues(const double *U, const AuxQnt *A, double *ap, double *am);


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
inline void complete_aux(AuxQnt *A)
{
  A->W2   =   1.0 + A->ux*A->ux + A->uy*A->uy + A->uz*A->uz;
  A->W    =   sqrt(A->W2);
  A->v2   =   1.0 - 1.0/A->W2;

  A->vx   =   A->ux/A->W;
  A->vy   =   A->uy/A->W;
  A->vz   =   A->uz/A->W;

  A->B2   =   A->Bx*A->Bx + A->By*A->By + A->Bz*A->Bz;
  A->Bv   =   A->Bx*A->vx + A->By*A->vy + A->Bz*A->vz;
  A->b0   =   A->W * A->Bv;
  A->bx   =   (A->Bx + A->b0 * A->ux) / A->W;
  A->by   =   (A->By + A->b0 * A->uy) / A->W;
  A->bz   =   (A->Bz + A->b0 * A->uz) / A->W;
  A->b2   =   (A->B2 + A->b0*A->b0) / A->W2;

  A->e    =  eos_sie(A->rho, A->pre);
  A->cs2  =  eos_cs2(A->rho, A->pre);
  A->h    =  1.0 + A->e + A->pre/A->rho;

  A->e_   =   A->e   + 0.5  * A->b2 / A->rho;
  A->p_   =   A->pre + 0.5  * A->b2;
  A->h_   =   1.0   + A->e_ + A->p_ / A->rho;
}
inline void aux_to_cons(const AuxQnt *A, double *U)
{
  U[ddd] = A->rho*A->W;
  U[tau] = A->rho*A->h_ * A->W2 - A->p_ - A->b0*A->b0 - U[ddd];
  U[Sx ] = A->rho*A->h_ * A->W2 * A->vx - A->b0*A->bx;
  U[Sy ] = A->rho*A->h_ * A->W2 * A->vy - A->b0*A->by;
  U[Sz ] = A->rho*A->h_ * A->W2 * A->vz - A->b0*A->bz;
  U[Bx ] = A->Bx;
  U[By ] = A->By;
  U[Bz ] = A->Bz;
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
  return 1;
}

int reconstruct(const AuxQnt *A, AuxQnt *Al, AuxQnt *Ar, double *Ul, double *Ur)
{
  const size_t S = stride[dimension];
  const size_t T = 2*S;

  Ar->rho = A[S].rho - .5*plm_minmod(A[ 0].rho, A[S].rho, A[T].rho);
  Al->rho = A[0].rho + .5*plm_minmod(A[-S].rho, A[0].rho, A[S].rho);

  Ar->pre = A[S].pre - .5*plm_minmod(A[ 0].pre, A[S].pre, A[T].pre);
  Al->pre = A[0].pre + .5*plm_minmod(A[-S].pre, A[0].pre, A[S].pre);

  Ar->ux = A[S].ux - .5*plm_minmod(A[ 0].ux, A[S].ux, A[T].ux);
  Al->ux = A[0].ux + .5*plm_minmod(A[-S].ux, A[0].ux, A[S].ux);

  Ar->uy = A[S].uy - .5*plm_minmod(A[ 0].uy, A[S].uy, A[T].uy);
  Al->uy = A[0].uy + .5*plm_minmod(A[-S].uy, A[0].uy, A[S].uy);

  Ar->uz = A[S].uz - .5*plm_minmod(A[ 0].uz, A[S].uz, A[T].uz);
  Al->uz = A[0].uz + .5*plm_minmod(A[-S].uz, A[0].uz, A[S].uz);

  Ar->Bx = A[S].Bx - .5*plm_minmod(A[ 0].Bx, A[S].Bx, A[T].Bx);
  Al->Bx = A[0].Bx + .5*plm_minmod(A[-S].Bx, A[0].Bx, A[S].Bx);

  Ar->By = A[S].By - .5*plm_minmod(A[ 0].By, A[S].By, A[T].By);
  Al->By = A[0].By + .5*plm_minmod(A[-S].By, A[0].By, A[S].By);

  Ar->Bz = A[S].Bz - .5*plm_minmod(A[ 0].Bz, A[S].Bz, A[T].Bz);
  Al->Bz = A[0].Bz + .5*plm_minmod(A[-S].Bz, A[0].Bz, A[S].Bz);

  complete_aux(Al);
  complete_aux(Ar);
  aux_to_cons(Al, Ul);
  aux_to_cons(Ar, Ur);

  return 0;
}

int Fiph(const double *U, double *fiph)
{
  const int S = stride[dimension];
  int i;
  for (i=0; i<S*8; ++i)
    {
      fiph[i] = 0;
    }
  for (i=S*8; i<stride[0]*8-S*16; i+=8)
    {
      double Fl[8], Fr[8];
      double Ul[8], Ur[8];
      AuxQnt Al, Ar;
      double epl, epr, eml, emr;
      reconstruct(&AuxlQuantArray[i/8], &Al, &Ar, Ul, Ur);

      rmhd_flux(Ul, &Al, Fl);
      rmhd_flux(Ur, &Ar, Fr);
      eigenvalues(Ul, &Al, &epl, &eml);
      eigenvalues(Ur, &Ar, &epr, &emr);

      hll_flux(Ul, Ur, Fl, Fr, eml, epl, emr, epr, &fiph[i]);
    }
  for (i=stride[0]*8-S*16; i<stride[0]*8; ++i)
    {
      fiph[i] = 0;
    }
  return 0;
}

int rmhd_flux(const double *U, const AuxQnt *A, double *F)
{
  switch (dimension)
    {
    case 1:
      F[ddd] = U[ddd] * A->vx;
      F[tau] = U[tau] * A->vx - A->b0*A->Bx / A->W + A->p_*A->vx;
      F[Sx ] = U[Sx ] * A->vx - A->bx*A->Bx / A->W + A->p_;
      F[Sy ] = U[Sy ] * A->vx - A->by*A->Bx / A->W;
      F[Sz ] = U[Sz ] * A->vx - A->bz*A->Bx / A->W;
      F[Bx ] = 0.0;
      F[By ] = A->vx*A->By - A->vy*A->Bx;
      F[Bz ] = A->vx*A->Bz - A->vz*A->Bx;
      break;
    case 2:
      F[ddd] = U[ddd] * A->vy;
      F[tau] = U[tau] * A->vy - A->b0*A->By / A->W + A->p_*A->vy;
      F[Sx ] = U[Sx ] * A->vy - A->bx*A->By / A->W;
      F[Sy ] = U[Sy ] * A->vy - A->by*A->By / A->W + A->p_;
      F[Sz ] = U[Sz ] * A->vy - A->bz*A->By / A->W;
      F[Bx ] = A->vy*A->Bx - A->vx*A->By;
      F[By ] = 0.0;
      F[Bz ] = A->vy*A->Bz - A->vz*A->By;
      break;
    case 3:
      F[ddd] = U[ddd] * A->vz;
      F[tau] = U[tau] * A->vz - A->b0*A->Bz / A->W + A->p_*A->vz;
      F[Sx ] = U[Sx ] * A->vz - A->bx*A->Bz / A->W;
      F[Sy ] = U[Sy ] * A->vz - A->by*A->Bz / A->W;
      F[Sz ] = U[Sz ] * A->vz - A->bz*A->Bz / A->W + A->p_;
      F[Bx ] = A->vz*A->Bx - A->vx*A->Bz;
      F[By ] = A->vz*A->By - A->vy*A->Bz;
      F[Bz ] = 0.0;
      break;
    }
  return 0;
}

int eigenvalues(const double *U, const AuxQnt *A, double *ap, double *am)
{
  double vi=0, bi=0;
  switch (dimension)
    {
    case 1:
      vi = A->vx; bi = A->bx;
      break;
    case 2:
      vi = A->vy; bi = A->by;
      break;
    case 3:
      vi = A->vz; bi = A->bz;
      break;
    }

  const double b0   =  A->b0;
  const double b2   =  A->b2;
  const double W4   =  A->W2*A->W2;
  const double cs2  =  A->cs2;
  const double V2   =  pow(vi,2.);
  const double V3   =  pow(vi,3.);
  const double V4   =  pow(vi,4.);

  const double K  =    A->rho*A->h * (1./cs2-1.) *    W4;
  const double L  =  -(A->rho*A->h +  b2/cs2)    * A->W2;

  const double A4 =    K    - L        -   b0*b0;
  const double A3 = -4*K*vi + L*vi*2   + 2*b0*bi;
  const double A2 =  6*K*V2 + L*(1-V2) +   b0*b0 - bi*bi;
  const double A1 = -4*K*V3 - L*vi*2   - 2*b0*bi;
  const double A0 =    K*V4 + L*V2     +   bi*bi;

  double r1, r2, r3, r4;
  int nr12, nr34;
  int nr = solve_quartic_equation(A4,A3,A2,A1,A0, &r1, &r2, &r3, &r4, &nr12, &nr34);

  double ap12 = (r1>r2) ? r1 : r2;
  double ap34 = (r3>r4) ? r3 : r4;

  double am12 = (r1<r2) ? r1 : r2;
  double am34 = (r3<r4) ? r3 : r4;

  *ap = (nr==2) ? ((nr12==2) ? ap12 : ap34) : ((ap12>ap34) ? ap12 : ap34);
  *am = (nr==2) ? ((nr12==2) ? am12 : am34) : ((am12<am34) ? am12 : am34);

  return 0;
}

int cons_to_aux(const double *U, AuxQnt *A, int use_guess, int verbose)
{
  static const double PRES_FLOOR = 1e-10;
  static const double ERROR_TOLR = 1e-6;
  static const int NEWTON_MAX_ITER = 25;

  const double gamf  = (lib_state.adiabatic_gamma - 1.0) / lib_state.adiabatic_gamma;
  const double D     = U[ddd];
  const double Tau   = U[tau];
  const double S2    = U[Sx]*U[Sx] + U[Sy]*U[Sy] + U[Sz]*U[Sz];
  const double B2    = U[Bx]*U[Bx] + U[By]*U[By] + U[Bz]*U[Bz];
  const double BS    = U[Bx]*U[Sx] + U[By]*U[Sy] + U[Bz]*U[Sz];
  const double BS2   = BS*BS;

  double W = (use_guess) ? sqrt(S2/(D*D) + 1) : A->W;
  double Z = (use_guess) ? D*W                : A->rho*A->h*A->W2;

  int use_pres_floor = 0;
  int soln_found = 0;
  int n_iter = 0;

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
              W =   A->W;
              Z = D*A->W;
            }
        }
      if (n_iter++ == NEWTON_MAX_ITER)
        {
          if (Pre < PRES_FLOOR)
            {
              n_iter = 0;
              use_pres_floor = 1;
              W =   A->W;
              Z = D*A->W;
            }
          else
            {
              printf("failure!\n");
              return 0;
            }
        }
    }
  lib_state.cons_to_prim_iter += n_iter;

  A->b0  =   BS * W / Z;
  A->rho =   D/W;
  A->pre =  (use_pres_floor) ? PRES_FLOOR : (D/W) * (Z/(D*W) - 1) * gamf;
  A->vx  =  (U[Sx] + A->b0*U[Bx]/W) / (Z+B2);
  A->vy  =  (U[Sy] + A->b0*U[By]/W) / (Z+B2);
  A->vz  =  (U[Sz] + A->b0*U[Bz]/W) / (Z+B2);

  A->Bx  =   U[Bx];
  A->By  =   U[By];
  A->Bz  =   U[Bz];

  A->v2  =   A->vx*A->vx + A->vy*A->vy + A->vz*A->vz;
  A->B2  =   A->Bx*A->Bx + A->By*A->By + A->Bz*A->Bz;
  A->Bv  =   A->Bx*A->vx + A->By*A->vy + A->Bz*A->vz;
  A->W2  =   1.0 / (1.0 - A->v2);
  A->W   =   sqrt(A->W2);

  A->ux  =  A->vx * A->W;
  A->uy  =  A->vy * A->W;
  A->uz  =  A->vz * A->W;

  A->bx  =   (A->Bx + A->b0 * A->ux) / A->W;
  A->by  =   (A->By + A->b0 * A->uy) / A->W;
  A->bz  =   (A->Bz + A->b0 * A->uz) / A->W;
  A->b2  =   (A->B2 + A->b0 * A->b0) / A->W2;

  A->e   =  eos_sie(A->rho, A->pre);
  A->cs2 =  eos_cs2(A->rho, A->pre);
  A->h   =  1.0 + A->e + A->pre/A->rho;
  A->e_  =  A->e   + 0.5  * A->b2 / A->rho;
  A->p_  =  A->pre + 0.5  * A->b2;
  A->h_  =  1.0   + A->e_ + A->p_ / A->rho;

  return 1;
}

int recover_aux_array(const double *U)
{
  static int first_call_ever = 1;
  int i;
  for (i=0; i<stride[0]*8; i+=8)
    {
      const double *Ui = &U[i];
      AuxQnt       *Ai = &AuxlQuantArray[i/8];

      if (Ui[ddd] < 0.0 || Ui[tau] < 0.0 ||
          Ai->pre < 0.0 || Ai->rho < 0.0)
        {
          report_nonphysical_failure(Ui, Ai);
	  return 1;
        }
      if (!cons_to_aux(Ui, Ai, first_call_ever, 0))
        {
          report_cons_to_aux_failure(Ui, Ai);
	  return 1;
        }
    }
  first_call_ever = 0;

  return 0;
}

int dUdt_1d(const double *U, double *L)
{
  recover_aux_array(U);
  dimension = 1;
  Fiph(U, FluxInterArray);

  size_t Si = stride[1]*8;

  int i;
  for (i=Si; i<stride[0]*8; ++i)
    {
      L[i] = -(FluxInterArray[i] - FluxInterArray[i-Si]) / dx;
    }
  return 0;
}

int report_nonphysical_failure(const double *U, const AuxQnt *A)
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

  printf("Pressure = %f\n", A->pre);
  printf("Density  = %f\n", A->rho);
  printf("e        = %f\n", A->e);
  printf("cs2      = %f\n", A->cs2);
  printf("vx       = %f\n", A->vx);
  printf("vy       = %f\n", A->vy);
  printf("vz       = %f\n", A->vz);
  return 0;
}


int report_cons_to_aux_failure(const double *U, const AuxQnt *A)
{
  printf("\n\n ## FATAL ERROR! ZMHD encountered primitive recovery failure!\n\n");
  printf("The conserved state which choked the solver was:\n");

  AuxQnt A_tmp = *A;
  report_nonphysical_failure(U,&A_tmp);

  printf("\n\nHere is the output of the solver as it is failing:\n\n");
  cons_to_aux(U, &A_tmp, 0, 1);
  return 0;
}


int cons_to_prim_array(const double *U, double *P)
{
  recover_aux_array(U);
  int i;
  for (i=0; i<stride[0]*8; i+=8)
    {
      AuxQnt *A = &AuxlQuantArray[i/8];
      P[i + rho] = A->rho;
      P[i + pre] = A->pre;
      P[i + vx ] = A->vx;
      P[i + vy ] = A->vy;
      P[i + vz ] = A->vz;
      P[i + Bx ] = A->Bx;
      P[i + By ] = A->By;
      P[i + Bz ] = A->Bz;
    }
  return 0;
}
int prim_to_cons_array(const double *P, double *U)
{
  int i;
  for (i=0; i<stride[0]*8; i+=8)
    {
      const double *Pi = &P[i];
      double       *Ui = &U[i];

      const double v2   =   Pi[vx]*Pi[vx] + Pi[vy]*Pi[vy] + Pi[vz]*Pi[vz];
      const double B2   =   Pi[Bx]*Pi[Bx] + Pi[By]*Pi[By] + Pi[Bz]*Pi[Bz];
      const double Bv   =   Pi[Bx]*Pi[vx] + Pi[By]*Pi[vy] + Pi[Bz]*Pi[vz];
      const double W    =   1.0 / sqrt(1.0 - v2);
      const double b0   =   W * Bv;
      const double b2   =   (B2 + b0*b0) / (W*W);
      const double bx   =   (Pi[Bx] + b0 * W*Pi[vx]) / W;
      const double by   =   (Pi[By] + b0 * W*Pi[vy]) / W;
      const double bz   =   (Pi[Bz] + b0 * W*Pi[vz]) / W;
      const double e    =   eos_sie(Pi[rho], Pi[pre]);
      const double p    =   Pi[pre];
      const double e_   =   e + 0.5 * b2 / Pi[rho];
      const double p_   =   p + 0.5 * b2;
      const double h_   =   1.0 + e_ + p_ / Pi[rho];

      Ui[ddd] = Pi[rho] * W;
      Ui[tau] = Pi[rho] * h_ * W*W - p_    - b0*b0 - Ui[ddd];
      Ui[Sx ] = Pi[rho] * h_ * W*W * Pi[vx] - b0*bx;
      Ui[Sy ] = Pi[rho] * h_ * W*W * Pi[vy] - b0*by;
      Ui[Sz ] = Pi[rho] * h_ * W*W * Pi[vz] - b0*bz;
      Ui[Bx ] = Pi[Bx ];
      Ui[By ] = Pi[By ];
      Ui[Bz ] = Pi[Bz ];
    }
  return 0;
}

int solve_quartic_equation(double  d4, double  d3, double  d2, double  d1, double d0,
                           double *r1, double *r2, double *r3, double *r4,
                           int *nr12, int *nr34)
{
  double a3 = d3/d4;
  double a2 = d2/d4;
  double a1 = d1/d4;
  double a0 = d0/d4;

  double au2 = -a2;
  double au1 = (a1*a3 - 4.0*a0) ;
  double au0 = 4.0*a0*a2 - a1*a1 - a0*a3*a3;

  double x1, x2, x3;
  int nr = solve_cubic_equation(1.0, au2, au1, au0, &x1, &x2, &x3);

  double u1;
  if (nr==1) u1 = x1;
  else u1 = (x1>x3) ? x1 : x3;

  double R2 = 0.25*a3*a3 + u1 - a2;
  double R = (R2>0.0) ? sqrt(R2) : 0.0;

  double D2, E2;
  if (R != 0.0)
    {
      double foo1 = 0.75*a3*a3 - R2 - 2.0*a2;
      double foo2 = 0.25*(4.0*a3*a2 - 8.0*a1 - a3*a3*a3) / R;
      D2 = foo1 + foo2;
      E2 = foo1 - foo2;
    }
  else
    {
      double foo1 = 0.75*a3*a3 - 2.0*a2;
      double foo2 = 2.0 * sqrt(u1*u1 - 4.0*a0);
      D2 = foo1 + foo2;
      E2 = foo1 - foo2;
    }

  if (D2 >= 0.0)
    {
      double D = sqrt(D2);
      *r1 = -0.25*a3 + 0.5*R - 0.5*D;
      *r2 = -0.25*a3 + 0.5*R + 0.5*D;
      *nr12 = 2;
    }
  else
    {
      *r1 = *r2 = -0.25*a3 + 0.5*R;
      *nr12 = 0;
    }

  if (E2 >= 0.0)
    {
      double E = sqrt(E2);
      *r3 = -0.25*a3 - 0.5*R - 0.5*E;
      *r4 = -0.25*a3 - 0.5*R + 0.5*E;
      *nr34 = 2;
    }
  else
    {
      *r3 = *r4 = -0.25*a3 - 0.5*R;
      *nr34 = 0;
    }

  return *nr12 + *nr34;
}

int solve_cubic_equation(double  c3, double  c2,  double c1, double c0,
                         double *x1, double *x2, double *x3)
{
  double a2 = c2/c3;
  double a1 = c1/c3;
  double a0 = c0/c3;

  double q = a1/3.0 - a2*a2/9.0;
  double r = (a1*a2 - 3.0*a0)/6.0 - a2*a2*a2 / 27.0;
  double delta = q*q*q + r*r;

  if (delta>0.0)
    {
      double s1 = r + sqrt(delta);
      s1 = (s1>=0.0) ? pow(s1,1./3.) : -pow(-s1,1./3.);

      double s2 = r - sqrt(delta);
      s2 = (s2>=0.0) ? pow(s2,1./3.) : -pow(-s2,1./3.);

      *x1 = (s1+s2) - a2/3.0;
      *x2 = *x3 = -0.5 * (s1+s2) - a2/3.0;

      return 1;
    }
  else if (delta < 0.0)
    {
      double theta = acos(r/sqrt(-q*q*q)) / 3.0;
      double costh = cos(theta);
      double sinth = sin(theta);
      double sq = sqrt(-q);

      *x1 = 2.0*sq*costh - a2/3.0;
      *x2 = -sq*costh - a2/3.0 - sqrt(3.) * sq * sinth;
      *x3 = -sq*costh - a2/3.0 + sqrt(3.) * sq * sinth;

      return 3;
    }
  else
    {
      double s = (r>=0.0) ? pow(r,1./3.) : -pow(-r,1./3.);
      *x1 = 2.0*s - a2/3.0;
      *x2 = *x3 = -s - a2/3.0;

      return 3;
    }
}
