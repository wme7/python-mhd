
/*------------------------------------------------------------------------------
 * FILE: integrate.c
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
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
int integrate_init(int N[4], double L[4], int num_comp, int num_dims);
int integrate_free();
int get_failure_mask(int *M);

int advance_state_fwd_euler   (double *P, double dt);
int advance_state_midpoint    (double *P, double dt);
int advance_state_RK3         (double *P, double dt);
int advance_state_ctu_hancock (double *P, double dt);

int prim_to_cons_array(double *P, double *U, int N);
int cons_to_prim_array(double *U, double *P, int N);

int NQ,ND;

/*------------------------------------------------------------------------------
 *
 * External Dependecies
 *
 */
int hll_flux(const double *pl, const double *pr, double *U, double *F,
             double s, int dim);
int hllc_flux(const double *pl, const double *pr, double *U, double *F,
             double s, int dim);

int flux_and_eval(const double *U, const double *P, double *F,
                  double *ap, double *am, int dim);

int prim_to_cons_point(const double *P, double *U);
int cons_to_prim_point(const double *U, double *P);

int constrained_transport_2d(double *Fx, double *Fy, int stride[4]);
int constrained_transport_3d(double *Fx, double *Fy, double *Fz, int stride[4]);

/*------------------------------------------------------------------------------
 *
 * Internal Dependecies
 *
 */
static int reconstruct_plm(const double *P0, double *Pl, double *Pr, int S);
static int intercell_flux_sweep(const double *P, double *F, int dim);

static int drive_sweeps_1d(const double *P, double *L);
static int drive_sweeps_2d(const double *P, double *L);
static int drive_sweeps_3d(const double *P, double *L);

static inline double plm_minmod(double ul, double u0, double ur);
static inline double MC_limiter(double ul, double u0, double ur);
static inline double harmonic_mean(double ul, double u0, double ur);
/*------------------------------------------------------------------------------
 *
 * Private Data
 *
 */
typedef double (*SlopeLimiter)(double, double, double);
typedef int (*RiemannSolver)(const double*, const double*, double*, double*, double, int);
typedef int (*SweepsDriver)(const double*, double*);

static SlopeLimiter slope_limiter = plm_minmod;
static RiemannSolver intercell_flux = hll_flux;
static SweepsDriver drive_sweeps;

#define MAXNQ 8 // Used for static array initialization
static int *failure_mask;
static int cons_to_prim_fail_fatal=0;
static int stride[4];
static double plm_theta=2.0;
static double dx,dy,dz;

/*------------------------------------------------------------------------------
 *
 * Inline functions
 *
 */
static inline double sign(double x)
{
  return (x>0)-(x<0);
}
static inline double max2(double a, double b)
{
  return (a>b)?a:b;
}
static inline double max3(double a, double b, double c)
{
  double ab=(a>b)?a:b;
  return (ab>c)?ab:c;
}
static inline double min3(double a, double b, double c)
{
  double ab=(a<b)?a:b;
  return (ab<c)?ab:c;
}
static inline double plm_minmod(double ul, double u0, double ur)
{
  const double a = plm_theta * (u0 - ul);
  const double b =     0.5   * (ur - ul);
  const double c = plm_theta * (ur - u0);
  return 0.25*fabs(sign(a)+sign(b))*(sign(a)+sign(c))*min3(fabs(a),fabs(b),fabs(c));
}
static inline double MC_limiter(double ul, double u0, double ur)
{
  const double qp = ur - u0;
  const double qm = u0 - ul;
  const double si = 0.5*(sign(qp) + sign(qm));
  return si*min3(2*fabs(qp), 2*fabs(qm), 0.5*fabs(ur-ul));
}
static inline double harmonic_mean(double ul, double u0, double ur)
{
  const double qp = ur - u0;
  const double qm = u0 - ul;
  return (sign(qp)==sign(qm)) ? 2*qp*qm/(ur-ul) : 0;
}

/*------------------------------------------------------------------------------
 *
 * Private Functions Definitions
 *
 */
static int reconstruct_plm(const double *P0, double *Pl, double *Pr, int S)
{
  int i,T=2*S;
  for (i=0; i<NQ; ++i)
    {
      Pr[i] = P0[S+i] - 0.5 * slope_limiter(P0[ 0+i], P0[S+i], P0[T+i]);
      Pl[i] = P0[0+i] + 0.5 * slope_limiter(P0[-S+i], P0[0+i], P0[S+i]);
    }
  return 0;
}
static int intercell_flux_sweep(const double *P, double *F, int dim)
{
  int i,S=stride[dim];
  double Pl[MAXNQ], Pr[MAXNQ];
  for (i=S; i<stride[0]-2*S; i+=NQ)
    {
      reconstruct_plm(&P[i], Pl, Pr, S);
      intercell_flux(Pl, Pr, 0, &F[i], 0.0, dim);
    }
  return 0;
}

static int drive_sweeps_1d(const double *P, double *L)
{
  double *F = (double*) malloc(stride[0]*sizeof(double));

  int i,sx=stride[1];
  intercell_flux_sweep(P,F,1);

  for (i=sx; i<stride[0]; ++i)
    {
      L[i] = -(F[i]-F[i-sx])/dx;
    }

  free(F);

  return 0;
}
static int drive_sweeps_2d(const double *P, double *L)
{
  double *F = (double*) malloc(stride[0]*sizeof(double));
  double *G = (double*) malloc(stride[0]*sizeof(double));

  int i,sx=stride[1],sy=stride[2];
  intercell_flux_sweep(P,F,1);
  intercell_flux_sweep(P,G,2);

  for (i=sx; i<stride[0]; ++i)
    {
      L[i] = -(F[i]-F[i-sx])/dx - (G[i]-G[i-sy])/dy;
    }

  constrained_transport_2d(F,G,stride);

  free(F);
  free(G);

  return 0;
}
static int drive_sweeps_3d(const double *P, double *L)
{
  double *F = (double*) malloc(stride[0]*sizeof(double));
  double *G = (double*) malloc(stride[0]*sizeof(double));
  double *H = (double*) malloc(stride[0]*sizeof(double));

  int i,sx=stride[1],sy=stride[2],sz=stride[3];
  intercell_flux_sweep(P,F,1);
  intercell_flux_sweep(P,G,2);
  intercell_flux_sweep(P,H,3);

  for (i=sx; i<stride[0]; ++i)
    {
      L[i] = -(F[i]-F[i-sx])/dx - (G[i]-G[i-sy])/dy - (H[i]-H[i-sz])/dz;
    }

  constrained_transport_3d(F,G,H,stride);

  free(F);
  free(G);
  free(H);

  return 0;
}

/*------------------------------------------------------------------------------
 *
 * Public Functions Definitions
 *
 */
int integrate_init(int N[4], double L[4], int num_comp, int num_dims)
{
  NQ = num_comp;
  ND = num_dims;

  stride[0] = N[1]*N[2]*N[3]*NQ;
  stride[1] =      N[2]*N[3]*NQ;
  stride[2] =           N[3]*NQ;
  stride[3] =                NQ;

  dx = L[1] / (N[1]-2*N[0]);
  dy = L[2] / (N[2]-2*N[0]);
  dz = L[3] / (N[3]-2*N[0]);

  switch (num_dims)
    {
    case 1: drive_sweeps = drive_sweeps_1d; break;
    case 2: drive_sweeps = drive_sweeps_2d; break;
    case 3: drive_sweeps = drive_sweeps_3d; break;
    }

  failure_mask = (int*) malloc(stride[0]/NQ*sizeof(int));

  return 0;
}
int integrate_free()
{
  free(failure_mask);
  return 0;
}
int get_failure_mask(int *M)
{
  memcpy(M, failure_mask, stride[0]/NQ*sizeof(int));
  return 0;
}
int set_riemann_solver(int solver)
{
  switch (solver)
    {
    case 0:
      intercell_flux = hll_flux;
      break;
    case 1:
      intercell_flux = hllc_flux;
      break;
    }
  return 0;
}

int cons_to_prim_array(double *U, double *P, int N)
{
  int i,failures=0;
  for (i=0; i<N*NQ; i+=NQ)
    {
      failure_mask[i/NQ] = 0;
      if (cons_to_prim_point(&U[i],&P[i]))
        {
	  failure_mask[i/NQ] = 1;
	  failures++;
        }
    }
  return failures;
}
int prim_to_cons_array(double *P, double *U, int N)
{
  int i;
  for (i=0; i<N*NQ; i+=NQ)
    {
      prim_to_cons_point(&P[i], &U[i]);
    }
  return 0;
}

int advance_state_fwd_euler(double *P, double dt)
{
  int i,failures=0;
  double *L  = (double*) malloc(stride[0]*sizeof(double));
  double *U0 = (double*) malloc(stride[0]*sizeof(double));

  prim_to_cons_array(P,U0,stride[0]/NQ);

  drive_sweeps(P,L);
  for (i=0; i<stride[0]; ++i) U0[i] += dt*L[i];
  failures += cons_to_prim_array(U0,P,stride[0]/NQ);

  free(L);
  free(U0);

  return failures;
}
int advance_state_midpoint(double *P, double dt)
{
  int i,failures=0;
  double *L  = (double*) malloc(stride[0]*sizeof(double));
  double *U0 = (double*) malloc(stride[0]*sizeof(double));
  double *U1 = (double*) malloc(stride[0]*sizeof(double));

  prim_to_cons_array(P,U0,stride[0]/NQ);

  drive_sweeps(P,L);
  for (i=0; i<stride[0]; ++i) U1[i] = U0[i] + 0.5*dt*L[i];
  failures += cons_to_prim_array(U1,P,stride[0]/NQ);

  drive_sweeps(P,L);
  for (i=0; i<stride[0]; ++i) U0[i] += dt*L[i];
  failures += cons_to_prim_array(U0,P,stride[0]/NQ);

  free(L);
  free(U0);
  free(U1);

  return failures;
}
int advance_state_RK3(double *P, double dt)
{
  int i,failures=0;
  double *L  = (double*) malloc(stride[0]*sizeof(double));
  double *U0 = (double*) malloc(stride[0]*sizeof(double));
  double *U1 = (double*) malloc(stride[0]*sizeof(double));

  prim_to_cons_array(P,U0,stride[0]/NQ);

  drive_sweeps(P,L);
  for (i=0; i<stride[0]; ++i) U1[i] =      U0[i]              +      dt*L[i];
  failures += cons_to_prim_array(U1,P,stride[0]/NQ);

  drive_sweeps(P,L);
  for (i=0; i<stride[0]; ++i) U1[i] = 3./4*U0[i] + 1./4*U1[i] + 1./4*dt*L[i];
  failures += cons_to_prim_array(U1,P,stride[0]/NQ);

  drive_sweeps(P,L);
  for (i=0; i<stride[0]; ++i) U0[i] = 1./3*U0[i] + 2./3*U1[i] + 2./3*dt*L[i];
  failures += cons_to_prim_array(U0,P,stride[0]/NQ);

  free(L);
  free(U0);
  free(U1);

  return failures;
}
int advance_state_ctu_hancock(double *P, double dt)
{
  const int num_dims=ND;
  const int D1=(num_dims>=1);
  const int D2=(num_dims>=2);
  const int D3=(num_dims>=3);

  double *U    = (double*) malloc(stride[0]*sizeof(double));
  double *F    = (double*) malloc(stride[0]*D1*sizeof(double));
  double *G    = (double*) malloc(stride[0]*D2*sizeof(double));
  double *H    = (double*) malloc(stride[0]*D3*sizeof(double));

  double *Ux   = (double*) malloc(stride[0]*D1*sizeof(double));
  double *Uy   = (double*) malloc(stride[0]*D2*sizeof(double));
  double *Uz   = (double*) malloc(stride[0]*D3*sizeof(double));

  double *Px   = (double*) malloc(stride[0]*D1*sizeof(double));
  double *Py   = (double*) malloc(stride[0]*D2*sizeof(double));
  double *Pz   = (double*) malloc(stride[0]*D3*sizeof(double));

  double *dPdx = (double*) malloc(stride[0]*D1*sizeof(double));
  double *dPdy = (double*) malloc(stride[0]*D2*sizeof(double));
  double *dPdz = (double*) malloc(stride[0]*D3*sizeof(double));

  int i,j,sx=stride[1],sy=stride[2],sz=stride[3],failures=0;
  prim_to_cons_array(P,U,stride[0]/NQ);

  /* Step 1
     ---------------------------------------------------------------------------------------
     Compute the slopes of primitive quantities in each direction.

     Note: These derivatives, computed at the beginning of the time step, are used for the
     reconstructed values in the Hancock and Godunov operators for both the predictor and
     corrector steps.
     ---------------------------------------------------------------------------------------
  */
  for (i=sx; i<stride[0]-sx; ++i)
    {
      if(D1) dPdx[i] = slope_limiter(P[i-sx], P[i], P[i+sx]);
      if(D2) dPdy[i] = slope_limiter(P[i-sy], P[i], P[i+sy]);
      if(D3) dPdz[i] = slope_limiter(P[i-sz], P[i], P[i+sz]);
    }

  /* Step 2
     ---------------------------------------------------------------------------------------
     Apply the Godunov operator to each faces in all directions.

     Notes:

     1) The whole array of intercell fluxes is computed at once, rather than one cell at a
     time. This avoids redundant evaluations of the Godunov operator by computing only one
     solution to the Riemann problem per face.

     2) For the predictor step along a given axis, only the intercell fluxes in the
     transverse directions are used. The Hancock operator is applied along the longitudinal
     direction.
     ---------------------------------------------------------------------------------------
  */

  for (i=0; i<stride[0]-sx; i+=NQ)
    {
      double Pl[MAXNQ], Pr[MAXNQ];

      for (j=0; j<NQ*D1; ++j)
        {
          Pr[j] = P[i+sx+j] - 0.5*dPdx[i+sx+j];
          Pl[j] = P[i   +j] + 0.5*dPdx[i   +j];
        }
      if(D1) intercell_flux(Pl, Pr, 0, &F[i], 0.0, 1);

      for (j=0; j<NQ*D2; ++j)
        {
          Pr[j] = P[i+sy+j] - 0.5*dPdy[i+sy+j];
          Pl[j] = P[i   +j] + 0.5*dPdy[i   +j];
        }
      if(D2) intercell_flux(Pl, Pr, 0, &G[i], 0.0, 2);

      for (j=0; j<NQ*D3; ++j)
        {
          Pr[j] = P[i+sz+j] - 0.5*dPdz[i+sz+j];
          Pl[j] = P[i   +j] + 0.5*dPdz[i   +j];
        }
      if(D3) intercell_flux(Pl, Pr, 0, &H[i], 0.0, 3);
    }

  /* Step 3
     ---------------------------------------------------------------------------------------
     Apply the Hancock operators and complete the corrector step.

     Notes: This update occurs for a given cell all at once. The Hancock operator is
     evaluated by summing the fluxes on the inner walls of the local cell, so there is no
     danger of redundant calculations.
     ---------------------------------------------------------------------------------------
  */

  for (i=0; i<stride[0]; i+=NQ)
    {
      double PL[MAXNQ], PR[MAXNQ]; // Capital L/R refers to left and right interior walls of
      double UL[MAXNQ], UR[MAXNQ]; // the local cell, whereas lower-case l/r refer to i_{+-1/2}

      double FL[MAXNQ], FR[MAXNQ];
      double GL[MAXNQ], GR[MAXNQ];
      double HL[MAXNQ], HR[MAXNQ];

      if (D1)
        {
          for (j=0; j<NQ; ++j)
            {
              PL[j] = P[i+j] - 0.5*dPdx[i+j]; // Primitive states on the inner x-facing
              PR[j] = P[i+j] + 0.5*dPdx[i+j]; // walls of the local cell.
            }

          prim_to_cons_point(PL,UL); // Corresponding conserved quantities and fluxes
          prim_to_cons_point(PR,UR); // in the x-direction.

          flux_and_eval(UL, PL, FL, 0, 0, 1);
          flux_and_eval(UR, PR, FR, 0, 0, 1);
        }
      if (D2)
        {
          for (j=0; j<NQ; ++j)
            {
              PL[j] = P[i+j] - 0.5*dPdy[i+j]; // Primitive states on the inner y-facing
              PR[j] = P[i+j] + 0.5*dPdy[i+j]; // walls of the local cell.
            }

          prim_to_cons_point(PL,UL); // Corresponding conserved quantities and fluxes
          prim_to_cons_point(PR,UR); // in the y-direction.

          flux_and_eval(UL, PL, GL, 0, 0, 2);
          flux_and_eval(UR, PR, GR, 0, 0, 2);
        }
      if (D3)
        {
          for (j=0; j<NQ; ++j)
            {
              PL[j] = P[i+j] - 0.5*dPdz[i+j]; // Primitive states on the inner z-facing
              PR[j] = P[i+j] + 0.5*dPdz[i+j]; // walls of the local cell.
            }

          prim_to_cons_point(PL,UL); // Corresponding conserved quantities and fluxes
          prim_to_cons_point(PR,UR); // in the z-direction.

          flux_and_eval(UL, PL, HL, 0, 0, 3);
          flux_and_eval(UR, PR, HR, 0, 0, 3);
        }

      switch (num_dims)
        {
        case 1:                    /**** 1D ****/
          for (j=0; j<NQ; ++j)
            {
              //                   Hancock (normal)
              // ==========================================
              Ux[i+j] = U[i+j] - ((FR[j]-FL[j])/dx)*0.5*dt;
              // ==========================================
            }
          break;
        case 2:                    /**** 2D ****/
          for (j=0; j<NQ; ++j)
            {
              //                   Hancock (normal)            Godunov (transverse)
              // ==================================================================
              Ux[i+j] = U[i+j] - ((FR[j]-FL[j])/dx + (G[i+j]-G[i-sy+j])/dy)*0.5*dt;
              Uy[i+j] = U[i+j] - ((GR[j]-GL[j])/dy + (F[i+j]-F[i-sx+j])/dx)*0.5*dt;
              // ==================================================================
            }
          break;
        case 3:                    /**** 3D ****/
          for (j=0; j<NQ; ++j)
            {
              //                   Hancock (normal)            Godunov (transverse)
              // ==========================================================================================
              Ux[i+j] = U[i+j] - ((FR[j]-FL[j])/dx + (G[i+j]-G[i-sy+j])/dy + (H[i+j]-H[i-sz+j])/dz)*0.5*dt;
              Uy[i+j] = U[i+j] - ((GR[j]-GL[j])/dy + (H[i+j]-H[i-sz+j])/dz + (F[i+j]-F[i-sx+j])/dx)*0.5*dt;
              Uz[i+j] = U[i+j] - ((HR[j]-HL[j])/dz + (F[i+j]-F[i-sx+j])/dx + (G[i+j]-G[i-sy+j])/dy)*0.5*dt;
              // ==========================================================================================
            }
          break;
        }
    }

  if(D1) memcpy(Px,P,stride[0]*sizeof(double));
  if(D2) memcpy(Py,P,stride[0]*sizeof(double));
  if(D3) memcpy(Pz,P,stride[0]*sizeof(double));

  if(D1) failures += cons_to_prim_array(Ux,Px,stride[0]/NQ); // Inversion of the predicted U-states
  if(D2) failures += cons_to_prim_array(Uy,Py,stride[0]/NQ); // for each direction.
  if(D3) failures += cons_to_prim_array(Uz,Pz,stride[0]/NQ);

  /* Step 4
     ---------------------------------------------------------------------------------------
     Complete the integration by applying the Godunov fluxes.

     Notes: The final derivative operator is obtained from the Godunov intercell fluxes,
     which in turn are computed using the cell-centered predicted state, reconstructed to
     zone edges using the derivatives from the beginning of the time step.
     ---------------------------------------------------------------------------------------
  */
  for (i=sx; i<stride[0]-sx; i+=NQ)
    {
      double Pl[MAXNQ], Pr[MAXNQ];

      for (j=0; j<NQ*D1; ++j)
        {
          Pr[j] = Px[i+sx+j] - 0.5*dPdx[i+sx+j];
          Pl[j] = Px[i   +j] + 0.5*dPdx[i   +j];
        }
      if(D1) intercell_flux(Pl, Pr, 0, &F[i], 0.0, 1);

      for (j=0; j<NQ*D2; ++j)
        {
          Pr[j] = Py[i+sy+j] - 0.5*dPdy[i+sy+j];
          Pl[j] = Py[i   +j] + 0.5*dPdy[i   +j];
        }
      if(D2) intercell_flux(Pl, Pr, 0, &G[i], 0.0, 2);

      for (j=0; j<NQ*D3; ++j)
        {
          Pr[j] = Pz[i+sz+j] - 0.5*dPdz[i+sz+j];
          Pl[j] = Pz[i   +j] + 0.5*dPdz[i   +j];
        }
      if(D3) intercell_flux(Pl, Pr, 0, &H[i], 0.0, 3);
    }

  switch (num_dims)
    {
      /*---------------------------------- 1D --------------------------------*/
    case 1:
      // No constrained transport
      for (i=sx; i<stride[0]; ++i)
        {
          U[i] -= dt*((F[i]-F[i-sx])/dx);
        }
      break;
      /*---------------------------------- 2D --------------------------------*/
    case 2:
      constrained_transport_2d(F,G,stride);
      for (i=sx; i<stride[0]; ++i)
        {
          U[i] -= dt*((F[i]-F[i-sx])/dx + (G[i]-G[i-sy])/dy);
        }
      break;
      /*---------------------------------- 3D --------------------------------*/
    case 3:
      constrained_transport_3d(F,G,H,stride);
      for (i=sx; i<stride[0]; ++i)
        {
          U[i] -= dt*((F[i]-F[i-sx])/dx + (G[i]-G[i-sy])/dy + (H[i]-H[i-sz])/dz);
        }
      break;
    }
  failures += cons_to_prim_array(U,P,stride[0]/NQ);

  free(U);

  if(D1) free(F);
  if(D2) free(G);
  if(D3) free(H);

  if(D1) { free(Ux); free(Px); free(dPdx); }
  if(D2) { free(Uy); free(Py); free(dPdy); }
  if(D3) { free(Uz); free(Pz); free(dPdz); }

  return failures;
}
