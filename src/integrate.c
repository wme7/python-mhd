
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
int integrate_init(int N[4], double L[4], int num_dims);
int advance_state_fwd_euler  (double *P, double dt);
int advance_state_midpoint   (double *P, double dt);
int advance_state_RK3        (double *P, double dt);


/*------------------------------------------------------------------------------
 *
 * External Dependecies
 *
 */
int hll_flux(const double *pl, const double *pr, double *U, double *F,
	     double s, int dim);

int prim_to_cons_array(const double *P, double *U, int N);
int cons_to_prim_array(const double *U, double *P, int N);

int constrained_transport_2d(double *Fx, double *Fy);
int constrained_transport_3d(double *Fx, double *Fy, double *Fz);


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


/*------------------------------------------------------------------------------
 *
 * Private Data
 *
 */
static int stride[4];
static double plm_theta=2.0;
static double dx,dy,dz;

static double (*slope_limiter)(double, double, double);
static int (*godunov_intercell_flux)(const double*, const double*, double*, double*,
				     double, int);
static int (*drive_sweeps)(const double*, double*);


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
/*------------------------------------------------------------------------------*/



/*------------------------------------------------------------------------------
 *
 * Public Functions Definitions
 *
 */
int integrate_init(int N[4], double L[4], int num_dims)
{
  stride[0] = N[1]*N[2]*N[3]*NQ;
  stride[1] =      N[2]*N[3]*NQ;
  stride[2] =           N[3]*NQ;
  stride[3] =                NQ;

  dx = L[1] / (N[1]-2*N[0]);
  dy = L[2] / (N[2]-2*N[0]);
  dz = L[3] / (N[3]-2*N[0]);

  slope_limiter = plm_minmod;
  godunov_intercell_flux = hll_flux;

  switch (num_dims)
    {
    case 1: drive_sweeps = drive_sweeps_1d; break;
    case 2: drive_sweeps = drive_sweeps_2d; break;
    case 3: drive_sweeps = drive_sweeps_3d; break;
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

  return 0;
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

  return 0;
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

  return 0;
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
  double Pl[NQ], Pr[NQ];

  for (i=0; i<S; ++i)
    {
      F[i] = 0.0;
    }
  for (i=S; i<stride[0]-S*2; i+=NQ)
    {
      reconstruct_plm(&P[i], Pl, Pr, S);
      godunov_intercell_flux(Pl, Pr, 0, &F[i], 0.0, dim);
    }
  for (i=stride[0]-S*2; i<stride[0]; ++i)
    {
      F[i] = 0.0;
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

  constrained_transport_2d(F,G);

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

  constrained_transport_3d(F,G,H);

  free(F);
  free(G);
  free(H);

  return 0;
}
