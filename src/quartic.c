
/*------------------------------------------------------------------------------
 * FILE: quartic.c
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *         Adapted from the nvwa code by Weiqun Zhang
 *
 *------------------------------------------------------------------------------
 */

#include <math.h>

/*------------------------------------------------------------------------------
 *
 * Public Interface
 *
 */
int quartic_equation_init(double d4, double d3, double d2, double d1, double d0);
int quartic_equation_solve_apprx(double *x);
int quartic_equation_solve_exact(double *r1, double *r2, double *r3, double *r4,
				 int *nr12, int *nr34);

/*------------------------------------------------------------------------------
 *
 * Private Data
 *
 */
static double d4,d3,d2,d1,d0;

/*------------------------------------------------------------------------------
 *
 * Internal Dependecies
 *
 */
static int solve_cubic_equation(double  c3, double  c2,  double c1, double c0,
				double *x1, double *x2, double *x3);

/*------------------------------------------------------------------------------
 *
 * Inline functions
 *
 */
static inline double evaluate_f(double x)
{
  double x2 = x*x;
  double x3 = x*x2;
  double x4 = x*x3;
  return d4*x4 + d3*x3 + d2*x2 + d1*x + d0;
}
static inline double evaluate_dfdx(double x)
{
  double x2 = x*x;
  double x3 = x*x2;
  return 4*d4*x3 + 3*d3*x2 + 2*d2*x + d1;
}

/*------------------------------------------------------------------------------
 *
 * Public Functions Definitions
 *
 */
int quartic_equation_init(double D4, double D3, double D2, double D1, double D0)
{
  d4=D4; d3=D3; d2=D2; d1=D1; d0=D0;
  return 0;
}

int quartic_equation_solve_apprx(double *x)
{
  static const double ERROR_TOLR = 1e-14;
  static const int NEWTON_MAX_ITER = 10;

  double f,g;
  int n_iter=0;

  f = evaluate_f(*x);
  while (fabs(f) > ERROR_TOLR)
    {
      f = evaluate_f(*x);
      g = evaluate_dfdx(*x);
      *x -= f/g;
      if (++n_iter > NEWTON_MAX_ITER) return 1;
    }

  return 0;
}
int quartic_equation_solve_exact(double *r1, double *r2, double *r3, double *r4,
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


/*------------------------------------------------------------------------------
 *
 * Private Functions Definitions
 *
 */
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
