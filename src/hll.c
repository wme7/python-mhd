
/*------------------------------------------------------------------------------
 * FILE: hll.c
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
void reset_max_lambda();
double get_max_lambda();
int hll_flux(const double *pl, const double *pr, double *U, double *F,
	     double s, int dim);

/*------------------------------------------------------------------------------
 *
 * External Dependecies
 *
 */
int flux_and_eval(const double *U, const double *P, double *F,
		  double *ap, double *am, int dim);
int prim_to_cons_point(const double *P, double *U);
int cons_to_prim_point(const double *U, double *P);


/*------------------------------------------------------------------------------
 *
 * Private Data
 *
 */
static double max_lambda=0.0;


/*------------------------------------------------------------------------------
 *
 * Public Functions Definitions
 *
 */
double get_max_lambda()
{
  return max_lambda;
}
void reset_max_lambda()
{
  max_lambda = 0.0;
}

int hll_flux(const double *pl, const double *pr, double *U, double *F,
	     double s, int dim)
{
  int i;
  double epl, epr, eml, emr;
  double Ul[NQ], Ur[NQ];
  double Pl[NQ], Pr[NQ];
  double Fl[NQ], Fr[NQ];

  memcpy(Pl,pl,NQ*sizeof(double));
  memcpy(Pr,pr,NQ*sizeof(double));

  prim_to_cons_point(Pl,Ul);
  prim_to_cons_point(Pr,Ur);

  flux_and_eval(Ul, Pl, Fl, &epl, &eml, dim);
  flux_and_eval(Ur, Pr, Fr, &epr, &emr, dim);

  double ap = (epl>epr) ? epl : epr;
  double am = (eml<emr) ? eml : emr;

  double ml = (fabs(am)<fabs(ap)) ? fabs(ap) : fabs(am);
  if (max_lambda < ml) max_lambda = ml;

  double F_hll[NQ], U_hll[NQ];
  for (i=0; i<NQ; ++i)
    {
      U_hll[i] = (ap*Ur[i] - am*Ul[i] +       (Fl[i] - Fr[i])) / (ap - am);
      F_hll[i] = (ap*Fl[i] - am*Fr[i] + ap*am*(Ur[i] - Ul[i])) / (ap - am);
    }

  if (U != 0) {
    if      (         s<=am ) for (i=0; i<NQ; ++i) U[i] = Ul   [i];
    else if ( am<s && s<=ap ) for (i=0; i<NQ; ++i) U[i] = U_hll[i];
    else if ( ap<s          ) for (i=0; i<NQ; ++i) U[i] = Ur   [i];
  }
  {
    if      (         s<=am ) for (i=0; i<NQ; ++i) F[i] = Fl   [i];
    else if ( am<s && s<=ap ) for (i=0; i<NQ; ++i) F[i] = F_hll[i];
    else if ( ap<s          ) for (i=0; i<NQ; ++i) F[i] = Fr   [i];
  }
  return 0;
}
