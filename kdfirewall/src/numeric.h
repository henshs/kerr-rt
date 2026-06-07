#ifndef __NUMERICS_H__
#define __NUMERICS_H__

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_spline.h>

#define ODESYS(f) int (*f)(double, const double[], double [], void *)
#define FUNC(f) double (*f)(double, void *p)

typedef struct 
{
    int dim;
    int n, nactual;
    double *wtab;
    double **ytab;
}ODESOLUTION;

extern void odesol_init(int dim, int n, ODESOLUTION *sol);
extern void odesol_free(ODESOLUTION *sol);
extern int ndsolve(ODESYS(f), void *p, double y0[],double w1, double w2, double eps, ODESOLUTION *sol);
extern int findroot(FUNC(f),void *p, double w1, double w2, double eps, double *r);

#endif