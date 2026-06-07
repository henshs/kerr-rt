#ifndef __KDISK_H__
#define __KDISK_H__

#include "numeric.h"
#include "kerr.h"
#include <cmath>
#include <gsl/gsl_spline.h>

extern double kdisk_freq_shift(double r, KERRPARAMS *p);
extern int kdisk_intersection(double rin, double rout,ODESOLUTION *sol, KERRPARAMS *p,double re[],int *cnt);

#endif
