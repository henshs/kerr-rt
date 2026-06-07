#ifndef __KD_FLUX_H__
#define __KD_FLUX_H__
#include "kd_rms.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <cmath>

struct KDFLUX{
    gsl_interp_accel *acc;
    gsl_spline *spline;
    double rin,rout;
};

void kd_flux_interpolate(double rin, double rout, int n, double a, double M,double Mdot,KDFLUX *kd);

void kd_flux_spline_free(KDFLUX *kd);

double kd_flux_spline_val(double r, KDFLUX *kd);


double kd_flux(double re,KDFLUX *kd);

#endif