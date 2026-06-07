#include "kd_flux.h"
#include "Noviflux_PageThorne74.h"
#include <gsl/gsl_errno.h>

void kd_flux_interpolate(double rin, double rout, int n, double a, double M,double Mdot,KDFLUX *kd)
{
    kd->rin=rin;
    kd->rout=rout;
    double *rtab,
           *fluxtab;
    double r,dr;

    rtab = new double[n];
    fluxtab = new double[n];
    
    dr=(1.1*rout-rin)/n;
    for(int i=0;i<n;i++)
    {
        r=rin+i*dr;
        fluxtab[i] = func_flux(r,a,M,Mdot);
        rtab[i]=r;
    }

    kd->acc = gsl_interp_accel_alloc();
    kd->spline = gsl_spline_alloc(gsl_interp_cspline,n);
    gsl_spline_init(kd->spline,rtab,fluxtab,n);

    delete rtab;
    delete fluxtab;
}

void kd_flux_spline_free(KDFLUX *kd)
{
    gsl_spline_free(kd->spline);
    gsl_interp_accel_free(kd->acc);
}


double kd_flux_spline_val(double r, KDFLUX *kd)
{
    return gsl_spline_eval(kd->spline,r,kd->acc);
}


double kd_flux(double re,KDFLUX *kd)
{

    gsl_set_error_handler_off();     

    double flux=kd_flux_spline_val(re,kd);
  
    return flux;
}
  
