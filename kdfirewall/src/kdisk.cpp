#include "kdisk.h"

double kdisk_freq_shift(double r, KERRPARAMS *p)
{
    double r2,a,a2,A;
    double delta,omega,omega2;
    double l,f;

    r2=r*r;
    a=p->a;
    a2=a*a;
    l=p->l;

    delta = r*(r-2)+a2;
    omega=1/(sqrt(r2*r)+a);
    A=(r2+a2)*(r2+a2)-a2*delta;

    f = 1 - 2*(1-a*omega)*(1-a*omega)/r - (r2+a2)*omega*omega;

    return sqrt(f)/(1-l*omega);
}

int  kdisk_intersection(double rin, double rout,ODESOLUTION *sol, KERRPARAMS *p, double re[],int *cnt)
{
    gsl_spline *uspline,
               *mspline;
    gsl_interp_accel *acc;

    acc = gsl_interp_accel_alloc();

    uspline = gsl_spline_alloc(gsl_interp_cspline,sol->nactual);
    mspline = gsl_spline_alloc(gsl_interp_cspline,sol->nactual);

    gsl_spline_init(uspline,sol->wtab,sol->ytab[2],sol->nactual);
    gsl_spline_init(mspline,sol->wtab,sol->ytab[3],sol->nactual);
    
    double w,w1,w2,dw,wa,wb,fa,fb,n;
    double wm,fm;
    n=100;
    w1=0;
    w2=sol->wtab[sol->nactual-1];
    dw=(w2-w1)/n;

    int istat=0;
    *cnt=0;
    for(int i=0;i<n-1;i++)
    {
        wa = w1 + i*dw;
        wb = wa + dw;

        fa = gsl_spline_eval(mspline,wa,acc);
        fb = gsl_spline_eval(mspline,wb,acc);
        if(fa*fb<0){
            double d=dw;
            while(d>1e-8){
                wm=0.5*(wa+wb);
                fm=gsl_spline_eval(mspline,wm,acc);
                if(fa*fm<0){
                    wb=wm;
                    fb=fm;
                }else{
                    wa=wm;
                    fa=fm;
                }

                d=fabs(wb-wa);
            }

            double rtmp;
            rtmp=gsl_spline_eval(uspline,wm,acc);
            rtmp=1/(rtmp);
            if(rtmp>rin && rtmp<rout){
                istat=1;
                re[*cnt]=rtmp;
                (*cnt)++;    
                if(*cnt==4) break;
            }
        }

    }

    gsl_spline_free(uspline);
    gsl_spline_free(mspline);
    gsl_interp_accel_free(acc);

    return istat;

}
