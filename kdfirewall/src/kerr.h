#ifndef __KERR_H__
#define __KERR_H__

#include "numeric.h"

struct KERRPARAMS
{
    double a,uh;
    double l,q;
    double ut;
    int u_turpoi;
    double uo,mo;
    double msgn,usgn;
};

extern void kerr_get_uh(KERRPARAMS *p);
extern void kerr_get_u_turpoi(KERRPARAMS *p);
extern void kerr_get_lq(double alpha,double beta,KERRPARAMS *p);
extern double kerr_get_wmax(double uo, double ue,KERRPARAMS *p);
extern int kerr_odesys(double w, const double y[], double f[], void *par);
extern void kerr_photon_ic(double y[],KERRPARAMS *p);


#endif