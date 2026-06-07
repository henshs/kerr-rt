#include <gsl/gsl_poly.h>
#include <cmath>
#include <cstdio>



double kerr_rms(double a, double sgn);

int main(int c, char *argv[])
{
    double a,sgn,r;

    if(c!=3){
        fprintf(stderr,"kerr_rms <a> <sgn>\n");
        return -1;
    }

    a=atof(argv[1]);
    sgn=atof(argv[2]);

    printf("%lf",kerr_rms(a,sgn));
}


double kerr_rms(double a, double sgn)
{
    double a2,Z1,Z2;

    a2=a*a;

    if(a<=1) Z1=1 + pow(1-a2,1.0/3.0)*(pow(1+a,1./3.) + pow(1-a,1./3.));
    else Z1=1 - pow(a2-1,1.0/3.0)*(pow(1+a,1./3.) - pow(a-1,1./3.));
    Z2=sqrt(3*a2+Z1*Z1);

    return 3 + Z2 - sgn*sqrt((3-Z1)*(3+Z1+2*Z2));
}