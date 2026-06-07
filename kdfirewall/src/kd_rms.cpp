#include "kd_rms.h"

double kd_rms(double a)
{
    double Z1,Z2,a2;
    double w1,w2;
    
    a2=a*a;

    w1=1-a;
    if(w1<0) w1=-pow(-w1,1.0/3.0);
    else w1=pow(w1,1.0/3.0);
    w2=1-a2;
    if(w2<0) w2=-pow(-w2,1.0/3.0);
    else w2=pow(w2,1.0/3.0);

    Z1=1+w2*(pow(1+a,1.0/3.0)+w1);
    Z2=sqrt(3*a2+Z1*Z1);

    if(a>=0) return 3+Z2-sqrt((3-Z1)*(3+Z1+2*Z2));
    
    return 3 + Z2 + sqrt((3-Z1)*(3+Z1+2*Z2));
}