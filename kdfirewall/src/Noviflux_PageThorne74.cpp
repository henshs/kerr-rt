#include "Noviflux_PageThorne74.h"

double func_flux(double rad,double spin,double M,double Mdot)
{
    double scriptLbyL0,multi_factor,y,F;

    y=sqrt(rad);

    scriptLbyL0 = func_scriptLbyL0(y,spin);

    multi_factor = func_factor(y,spin);

    double k=9.1721376255e+28*(Mdot/M);

    F = k*(1.0/(4.0*pi*rad))*scriptLbyL0*multi_factor; // in erg/(cm^2.sec)

    return F;

}

double func_rms(double spin)
{
    double spin2,Z1,Z2,rms;

    spin2=pow(spin,2.0);

    Z1=1+pow(1-spin2,1.0/3.0)*(pow(1+spin,1.0/3.0)+pow(1-spin,1.0/3.0));
    Z2=sqrt(3*spin2+pow(Z1,2.0));

    rms=3+ Z2 - sqrt((3-Z1)*(3+Z1+2*Z2));

    return rms;

}

double func_factor(double y, double spin)
{
    double y2,y3;

    y2=y*y;
    y3=y*y*y;

    return 1.5/(y2*(y3-3.0*y+2.0*spin));
}

double func_y0(double spin)
{
    return sqrt(func_rms(spin));
}

double func_y1(double spin)
{

    return +2*cos((1.0/3.0)*acos(spin)-pi/3.0);
}

double func_y2(double spin)
{

    return +2*cos((1.0/3.0)*acos(spin)+pi/3.0);
}

double func_y3(double spin)
{

    return -2*cos((1.0/3.0)*acos(spin));
}

double func_f0(double y, double spin)
{
    double y0=func_y0(spin);

    return y-y0-1.5*spin*log(y/y0);
}

double func_f1(double y, double spin)
{
    double y0,y1,y2,y3;
    y0=func_y0(spin);
    y1=func_y1(spin);
    y2=func_y2(spin);
    y3=func_y3(spin);

    return ((3.0*pow((y1-spin),2.0))/(y1*(y1-y2)*(y1-y3)))*log((y-y1)/(y0-y1));
}

double func_f2(double y, double spin)
{
    double y0,y1,y2,y3;
    y0=func_y0(spin);
    y1=func_y1(spin);
    y2=func_y2(spin);
    y3=func_y3(spin);

    return ((3.0*pow((y2-spin),2.0))/(y2*(y2-y1)*(y2-y3)))*log((y-y2)/(y0-y2));
}

double func_f3(double y, double spin)
{
    double y0,y1,y2,y3;
    y0=func_y0(spin);
    y1=func_y1(spin);
    y2=func_y2(spin);
    y3=func_y3(spin);

    return ((3.0*pow((y3-spin),2.0))/(y3*(y3-y1)*(y3-y2)))*log((y-y3)/(y0-y3));
}

double func_scriptLbyL0(double y, double spin)
{
    double f0,f1,f2,f3;

    f0=func_f0(y,spin);
    f1=func_f1(y,spin);
    f2=func_f2(y,spin);
    f3=func_f3(y,spin);

    return (f0-f1-f2-f3);
}
