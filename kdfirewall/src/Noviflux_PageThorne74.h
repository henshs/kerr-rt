#ifndef __NOVIFLUXABRAMOWICZ_H__
#define __NOVIFLUXABRAMOWICZ_H__

#include <cmath>
using namespace std;

#define pi M_PI


extern double func_flux(double rad,double spin,double M,double Mdot);
extern double func_rms(double spin);
extern double func_factor(double y, double spin);
extern double func_y0(double spin);
extern double func_y1(double spin);
extern double func_y2(double spin);
extern double func_y3(double spin);
extern double func_f0(double y, double spin);
extern double func_f1(double y, double spin);
extern double func_f2(double y, double spin);
extern double func_f3(double y, double spin);
extern double func_scriptLbyL0(double y, double spin);

#endif
