#ifndef __KDIMAGE_H__
#define __KDIMAGE_H__

#include <cstdio>
#include <cstdlib>
#include <cmath>
using namespace std;

struct KDIMAGE
{
    int nalpha,nbeta;
    double alphamin,alphamax,betamin,betamax;
    double *fluxtab;
    double *gtab;
    double *alphatab;
    double *betatab;
    
};

extern void kdimage_init(FILE *f, KDIMAGE *kd);
extern void kdimage_print(KDIMAGE *kd);
extern void kdimage_free(KDIMAGE *kd);


#endif
