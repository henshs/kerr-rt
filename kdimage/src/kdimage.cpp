#include "kdimage.h"

void kdimage_init(FILE *f, KDIMAGE *kd)
{
    double alphamin,alphamax,betamin,betamax;
    int nalpha,nbeta;
    double alpha,beta,flux,g,re;
    int ialpha,ibeta,reflex;

    fscanf(f,"%d %d %lf %lf %lf %lf",&nalpha,&nbeta,&alphamin,&alphamax,&betamin,&betamax);
    kd->alphamin=alphamin;
    kd->alphamax=alphamax;
    kd->betamin=betamin;
    kd->betamax=betamax;
    kd->nalpha=nalpha;
    kd->nbeta=nbeta;

    kd->fluxtab = new double[nalpha*nbeta];

   
    for(int i=0;i<nalpha*nbeta;i++) kd->fluxtab[i]=0;

    while(1)
    {
        int sts=fscanf(f,"%d %d %lf %lf %lf %lf %lf",&ialpha,&ibeta,&alpha,&beta,&g,&re,&flux);
        if(sts==EOF) break;

        int indx = ibeta*kd->nalpha+ialpha;
	
        kd->fluxtab[indx] = flux;
       
    
    }
}

void kdimage_print(KDIMAGE *kd)
{
    double alpha,beta,dalpha,dbeta;

    dalpha=(kd->alphamax-kd->alphamin)/kd->nalpha;
    dbeta=(kd->betamax-kd->betamin)/kd->nbeta;

    for(int ialpha=0;ialpha<kd->nalpha;ialpha++)
    {
        alpha=kd->alphamin+ialpha*dalpha;
        for(int ibeta=0;ibeta<kd->nbeta;ibeta++)
        {
            beta=kd->betamin+ibeta*dbeta;
            int indx=ibeta*kd->nalpha+ialpha;
            printf("%lf %lf %lf\n",alpha,beta,kd->fluxtab[indx]);
        }
    }
}

void kdimage_free(KDIMAGE *kd)
{
    delete kd->fluxtab;
}
