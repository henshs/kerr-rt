#include "numeric.h"
#include "kerr.h"
#include "kdisk.h"
#include "kd_flux.h"
#include "kd_rms.h"
#include "kd_photon.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace std;
int  get_color(double r, double rmin,double rmax,int n);
void print_geodesics(ODESOLUTION *sol);

int main(int argc, char *argv[])
{

    gsl_set_error_handler_off();
    vector<KD_PHOTON> photon_list;
    KD_PHOTON photon;
    KD_HEADER kd_header;
    KDFLUX kd;

    double k1[4];
    double k2[4];
    double re[16];
    double opacity[16];
    int recnt;

    ODESOLUTION sol;
    KERRPARAMS p;
    double alpha,beta,
           alphamin,alphamax,
           betamin,betamax,
           dalpha,dbeta, wmax,
           rmin,rmax,rms;

    int nalpha,nbeta;

    if(argc!=2){
        fprintf(stderr,"params: <init.tab>\n");
        return -1;
    }

    FILE *fin=fopen(argv[1],"rt");
    char s[16];
    fscanf(fin,"%s %lf",s,&(p.a));
    fscanf(fin,"%s %lf %lf",s,&(p.uo),&(p.mo));p.mo=cos(p.mo*M_PI/180.0);
    fscanf(fin,"%s %lf %lf %lf %lf",s,&alphamin,&alphamax,&betamin,&betamax);
    fscanf(fin,"%s %d %d",s,&nalpha,&nbeta);
    fscanf(fin,"%s %lf %lf",s,&rmin,&rmax);//Keplerian disk inner and outer edges
    fclose(fin);


    rms = kd_rms(p.a); 
    fprintf(stderr,"rms=%lf\n",rms);

    kerr_get_uh(&p);

    kd_flux_interpolate(rmin,rmax,100,p.a,6.5*pow(10.0,9.0),2.0*pow(10.0,-5.0),&kd); 
    fprintf(stderr,"rmin=%lf rmax=%lf\n",rmin,rmax);

    dalpha = (alphamax-alphamin)/nalpha;
    dbeta = (betamax - betamin)/nbeta;

    for(int i=0;i<16;i++) opacity[i]=(double)(1)/(1<<i);
    odesol_init(4,2048,&sol);

    for(int i=0;i<nalpha;i++)
    {
        alpha=alphamin + i*dalpha;
        for(int j=0;j<nbeta;j++)
        {
            beta=betamin+j*dbeta;
            
            kerr_get_lq(alpha,beta,&p);

            kerr_get_u_turpoi(&p);
            //double re;
            int istat=0;
            recnt = 0;

            if(!(p.u_turpoi==1 && 1/p.ut>1.5*rmax))
            {

                //fprintf(stderr,"after turpoi if");
            
                if(p.u_turpoi==1) wmax=kerr_get_wmax(p.uo,1/(1.5*rmax),&p);
                else wmax=kerr_get_wmax(p.uo,0.9999*p.uh,&p);

                
                double y0[4];
                p.usgn=1.0;
                kerr_photon_ic(y0,&p);

                
                ndsolve(kerr_odesys,&p,y0,0,wmax,1e-12,&sol);
                
                //fprintf(stderr,"after ndsolve");
                
               // print_geodesics(&sol);
                istat=kdisk_intersection(rmin,rmax,&sol,&p,re,&recnt);
                //if(recnt>0) fprintf(stderr,"recnt=%d\n",recnt);
                
            }
            /* OLD VERSION OF STORING RESULTS*/
            
            double g=0;
	        double reok=0;
            int color=0;
            double I0=0;
            int ni;

	        if(istat==1 && re[0]>rmin && re[0]<rmax )
            {
                
                  g=kdisk_freq_shift(re[0],&p);
		          reok=re[0];
                  I0=g*g*g*g*kd_flux(reok,&kd); //(1+0.5*RAND())
                
                  //fprintf(stderr,"after I0");
            }

            printf("%d %d %lf %lf %lf %lf %lf\n",i,j,alpha,beta,g,reok,I0);
            
        }
        printf("\n");

        fprintf(stderr,"Finished %.2lf %% \r",100.0*i/(nalpha-1));
    }
    fprintf(stderr,"\n");

    odesol_free(&sol);

    kd_flux_spline_free(&kd);

    return 0;
}


int  get_color(double r, double rmin,double rmax,int n)
{
    double dr=(rmax-rmin)/n;
    return (int)((r-rmin)/dr);
}

void print_geodesics(ODESOLUTION *sol)
{
    int n=sol->nactual;

    for(int i=0;i<n;i++)
    {
        double r,th;

        r=sol->ytab[2][i];
        r=1/r;
        th=acos(sol->ytab[3][i]);

        printf("%.10lf %.10lf %.10lf %.10lf\n",r*sin(th),r*cos(th),r,th*180.0/M_PI);
    }
}
