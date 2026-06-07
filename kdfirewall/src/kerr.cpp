#include "kerr.h"

void kerr_get_uh(KERRPARAMS *p)
{
    double a,a2;
    a=p->a;
    a2=a*a;

    p->uh = (a<=1.0 ? 1/(1+sqrt(1-a2)) : 1e3);
}

void kerr_get_u_turpoi(KERRPARAMS *p)
{
    gsl_poly_complex_workspace *ws = gsl_poly_complex_workspace_alloc(5);
    double b[5];
    double z[8];
    double a,a2,l,l2,q;
    int n;

    a=p->a;a2=a*a;
    l=p->l;l2=l*l;
    q=p->q;

    double uh = a<=1.0 ? 1/(1+sqrt(1-a2)) : 1e3;


    if(fabs(q)>1e-8)
    {
        n=4;
        b[0]=1.0;
        b[1]=0.0;
        b[2]=a2-l2-q;
        b[3]=2*((a-l)*(a-l)+q);
        b[4]=-a2*q;

        gsl_poly_complex_solve(b,5,ws,z);

        gsl_poly_complex_workspace_free(ws);

    }else{
        b[0]=1.0/(2*((a-l)*(a-l)));
        b[1]=0.0;
        b[2]=(a2-l2)/(2*((a-l)*(a-l)));
        b[3]=1.0;
        gsl_complex z1,z2,z3;
        n=gsl_poly_complex_solve_cubic(b[2],b[1],b[0],&z1,&z2,&z3);
        z[0]=GSL_REAL(z1);z[1]=GSL_IMAG(z1);
        z[2]=GSL_REAL(z2);z[3]=GSL_IMAG(z2);
        z[4]=GSL_REAL(z3);z[5]=GSL_IMAG(z3);
    }
    
    double ut=uh;
    int sts=0;

    for(int i=0;i<n;i++)
    {
        if(fabs(z[2*i+1])<1e-8 && z[2*i]>0 && z[2*i]<ut){
            sts=1;
            ut=z[2*i];
        }
    } 

    p->ut=ut;
    p->u_turpoi = sts;
}

void kerr_get_lq(double alpha,double beta,KERRPARAMS *p)
{
    double a,a2,mo2;

    a=p->a;a2=a*a;
    mo2=p->mo*p->mo;
    p->l=-alpha*sqrt(1-mo2);
    p->q=beta*beta + (alpha*alpha-a2)*mo2;
    p->msgn=GSL_SIGN(beta);
}

double funcU(double u, void *par){
    KERRPARAMS *p=(KERRPARAMS *)par;
    double a,a2,l,l2,q,U,u2;

    a=p->a;a2=a*a;
    l=p->l;l2=l*l;
    q=p->q;

    u2=u*u;

    U=1 + u2*(a2-l2-q + u*( 2*((a-l)*(a-l)+q) - a2*q*u));

    return 1/sqrt(fabs(U));
}

double kerr_get_wmax(double uo, double ue,KERRPARAMS *p)
{
    const int n=1000;
    const double eps=1e-8;
    gsl_integration_workspace *ws = gsl_integration_workspace_alloc(n);
    gsl_function F;
    double res,err,tmp;

    F.function=funcU;
    F.params = p;

    if(p->u_turpoi==0){
        gsl_integration_qags(&F,uo,ue,0,eps,n,ws,&res,&err);
    }else{
        double ut=0.999999999*p->ut;
        gsl_integration_qags(&F,uo,ut,0,eps,n,ws,&res,&err);
        if(fabs(p->ut-ue)>eps) gsl_integration_qags(&F,ue,ut,0,eps,n,ws,&tmp,&err);
        else tmp=0;
        res=res+tmp;
    }   

    gsl_integration_workspace_free(ws);

    return res;
}

int kerr_odesys(double w, const double y[], double f[], void *par)
{
    KERRPARAMS *p=(KERRPARAMS *)par;
    double l,l2,q,a,a2,dU,dM;
    double u,u2,m,m2,ku,km;

    a=p->a;a2=a*a;
    l=p->l;l2=l*l;
    q=p->q;

    ku=y[0];
    km=y[1];
    u=y[2];u2=u*u;
    m=y[3];m2=m*m;



    dU=2*u*(a2-l2-q + u*( 3*((a-l)*(a-l)+q) - 2*a2*q*u ));
    dM=2*m*(a2-l2-q - 2*a2*m2);

    f[0]=0.5*dU;
    f[1]=0.5*dM;
    f[2]=ku;
    f[3]=km;

    return GSL_SUCCESS;
}

void kerr_photon_ic(double y[],KERRPARAMS *p)
{   
    double uo,uo2,mo,mo2,l,l2,a,a2,q;

    uo=p->uo;uo2=uo*uo;
    mo=p->mo;mo2=mo*mo;
    l=p->l;l2=l*l;
    q=p->q;
    a=p->a;a2=a*a;
    
    y[0]=p->usgn*sqrt(fabs( 1.0 + uo2*(a2-l2-q+ 2*((a-l)*(a-l)+q)*uo - a2*a*uo2) ));
    y[1]=p->msgn*sqrt(fabs(p->q + mo2*(a2-l2-q - a2*mo2)));
    y[2]=p->uo;
    y[3]=p->mo;
}