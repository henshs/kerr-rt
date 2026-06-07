#include "numeric.h"

void odesol_init(int dim, int n, ODESOLUTION *sol)
{
    sol->dim = dim;
    sol->n = n;
    sol->wtab = new double[n];
    sol->ytab = new double*[dim];
    for(int i=0;i<dim;i++) sol->ytab[i] = new double[n];
}

void odesol_free(ODESOLUTION *sol)
{
    delete sol->wtab;
    for(int i=0;i<sol->dim;i++) delete sol->ytab[i];
    delete sol->ytab;
}

int ndsolve(ODESYS(f), void *p, double y0[],double w1, double w2, double eps, ODESOLUTION *sol)
{
    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;

    gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, sol->dim);
    gsl_odeiv2_control * c  = gsl_odeiv2_control_y_new (eps, 0.0);
    gsl_odeiv2_evolve * e   = gsl_odeiv2_evolve_alloc (sol->dim);

    gsl_odeiv2_system sys = {f, NULL, (size_t)sol->dim, p};

    double t = w1, t1 = w2;
    double h = eps;
    double *y;

    y = new double[sol->dim];
    for(int i=0;i<sol->dim;i++){
         y[i]=y0[i];
         sol->ytab[i][0]=y0[i];
    }
    sol->wtab[0]=w1;

    int status,cnt=1;

    while (t < t1)
    {
      status = gsl_odeiv2_evolve_apply (e, c, s,&sys,&t, t1,&h, y);

      if (status != GSL_SUCCESS)
          break;

      sol->wtab[cnt]=t;
      for(int i=0;i<sol->dim;i++){
          sol->ytab[i][cnt]=y[i];
      }
      cnt++;
      if(cnt==sol->n) break;
    }

    sol->nactual=cnt;
    
    delete y;

    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);
  
    return 0;
}

int findroot(FUNC(f),void *p, double w1, double w2, double eps, double *r)
{
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double rv;
    double x_lo = w1, x_hi = w2;
    gsl_function F;
  
    F.function = f;
    F.params = p;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);

    do{
      iter++;
      status = gsl_root_fsolver_iterate (s);
      rv = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi, 0, eps);
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    *r = rv;

    return status;
}