#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "allvars.h"
#include "proto.h"

/* this function calculates the 
 * main integral:
 * M_out(r) = int_0^{inf}4*pi*x*x*rho_nfw(x)*g(r, x)dx
*/

double mass_after_kick(double r)
{
   double x_lower, x_upper;       // interval where gfunc deviates from a step-function
   int status;
   int iter = 0, max_iter = 100;
   double x_lo, x_hi;
   const gsl_root_fsolver_type *T;
   T = gsl_root_fsolver_brent;
   gsl_root_fsolver *x_min;
   gsl_root_fsolver *x_max;
   x_min = gsl_root_fsolver_alloc(T);
   x_max = gsl_root_fsolver_alloc(T);

   gsl_function RT;

// evaluate integration lower bound
   if(r <= rmax(1E-10))                         // no root can be found
     x_lower = 1E-10; 
   else                                          // call gsl function 
    {
      RT.function = &rmax_root_kernel;
      RT.params   = &r;
      gsl_root_fsolver_set(x_min, &RT, 1E-10, r);
      iter = 0; 
      do
       {
        iter++;
        status = gsl_root_fsolver_iterate(x_min);
        x_lower = gsl_root_fsolver_root(x_min);
        x_lo = gsl_root_fsolver_x_lower(x_min);
        x_hi = gsl_root_fsolver_x_upper(x_min);
        status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);
       } 
       while(status == GSL_CONTINUE && iter < max_iter);
     }

// evaluate integration upper bound
    
   RT.function = &rmin_root_kernel;
   RT.params   = &r;
   gsl_root_fsolver_set(x_max, &RT, r, 10*Rmax);

   iter = 0; 
   do
    {
     iter++;
     status = gsl_root_fsolver_iterate(x_max);
     x_upper = gsl_root_fsolver_root(x_max);
     x_lo = gsl_root_fsolver_x_lower(x_max);
     x_hi = gsl_root_fsolver_x_upper(x_max);
     status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);
    } 
    while(status == GSL_CONTINUE && iter < max_iter);
   gsl_root_fsolver_free(x_min);
   gsl_root_fsolver_free(x_max);

// integration

   gsl_integration_cquad_workspace *w
     = gsl_integration_cquad_workspace_alloc(100);
   double res, err;
   size_t neval;

   gsl_function F;
   F.function = &integ_mass;
   F.params = &r;
   
   gsl_integration_cquad(&F, x_lower, x_upper, 0, 1e-7, w, &res, &err, &neval);
   gsl_integration_cquad_workspace_free(w);

   return res + mass_nfw(x_lower);
}

/*
 * integral kernel for mass_out
 * note that the density profile
 * in this kernel is the input
 * density profile at first
*/
double integ_mass(double x, void *params)
{
   double r = *(double *) params;

   return 4.0*PI*x*x*rho_nfw(x)*gfunc(r, x);
}
/*
 * kernel function for solving 
 * the function
 * rmin(x) - r = 0
*/
double rmin_root_kernel(double x, void *params)
{
   double r = *(double *) params;
   return rmin(x) - r;
}

/*
 * kernel function for solving 
 * the function
 * rmax(x) - r = 0
*/
double rmax_root_kernel(double x, void *params)
{
   double r = *(double *) params;
   return rmax(x) - r;
}
