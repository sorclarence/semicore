#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "allvars.h"
#include "proto.h"

/*
*  Unit system of this code lines with
   Gadget's unit system.
*/

void run()
{
   int k, n;

  for(n = 0; n < N_TP; n++)
    {
       printf("Step: %d\n", n);

       spline_alloc();
      // set up spline for cumulative mass
       printf("---> set up cubic spline for cumulative mass...\n");
       mass_spline_setup();
       printf("-----> done.\n");

      // set up spline for rmin and rmax
       printf("---> set up cubic spline for rmin and rmax...\n");
       rm_spline_setup();
       printf("-----> done!\n");

      //update mass profile 
       printf("---> updating mass profile...\n");
       for(k = 0; k < N_SP; k++)
         {
            SpMass_afk[k] = mass_after_kick(SpR[k]);
            SpMassMom[k] -= DeltaFrac*mass_nfw(SpR[k]);
            SpMassDau[k] += DeltaFrac*SpMass_afk[k];
            SpMass[k] = SpMassMom[k] + SpMassDau[k];
            SpMassOverRsqu[k] = SpMass[k]/SpR[k]/SpR[k];
         }
       printf("-----> done!\n");
       spline_free();
       printf("   Step %d done.\n", n);
     }
}

/*
 * this function sets up
 * the spline for cumulative mass.
*/
void mass_spline_setup()
{
   // set up spline for cumulative mass
   gsl_spline_init(sp_mass, SpR, SpMass, N_SP);
   power_law_setup(SpR[0], SpR[1], SpMass[0], SpMass[1]);
   K_Mass_in = PowerSlop;
   A_Mass_in = PowerNorm;

   power_law_setup(SpR[N_SP-2], SpR[N_SP-1], SpMass[N_SP-2], SpMass[N_SP-1]);
   K_Mass_out = PowerSlop;
   A_Mass_out = PowerNorm;

   // set up for potential integration 
   gsl_spline_init(sp_mass_over_rsqu, SpR, SpMassOverRsqu, N_SP);
}

/*
 * this function returns
 * cumulative mass at radius r.
 * main part is a cubic spline.
 * for inner halo part and outer halo
 * part, power law interpolation is
 * deployed.
*/
double Mass(double r)
{
/* inner power law interpolation */
  if(r <= SpR[0])
     return A_Mass_in*pow(r, K_Mass_in);
/* outer power law interpolation */
  else if(r >= SpR[N_SP-1])
     return A_Mass_out*pow(r, K_Mass_out);
/* intermediate part */
  else
    return gsl_spline_eval(sp_mass, r, acc_mass);
}

/*
* this function returns
* gravitational potential
* of M(r).
*/

double Phi(double r)
{
   double phi, res, err;
   double rmin = SpR[0];
   double rmax = SpR[N_SP-1];
   double param;
   size_t neval;

/*
   Integration check
*/
    if(K_Mass_out - 1.0 >= 0)
     {
       printf("Potential integration diverges. We better stops.\n");
       exit(0);
     }

    if(r >= rmax)
       phi = G*A_Mass_out/(K_Mass_out-1.0)*pow(r, K_Mass_out-1.0);
    else if(r >= rmin)
     {
       gsl_integration_cquad_workspace *w
         = gsl_integration_cquad_workspace_alloc(100);

       gsl_function F;
       F.function = &integ_kernel_phi;
       F.params = &param;

       gsl_integration_cquad(&F, r, rmax, 0, 1E-7, w, &res, &err, &neval);
       gsl_integration_cquad_workspace_free(w);

       phi = Phi(rmax)-G*res;
     }
    else
       phi = Phi(rmin)-G*A_Mass_in/(K_Mass_in-1.0)*(pow(rmin, K_Mass_in-1.0)-pow(r, K_Mass_in-1.0));

    return phi;
}

/*
 * this function returns 
 * effective potential 
*/
double Veff(double r, double r0)
{
    double j2 = G*r0*Mass(r0);                  // specific angular momentum squared at radius r0
    double res = j2/(2.0 * r * r) + Phi(r);

    return res;
}

/*
 * this function sets up
 * spline for rmin and rmax
 * NOTE: This function uses
 * 1.0E10 as a proxy of infinity 
 */

void rm_spline_setup()
{
   int k;

   int status;
   int iter = 0, max_iter = 100;
   double r0;
   double x_lo, x_hi;
   double rmin, rmax;
   double norm, power;   // power-law normalization and power
   const gsl_root_fsolver_type *T;
   T = gsl_root_fsolver_brent;
   gsl_root_fsolver *s_min;
   gsl_root_fsolver *s_max;
   s_min = gsl_root_fsolver_alloc(T);
   s_max = gsl_root_fsolver_alloc(T);
 
   double gap;

   gsl_function F;
   F.function = &roots_kernel;

   for(k = 0; k < N_SP; k++)
    {
       r0 = SpR[k];
       F.params = &r0;
     
       iter = 0; 
       gsl_root_fsolver_set(s_min, &F, 1E-10, r0);
       do
        {
         iter++;
         status = gsl_root_fsolver_iterate(s_min);
         rmin = gsl_root_fsolver_root(s_min);
         x_lo = gsl_root_fsolver_x_lower(s_min);
         x_hi = gsl_root_fsolver_x_upper(s_min);
         status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);
        } 
       while(status == GSL_CONTINUE && iter < max_iter);
       SpRmin[k] = rmin;
      
       gap = Veff(10*Rmax, r0) - Veff(r0, r0) - 0.5*VK*VK; 
       if(gap > 0.0)                                       // decide whether rmax as a root exists
        {   
          iter = 0;
          gsl_root_fsolver_set(s_max, &F, r0, 10*Rmax);
          do
           {
            iter++;
            status = gsl_root_fsolver_iterate(s_max);
            rmax = gsl_root_fsolver_root(s_max);
            x_lo = gsl_root_fsolver_x_lower(s_max);
            x_hi = gsl_root_fsolver_x_upper(s_max);
            status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);
           } 
          while(status == GSL_CONTINUE && iter < max_iter);
          SpRmax[k] = rmax;
        }
       else                                   // using last two points for extrapolation 
        {
          if(k < 2)
           {  printf("\n---------------------------------\n");
              printf("TOO LARGE VK TO FORM HALO\n");
              printf("---------------------------------\n");
              exit(0);
           }
          else
           {
              power_law_setup(SpR[k-2], SpR[k-1], SpRmax[k-2], SpRmax[k-1]);
              norm = PowerNorm;
              power = PowerSlop;
              SpRmax[k] = norm*pow(SpR[k], power);
           }
        }
    }
   gsl_root_fsolver_free(s_min);
   gsl_root_fsolver_free(s_max);


   power_law_setup(SpR[0], SpR[1], SpRmin[0], SpRmin[1]);
   K_rmin_in = PowerSlop;
   A_rmin_in = PowerNorm;

   power_law_setup(SpR[0], SpR[1], SpRmax[0], SpRmax[1]);
   K_rmax_in = PowerSlop;
   A_rmax_in = PowerNorm;

   power_law_setup(SpR[N_SP-2], SpR[N_SP-1], SpRmin[N_SP-2], SpRmin[N_SP-1]);
   K_rmin_out = PowerSlop;
   A_rmin_out = PowerNorm;

   power_law_setup(SpR[N_SP-2], SpR[N_SP-1], SpRmax[N_SP-2], SpRmax[N_SP-1]);
   K_rmax_out = PowerSlop;
   A_rmax_out = PowerNorm;

/* set up spline for Rmin and Rmax */
   gsl_spline_init(sp_rmin, SpR, SpRmin, N_SP);
   gsl_spline_init(sp_rmax, SpR, SpRmax, N_SP);
}

/*
* this function returns 
* rmin for a given r0. 
* power-law is used to
* to interpolate the results
* at inner halo region
* and outer halo region.
*/
double rmin(double r0)
{
  double res;
/* inner power law interpolation */
  if(r0 <= SpR[0])
     res = A_rmin_in*pow(r0, K_rmin_in);
/* outer power law interpolation */
  else if(r0 >= SpR[N_SP-1])
     res = A_rmin_out*pow(r0, K_rmin_out);
/* intermediate part */
  else
     res = gsl_spline_eval(sp_rmin, r0, acc_rmin);
  
  return res;
}

/*
* this function returns 
* rmax for a given r0. 
* power-law is used to
* to interpolate the results
* at inner halo region
* and outer halo region.
*/
double rmax(double r0)
{
  double res;

/* inner power law interpolation */
  if(r0 <= SpR[0])
     res = A_rmax_in*pow(r0, K_rmax_in);
/* outer power law interpolation */
  else if(r0 >= SpR[N_SP-1])
     res = A_rmax_out*pow(r0, K_rmax_out);
/* intermediate part */
  else
     res = gsl_spline_eval(sp_rmax, r0, acc_rmax);

  return res;
}

/* the function below calculates 
 * the g-function
*/
double gfunc(double r, double r0)
{
   double g_res;
   double rm, rx;
   double x;             // ratio between rmax and rmin
   double a;             // semi-major axis
   double e;             // eccentricity

   double theta, cos_theta;

   double result, error;
   double res, err;
   size_t neval;

   rm = rmin(r0);
   rx = rmax(r0);
   x = rx/rm;
   a = (rx + rm) / 2.0;
   e = (x - 1.0) / (x + 1.0);
 
   if(r < rm)
      g_res = 0;
   else if (r > rx)
      g_res = 1.0;
   else
    {
      if((1.0-e)<1E-4)        // for eccentricity closing to 1.0
        g_res = 0;
      else
       {
          cos_theta = (a * (1.0 - e*e) - r) / (e*r);
          theta = acos(cos_theta);
     
          gsl_integration_cquad_workspace *w
            = gsl_integration_cquad_workspace_alloc(100);
          gsl_function F;
          F.function = &integ_ellipse;
          F.params = &e;
         
          gsl_integration_cquad(&F, 0, theta, 0, 1e-7, w, &result, &error, &neval);
/*
 * This integration will diverge if eccentricity
 * e is very close to 1.0
*/
          gsl_integration_cquad(&F, 0, PI, 0, 1e-7, w, &res, &err, &neval);   
 
          gsl_integration_cquad_workspace_free(w);
          g_res =  result/res;
       }
    }
  return g_res;
}


/*
 * integration kernel function
 * for Phi(r)
*/
double integ_kernel_phi(double r, void *params)
{
/* inner power law interpolation */
  if(r <= SpR[0])
     return A_Mass_in*pow(r, K_Mass_in-2.0);
/* outer power law interpolation */
  else if(r >= SpR[N_SP-1])
     return A_Mass_out*pow(r, K_Mass_out-2.0);
/* intermediate part */
  else
    return gsl_spline_eval(sp_mass_over_rsqu, r, acc_mass_over_rsqu);
}

/*
 * the two roots of this function
 * are rmin and rmax.
*/
double roots_kernel(double r, void *params)
{
  double r0 = *(double *) params;
  double veff, veff_min;
 
  veff     = Veff(r, r0);
  veff_min = Veff(r0, r0);

  return veff - veff_min - 0.5*VK*VK;
}

/*
 * Integral kernel for g-function
*/
double integ_ellipse(double phi, void *params)
{
  double ecc = *(double *) params;        // eccentricity of the ellipse 
  double result = (1.0 + ecc*cos(phi)) * (1.0 + ecc*cos(phi));

  return 1.0 / result;
}
