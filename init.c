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

void init()
{
  int k;
  double dr;

/* click the timer */
  T_start = time(NULL);

  printf("Initializing the code...\n");
/* time evolution related */
  TimeSpan = age_of_the_system();
  DeltaFrac = (1.0-exp(-log(2.0)*TimeSpan/TAU))/(double) N_TP;

/* initialize NFW profile */
  Rho_crit_0 = 3.0 * H0 * H0 / (8.0 * PI * G);

  // Bryan & Norman (1998) Delta_c 
  Delta_c = 18.0*PI*PI - 82.0*OMEGA_LAMBDA_0 
                       - 39.0*OMEGA_LAMBDA_0*OMEGA_LAMBDA_0;

  // Navarro et al. (1996)
  Rhos = Rho_crit_0 * Delta_c * pow(C_NFW, 3)/3.0/(log(1.0 + C_NFW) - C_NFW/(1 + C_NFW));

  // According to definition 
  Rs   = pow(3.0*M_NFW/(4.0*PI*Delta_c*Rho_crit_0), 1.0/3.0)/C_NFW;

  // Halo Virial Mass 
  Rvir = C_NFW * Rs;
 
/* set up spline beginning and ending nodes */ 
  Rmin = 1E-3*Rvir;
  Rmax = Rvir;
  dr = (Rmax - Rmin) / (double) N_SP;

/* set up nodes positions */
  for(k = 0; k < N_SP; k++)
   {
      SpR[k] = Rmin + k*dr;

/* initialize input cumulative mass profile */
      SpMassMom[k] = mass_nfw(SpR[k]);
      SpMassDau[k] = 0;
      SpMass[k] = SpMassMom[k] + SpMassDau[k];
      SpMassOverRsqu[k] = SpMass[k]/SpR[k]/SpR[k];
   }
   printf("Initialization finished.\n");
}

//-------------------------------------------------
// density profile related analytical functions 
//-------------------------------------------------
/*
 * here writes the nfw profile
*/
double rho_nfw(double r)
{
   double x = r/Rs;
   return Rhos/x/(1.0+x)/(1.0+x);
}

double rho_nfw_deriv(double r)
{
   double x = r/Rs;
   double term1 = -Rhos/Rs*pow(x, -2)*pow(1.0+x, -2);
   double term2 = -2.0*Rhos/Rs*pow(x, -1)*pow(1.0+x, -3);
 
   return term1 + term2;
}

/*
  here writes Burkert profile
*/
double rho_bur(double r)
{
   double x = r/Rs;
   return Rhos/(1.0+x)/(1.0+x*x);
}

/*
 * here writes the cumulative mass profile
 * of nfw density profile
*/
double mass_nfw(double r)
{
   double x = r/Rs;
   return 4.0*PI*Rhos*pow(Rs, 3)*(log(1.0+x) - x/(1.0+x)); 
}

/* this function calculates the time span of halo evolution 
 *  Unit is Gyr.
*/
double age_of_the_system()
{
  gsl_function F;
  gsl_integration_workspace *workspace;
  double result0, result1, abserr;
  double a0 = 1e-8;
  double b0 = TIMEBEGIN;
  double b1 = TIMEEND;

  workspace = gsl_integration_workspace_alloc(1000);
  F.function = &interg_a;

  gsl_integration_qag(&F,a0, b0, 0, 1.0e-8, 1000, GSL_INTEG_GAUSS61, workspace, &result0, &abserr);
  gsl_integration_qag(&F,a0, b1, 0, 1.0e-8, 1000, GSL_INTEG_GAUSS61, workspace, &result1, &abserr);

  gsl_integration_workspace_free(workspace);

  return (result1 - result0) / H0 * 0.978 / HUBBLE;   // convert to Gyr
}

double interg_a(double a, void *param)
{
  return 1.0 / sqrt((1.0-OMEGA_LAMBDA_0) / a + OMEGA_LAMBDA_0*a*a);
}

/* this function calculates
 * the power index k and 
 * the power normalization a
 * for a power-law function:
 * y = ax^k
*/

void power_law_setup(double x1, double x2, double y1, double y2)
{
// clean up for calculation 
    PowerSlop = 0;
    PowerNorm = 0;

    PowerSlop = (log10(y1) - log10(y2)) / (log10(x1) - log10(x2));
    PowerNorm = pow(10, log10(y1) - PowerSlop*log10(x1));
}
