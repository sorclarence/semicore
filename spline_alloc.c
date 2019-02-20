#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "allvars.h"
#include "proto.h"

void spline_alloc()
{
  acc_mass = gsl_interp_accel_alloc();
  acc_mass_over_rsqu = gsl_interp_accel_alloc();
  acc_rmin = gsl_interp_accel_alloc();
  acc_rmax = gsl_interp_accel_alloc();

  acc_mass_afk = gsl_interp_accel_alloc();

  sp_mass = gsl_spline_alloc(gsl_interp_cspline, N_SP); 
  sp_mass_over_rsqu = gsl_spline_alloc(gsl_interp_cspline, N_SP);
  sp_rmin = gsl_spline_alloc(gsl_interp_cspline, N_SP); 
  sp_rmax = gsl_spline_alloc(gsl_interp_cspline, N_SP); 

  sp_mass_afk = gsl_spline_alloc(gsl_interp_cspline, N_SP);
}

void spline_free()
{
  gsl_spline_free(sp_mass);
  gsl_spline_free(sp_mass_over_rsqu);
  gsl_spline_free(sp_rmin);
  gsl_spline_free(sp_rmax);
  gsl_spline_free(sp_mass_afk);

  gsl_interp_accel_free(acc_mass);
  gsl_interp_accel_free(acc_mass_over_rsqu);
  gsl_interp_accel_free(acc_rmin);
  gsl_interp_accel_free(acc_rmax);
  gsl_interp_accel_free(acc_mass_afk);
}
