#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "define.h"
#include "allvars.h"

/* Timer */
time_t T_start;
time_t T_end;

/* Evolution time */
double TimeSpan;            // time span for system evolution 
double DeltaFrac;           // Decay fraction unit 

/* cosmology */
double Rho_crit_0;          // Nowadays critical density 
double Delta_c;             // spherical collapse overdensity

/* input nfw model */
double Rs;                  // characteristic radius
double Rhos;                // characteristic density
double Rvir;                // halo virial mass 

/* spline related */
double Rmin;                        // minimum interpolation radius
double Rmax;                        // maximum interpolation radius
double K_Mass_in;                   // cumulative mass profile inner slope
double K_Mass_out;                  // cumulative mass profile outer slope
double A_Mass_in;                   // cumulative mass profile inner power law normalization
double A_Mass_out;                  // cumulative mass profile outer power law normalization

/* after kick related */
//------------------------------------------------------------------------------------------------------------
double SpMass_afk[N_SP];                // after kick cumulative mass

gsl_interp_accel *acc_mass_afk;
gsl_spline       *sp_mass_afk;          // spline for after-kick cumulative mass
//-------------------------------------------------------------------------------------------------------------

double K_rmin_in;                   // rmin interpolation inner slope
double K_rmin_out;                  // rmin interpolation outer slope
double K_rmax_in;                   // rmax interpolation inner slope
double K_rmax_out;                   // rmax interpolation outer slope
double A_rmin_in;                   // rmin interpolation inner power normalization
double A_rmin_out;                  // rmin interpolation outer power normalization
double A_rmax_in;                   // rmax interpolation inner power normalization
double A_rmax_out;                   // rmax interpolation outer power normalization

double SpR[N_SP];                   // radius points
double SpRhoMom[N_SP];              // density points for mother particle
double SpRhoDau[N_SP];              // density points for daughter particle
double SpMass[N_SP];                // total cumulative mass
double SpMassMom[N_SP];             // mother cumulative mass
double SpMassDau[N_SP];             // daughter cumulative mass
double SpMassOverRsqu[N_SP];        // total cumulative mass over r square
double SpRmin[N_SP];                // rmin for an ellipse orbit
double SpRmax[N_SP];                // rmax for an ellipse orbit

gsl_interp_accel *acc_mass;
gsl_interp_accel *acc_mass_over_rsqu;
gsl_interp_accel *acc_rmin;
gsl_interp_accel *acc_rmax;

gsl_spline *sp_mass;               // spline for cumulative mass
gsl_spline *sp_mass_over_rsqu;
gsl_spline *sp_rmin;                // spline for rmin
gsl_spline *sp_rmax;                // spline for rmax

/* variables for power law function */
double PowerSlop;                   // power index
double PowerNorm;                   // power normalization 
