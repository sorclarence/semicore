#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include "define.h"

/* Timer */
extern time_t T_start;
extern time_t T_end;

/* Time Evolution Related*/
extern double TimeSpan;            // time span for system evolution 
extern double DeltaFrac;           // Decay faction unit 

/* cosmology */
extern double Rho_crit_0;          // Nowadays critical density S
extern double Delta_c;             // spherical collapse halo overdensity 

/* input nfw model */
extern double Rs;                  // characteristic radius
extern double Rhos;                // characteristic density
extern double Rvir;                // halo virial mass


/* spline related */
extern double Rmin;                        // minimum interpolation radius
extern double Rmax;                        // maximum interpolation radius

extern double K_Mass_in;                     // cumulative mass profile inner slope
extern double K_Mass_out;                    // cumulative mass profile outer slope
extern double A_Mass_in;                     // cumulative mass profile inner power law normalization
extern double A_Mass_out;                    // cumulative mass profile outer power law normalization

/* after kick related */
//------------------------------------------------------------------------------------------------------------
extern double SpMass_afk[N_SP];                // cumulative mass

extern gsl_interp_accel *acc_mass_afk;
extern gsl_spline *sp_mass_afk;                // spline for after-kick cumulative mass
//-------------------------------------------------------------------------------------------------------------

extern double K_rmin_in;                   // rmin interpolation inner slope
extern double K_rmin_out;                  // rmin interpolation outer slope
extern double K_rmax_in;                   // rmax interpolation inner slope
extern double K_rmax_out;                   // rmax interpolation outer slope
extern double A_rmin_in;                   // rmin interpolation inner power normalization
extern double A_rmin_out;                  // rmin interpolation outer power normalization
extern double A_rmax_in;                   // rmax interpolation inner power normalization
extern double A_rmax_out;                   // rmax interpolation outer power normalization

extern double SpR[N_SP];                   // radius points
extern double SpRhoMom[N_SP];              // density points for mother particle
extern double SpRhoDau[N_SP];              // density points for daughter particle
extern double SpMass[N_SP];                // total cumulative mass
extern double SpMassMom[N_SP];             // mother cumulative mass
extern double SpMassDau[N_SP];             // daughter cumulative mass
extern double SpMassOverRsqu[N_SP];         // total cumulative mass over r^2 
extern double SpRmin[N_SP];                // rmin for an ellipse orbit
extern double SpRmax[N_SP];                // rmax for an ellipse orbit

extern gsl_interp_accel *acc_mass;
extern gsl_interp_accel *acc_mass_over_rsqu;
extern gsl_interp_accel *acc_rmin;
extern gsl_interp_accel *acc_rmax;

extern gsl_spline *sp_mass;                // spline for cumulative mass
extern gsl_spline *sp_mass_over_rsqu;  
extern gsl_spline *sp_rmin;                // spline for rmin
extern gsl_spline *sp_rmax;                // spline for rmax

/* variables for power law function */
extern double PowerSlop;                   // power index 
extern double PowerNorm;                   // power law nomalization 
