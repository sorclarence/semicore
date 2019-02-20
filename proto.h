#ifdef ALLVAR_H
   #include "allvars.h"
#endif

void init();
void power_law_setup(double x1, double x2, double y1, double y2);
void run();
double rho_nfw(double r);
double rho_nfw_deriv(double r);
double rho_bur(double r);
double Mass(double r);
double mass_nfw(double r);
double Phi(double r);
double Veff(double r, double r0);
double integ_kernel_phi(double r, void *params);
double roots_kernel(double r0, void *params);
double rmin_root_kernel(double x, void *params);
double rmax_root_kernel(double x, void *params);
void   mass_spline_setup();
void   rm_spline_setup();
double rmin(double r0);
double rmax(double r0);
double integ_ellipse(double phi, void *params);
double integ_mass(double r0, void *params);
double gfunc(double r, double r0);
double mass_after_kick(double r);
void finalize();
double age_of_the_system();
double interg_a(double a, void *param);
void   spline_alloc();
void   spline_free();
void acknowledge();
