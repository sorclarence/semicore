/* CONSTANTS */
#define PI                          3.1415926535897932384626433832795028842              /* pi */
#define G                           43007.1     //kpc*(km/s)^2/10^10Msun
#define H0                          0.1        // h km/s/kpc

/* DDM PARAMETERS */
#define TIMEBEGIN                   0.01           // starting time scale factor
#define TIMEEND                     1.0            // ending time scale factor 
#define VK                          20.0           // km/s
#define TAU                         1.0            // Gyr

/* Flat cosmology */
#define OMEGA_LAMBDA_0              0.6834           // Nowadays dark energy density parameter
#define HUBBLE                      0.6727         // dimensionless Hubble constant 

/*
 * Input NFW model 
*/

#define C_NFW                       23.6            // nfw concentration
#define M_NFW                       0.505          // nfw virial mass 10^10 / h solar mass

/*
 * SPLINE RELATED
*/
#define N_SP                         1000           // number of data points for spline construction
#define N_TP                           10          // number of timestep
