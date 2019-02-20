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
void finalize()
{
    int k;
    char output_fname[300];
    FILE *fp;
   
    /* output result */
    sprintf(output_fname, "/home/x3/Projects/ddm/adiabaticExpan/data/l6n256/ddm_denspro/convergence_test/ddm_masspro_vk%.1lf_NSP%d", VK, N_SP);
    fp = fopen(output_fname,"w");
    fprintf(fp, "      radius          mass\n");
    fprintf(fp, "==========================================================================\n");
    for(k = 0; k < N_SP; k++)
     {

       fprintf(fp, "     %.4E        %.6E \n", SpR[k], SpMass[k]);
       fflush(fp);
     }
    fclose(fp);
}
