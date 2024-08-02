#include "cycles.h"

#define SAND_POROSITY 0.339
#define SAND_PSIS 0.07
#define SAND_KS 5.8E-5
#define SAND_B 1.7

#define CLAY_POROSITY 0.468
#define CLAY_PSIS 0.37
#define CLAY_KS 1.7E-7
#define CLAY_B 7.6

void Initialize(cycles_struct *cycles, N_Vector CV_Y, void **cvode_mem)
{
    int             kz;
    soil_struct    *soil = &cycles->soil;
    phystate_struct *phys = &cycles->phys;
    ctrl_struct    *ctrl = &cycles->ctrl;
    wstate_struct  *ws = &cycles->ws;

    for (kz = 0; kz < NSOIL; kz++)
    {
        soil->porosity[kz] = SAND_POROSITY;
        soil->air_entry_pot[kz] = SAND_PSIS;
        soil->b[kz] = SAND_B;
        soil->ksat[kz] = SAND_KS;
        soil->dsat[kz] = soil->b[kz] * soil->ksat[kz] * soil->air_entry_pot[kz] / soil->porosity[kz];
        //soil->porosity[kz] = 0.464;
        //soil->air_entry_pot[kz] = 0.62;
        //soil->b[kz] = 8.72;
        //soil->ksat[kz] = 0.20E-5;
        //soil->dsat[kz] = soil->b[kz] * soil->ksat[kz] * soil->air_entry_pot[kz] / soil->porosity[kz];
    }

    soil->porosity[0] = CLAY_POROSITY;
    soil->air_entry_pot[0] = CLAY_PSIS;
    soil->b[0] = CLAY_B;
    soil->ksat[0] = CLAY_KS;

    soil->slope = 0.1;

    ctrl->stepsize = STEPSIZE;

    phys->zsoil[0] = -0.1;
    phys->zsoil[1] = -0.24;
    phys->zsoil[2] = -0.84;
    phys->zsoil[3] = -1.84;

    phys->soil_depth[0] = 0.1;
    phys->soil_depth[1] = 0.14;
    phys->soil_depth[2] = 0.6;
    phys->soil_depth[3] = 1.0;

    ws->smc[0] = 0.3;
    ws->smc[1] = 0.3;
    ws->smc[2] = 0.2;
    ws->smc[3] = 0.25;
    //ws->smc[0] = 0.3568;
    //ws->smc[1] = 0.3601;
    //ws->smc[2] = 0.3669;
    //ws->smc[3] = 0.3793;

    phys->wetting_front_depth = 0.0;

    for (kz = 0; kz < NSOIL; kz++)
    {
        ws->potential[kz] = SoilWaterPot(soil->porosity[kz], soil->air_entry_pot[kz], soil->b[kz], ws->smc[kz]);

        phys->smc_below[kz] = ws->smc[kz];
        phys->smc_wf[kz] = 0.0;
    }

    // Allocate memory for solver
    *cvode_mem = CVodeCreate(CV_BDF);
    if (*cvode_mem == NULL)
    {
        printf("Error in allocating memory for solver.\n");
        exit(EXIT_FAILURE);
    }

    // Initialize CVode variables
    for (kz = 0; kz < NSOIL; kz++)
    {
        NV_Ith_S(CV_Y, kz) = ws->potential[kz];
    }

    return;
}
