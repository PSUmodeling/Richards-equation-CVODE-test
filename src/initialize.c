#include "cycles.h"

void Initialize(cycles_struct *cycles, N_Vector CV_Y, void **cvode_mem)
{
    int             kz;
    soil_struct    *soil = &cycles->soil;
    phystate_struct *phys = &cycles->phys;
    ctrl_struct    *ctrl = &cycles->ctrl;
    wstate_struct  *ws = &cycles->ws;

    for (kz = 0; kz < NSOIL; kz++)
    {
        soil->porosity[kz] = 0.464;
        soil->air_entry_pot[kz] = 0.62;
        soil->b[kz] = 8.72;
        soil->ksat[kz] = 0.20E-5;
        soil->dsat[kz] = soil->b[kz] * soil->ksat[kz] * soil->air_entry_pot[kz] / soil->porosity[kz];
    }
    soil->slope = 0.1;

    ctrl->stepsize = STEPSIZE;

    for (kz = 0; kz < NSOIL; kz++)
    {
        phys->zsoil[kz] = -0.1 * (kz + 1);
        phys->soil_depth[kz] = 0.1;
        ws->smc[kz] = (kz < 2) ? 0.3568 : 0.464;
        ws->potential[kz] = SoilWaterPot(soil->porosity[kz], soil->air_entry_pot[kz], soil->b[kz], ws->smc[kz]);
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
