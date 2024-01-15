#include "cycles.h"

void Initialize(cycles_struct *cycles, N_Vector CV_Y, void **cvode_mem)
{
    int             kx;
    control_struct *control = &cycles->control;

    control->stepsize = STEPSIZE;

    for (kx = 0; kx < number_of_columns; kx++)
    {
        int             kz;
        soil_struct    *soil;
        wstate_struct  *ws;
        wflux_struct   *wf;

        soil = &cycles->grid[kx].soil;
        ws = &cycles->grid[kx].ws;
        wf = &cycles->grid[kx].wf;

        for (kz = 0; kz < number_of_layers; kz++)
        {
            soil->dsat[kz] = soil->b[kz] * soil->ksatv[kz] * soil->air_entry_potential[kz] / soil->porosity[kz];
            ws->potential[kz] = SoilWaterPotential(soil->porosity[kz], soil->air_entry_potential[kz], soil->b[kz], ws->smc[kz]);

            // Initialize CVode variables
            NV_Ith_S(CV_Y, INDEX(kx, kz)) = ws->potential[kz];

            wf->lateral[kz] = 0.0;
        }
    }

    // Allocate memory for solver
    *cvode_mem = CVodeCreate(CV_BDF);
    if (*cvode_mem == NULL)
    {
        printf("Error in allocating memory for solver.\n");
        exit(EXIT_FAILURE);
    }

    return;
}
