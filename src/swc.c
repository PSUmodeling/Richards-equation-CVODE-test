#include "cycles.h"

void SWC(int cstep, cycles_struct *cycles, void *cvode_mem, N_Vector CV_Y)
{
    double              t;
    int                 kx;
    forcing_struct     *forcing = &cycles->forcing;

    t = (cstep + 1) * cycles->control.stepsize;

    for (kx = 0; kx < number_of_columns; kx++)
    {
        int                 kz;
        double              total_depth = 0.0;
        wflux_struct       *wf = &cycles->grid[kx].wf;
        phystate_struct    *phys = &cycles->grid[kx].phys;

        wf->infiltration = forcing->infiltration[cstep];
        wf->soil_evaporation = forcing->soil_evaporation[cstep];

        for (kz = 0; kz < number_of_columns; kz++)
        {
            total_depth += phys->soil_depth[kz];
        }

        for (kz = 0; kz < number_of_columns; kz++)
        {
            wf->uptake[kz] = forcing->uptake[cstep] * phys->soil_depth[kz] / total_depth;
        }
    }

    // Solve Richards equation ODE using CVode
    SolveCVode((realtype)t, cvode_mem, CV_Y);

    // Get the solution, and constrain it to be between 0.02 and porosity
    for (kx = 0; kx < number_of_columns; kx++)
    {
        int                 kz;
        double              total_depth = 0.0;
        double              wplus = 0.0;
        double              runoff3 = 0.0;
        wstate_struct      *ws = &cycles->grid[kx].ws;
        soil_struct        *soil = &cycles->grid[kx].soil;
        phystate_struct    *phys = &cycles->grid[kx].phys;

        for (kz = 0; kz < number_of_layers; kz++)
        {
            ws->smc[kz] =
                SoilWaterContent(soil->porosity[kz], soil->air_entry_potential[kz], soil->b[kz], NV_Ith_S(CV_Y, kz)) +
                wplus / phys->soil_depth[kz];

            if (ws->smc[kz] > 0.02)
            {
                if (ws->smc[kz] > soil->porosity[kz])
                {
                    wplus = (ws->smc[kz] - soil->porosity[kz]) * phys->soil_depth[kz];
                    ws->smc[kz] = soil->porosity[kz];
                }
                else
                {
                    wplus = 0.0;
                }
            }
            ws->smc[kz] = MAX(ws->smc[kz], 0.02);
        }

        runoff3 = wplus;
    }

}
