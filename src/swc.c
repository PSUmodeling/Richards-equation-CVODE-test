#include "cycles.h"

void SWC(int cstep, cycles_struct *cycles, void *cvode_mem, N_Vector CV_Y)
{
    double              t;
    double              total_depth = 0.0;
    double              wplus = 0.0;
    double              runoff3 = 0.0;
    int                 kz;
    ctrl_struct        *ctrl = &cycles->ctrl;
    wflux_struct       *wf = &cycles->wf;
    wstate_struct      *ws = &cycles->ws;
    phystate_struct    *phys = &cycles->phys;
    forc_struct        *forc = &cycles->forc;
    soil_struct        *soil = &cycles->soil;

    t = (cstep + 1) * ctrl->stepsize;

    wf->infil = forc->infil[cstep];
    wf->soil_evap = forc->soil_evap[cstep];
    for (kz = 0; kz < NSOIL; kz++)
    {
        total_depth += phys->soil_depth[kz];
    }
    for (kz = 0; kz < NSOIL; kz++)
    {
        wf->uptake[kz] = forc->uptake[cstep] * phys->soil_depth[kz] / total_depth;
    }

    // Solve Richards equation ODE using CVode
    SolveCVode((realtype)t, cvode_mem, CV_Y);

    // Get the solution, and constrain it to be between 0.02 and porosity
    for (kz = 0; kz < NSOIL; kz++)
    {
        ws->smc[kz] = SoilWaterContent(soil->porosity[kz], soil->air_entry_pot[kz], soil->b[kz], NV_Ith_S(CV_Y, kz)) +
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
