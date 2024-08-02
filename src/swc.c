#include "cycles.h"

void SWC(int cstep, cycles_struct *cycles, void *cvode_mem, N_Vector CV_Y)
{
    double              t;
    double              total_depth = 0.0;
    double              wplus = 0.0;
    double              runoff3 = 0.0;
    double              surface_water = 0.0;
    int                 kz;
    int                 n;
    double              wetting_front_thickness[NSOIL];
    double              infiltration_capacity;
    ctrl_struct        *ctrl = &cycles->ctrl;
    wflux_struct       *wf = &cycles->wf;
    wstate_struct      *ws = &cycles->ws;
    phystate_struct    *phys = &cycles->phys;
    forc_struct        *forc = &cycles->forc;
    soil_struct        *soil = &cycles->soil;

    t = (cstep + 1) * ctrl->stepsize;

    wf->precipitation = forc->precipitation[cstep];
    wf->soil_evap = forc->soil_evap[cstep];
    for (kz = 0; kz < NSOIL; kz++)
    {
        total_depth += phys->soil_depth[kz];
    }
    for (kz = 0; kz < NSOIL; kz++)
    {
        wf->uptake[kz] = forc->uptake[cstep] * phys->soil_depth[kz] / total_depth;
    }

    if (wf->precipitation + surface_water / ctrl->stepsize == 0.0)
    {
        wf->infil = 0.0;
        phys->wetting_front_depth = 0.0;
        for (kz = 0; kz < NSOIL; kz++)
        {
            phys->smc_below[kz] = ws->smc[kz];
            phys->smc_wf[kz] = 0.0;
        }
    }
    else
    {
        phys->wetting_front_depth = (phys->wetting_front_depth == 0.0) ?
            DryDepth(ctrl->stepsize, soil, phys, ws->smc[0]) : phys->wetting_front_depth;

        WettingFrontPosition(phys, &n, wetting_front_thickness);

        infiltration_capacity = InfiltrationCapacity(soil, phys, n, wetting_front_thickness);

        if (wf->precipitation + surface_water / ctrl->stepsize > infiltration_capacity)
        {
            wf->infil = infiltration_capacity;
            surface_water += (wf->precipitation - infiltration_capacity) * ctrl->stepsize;
        }
        else
        {
            wf->infil = wf->precipitation + surface_water / ctrl->stepsize;
            surface_water = 0.0;
        }
    }

    wf->runoff = surface_water;
    surface_water = 0.0;

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

    if (phys->wetting_front_depth > 0.0)
    {
        for (kz = 0; kz <= n; kz++)
        {
            if (kz < n)
            {
                phys->smc_wf[kz] = ws->smc[kz];
                phys->smc_below[kz] = 0.0;
            }
            else if (kz == n)
            {
                phys->smc_wf[kz] = SoilMoistureAtWettingFront(kz, wetting_front_thickness[kz], phys->soil_depth[kz], phys->smc_below[kz], ws->smc[kz], soil->porosity[kz]);
            }
            else
            {
                phys->smc_wf[kz] = 0.0;
                phys->smc_below[kz] = ws->smc[kz];
            }
        }

        phys->wetting_front_depth += WettingFrontAdvanceRate(n, phys->wetting_front_depth, wetting_front_thickness, phys->smc_below[n], phys->smc_wf, soil) * ctrl->stepsize;
        phys->wetting_front_depth = MAX(phys->wetting_front_depth, 0.0);
    }
}
