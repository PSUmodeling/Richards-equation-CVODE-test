#include "cycles.h"

int Ode(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *cycles_data)
{
    int             kz;
    double         *y;
    double         *dy;
    double          wdf, wcnd;
    double          wdf2, wcnd2;
    double          dpsidz, dpsidz2;
    double          wpot[NSOIL];
    cycles_struct  *cycles;
    soil_struct    *soil;
    phystate_struct *phys;
    wstate_struct  *ws;
    wflux_struct   *wf;

    y = NV_DATA_S(CV_Y);
    dy = NV_DATA_S(CV_Ydot);
    cycles = (cycles_struct *)cycles_data;
    soil = &cycles->soil;
    phys = &cycles->phys;
    ws = &cycles->ws;
    wf = &cycles->wf;

    for (kz = 0; kz < NSOIL; kz++)
    {
        ws->smc[kz] = y[kz];
        wpot[kz] = SoilWaterPot(soil->porosity[kz], soil->air_entry_pot[kz], soil->b[kz], ws->smc[kz]);
    }

    // Top layer
    WDfCnd(ws->smc[0], soil->porosity[0], soil->dsat[0], soil->ksat[0], soil->b[0], &wdf, &wcnd);

    dpsidz = (wpot[0] - wpot[1]) / (0.5 * (phys->soil_depth[0] + phys->soil_depth[1]));

    dy[0] = (-wcnd * dpsidz - wcnd + wf->infil - wf->soil_evap - wf->uptake[0]) / phys->soil_depth[0];

    // All other layers
    for (kz = 1; kz < NSOIL; kz++)
    {
        double          slopx;

        // Slope of bottom layer is introduced
        slopx = (kz == NSOIL - 1) ? soil->slope : 1.0;

        // Retrieve the soil water diffusivity and hydraulic conductivity for this layer
        WDfCnd(ws->smc[kz], soil->porosity[kz], soil->dsat[kz], soil->ksat[kz], soil->b[kz], &wdf2, &wcnd2);

        dpsidz2 = (kz == NSOIL - 1) ?
            0.0 : (wpot[kz] - wpot[kz + 1]) / (0.5 * (phys->soil_depth[kz] + phys->soil_depth[kz + 1]));

        dy[kz] = (-wcnd2 * dpsidz2 - slopx * wcnd2 + wcnd * dpsidz + wcnd - wf->uptake[kz]) / phys->soil_depth[kz];

        if (kz == NSOIL - 1)
        {
            wf->drainage = slopx * wcnd2;
        }
        else
        {
            wdf = wdf2;
            wcnd = wcnd2;
            dpsidz = dpsidz2;
        }
    }

    return 0;
}
