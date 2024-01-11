#include "cycles.h"

int Ode(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *cycles_data)
{
    int             kx;
    double         *y;
    double         *dy;
    cycles_struct  *cycles;

    y = NV_DATA_S(CV_Y);
    dy = NV_DATA_S(CV_Ydot);
    cycles = (cycles_struct *)cycles_data;


    for (kx = 0; kx < number_of_columns; kx++)
    {
        int             kz;
        double          wdf, wcnd;
        double          wdf2, wcnd2;
        double          dpsidz, dpsidz2;
        soil_struct    *soil = &cycles->grid[kx].soil;
        phystate_struct *phys = &cycles->grid[kx].phys;
        wstate_struct  *ws = &cycles->grid[kx].ws;
        wflux_struct   *wf = &cycles->grid[kx].wf;

        for (kz = 0; kz < number_of_layers; kz++)
        {
            ws->potential[kz] = y[INDEX(kx, kz)];
            ws->smc[kz] = SoilWaterContent(soil->porosity[kz], soil->air_entry_potential[kz], soil->b[kz], ws->potential[kz]);
            phys->retention_capacity[kz] = RetentionCapacity(soil->porosity[kz], soil->air_entry_potential[kz], soil->b[kz], ws->potential[kz]);
        }

        // Top layer
        WDfCnd(ws->smc[0], soil->porosity[0], soil->dsat[0], soil->ksatv[0], soil->b[0], &wdf, &wcnd);

        dpsidz = (ws->potential[0] - ws->potential[1]) / (0.5 * (phys->soil_depth[0] + phys->soil_depth[1]));

        dy[INDEX(kx, 0)] = (-wcnd * dpsidz - wcnd + wf->infiltration - wf->soil_evaporation - wf->uptake[0]) / phys->soil_depth[0] /
            RetentionCapacity(soil->porosity[0], soil->air_entry_potential[0], soil->b[0], ws->potential[0]);

        // All other layers
        for (kz = 1; kz < number_of_layers; kz++)
        {
            // Retrieve the soil water diffusivity and hydraulic conductivity for this layer
            WDfCnd(ws->smc[kz], soil->porosity[kz], soil->dsat[kz], soil->ksatv[kz], soil->b[kz], &wdf2, &wcnd2);

            dpsidz2 = (kz == number_of_layers - 1) ?
                0.0 : (ws->potential[kz] - ws->potential[kz + 1]) / (0.5 * (phys->soil_depth[kz] + phys->soil_depth[kz + 1]));

            dy[INDEX(kx, kz)] = (-wcnd2 * dpsidz2 - ((kz == number_of_layers) ? 0.1 : 1.0) * wcnd2 + wcnd * dpsidz + wcnd - wf->uptake[kz]) / phys->soil_depth[kz] /
                RetentionCapacity(soil->porosity[kz], soil->air_entry_potential[kz], soil->b[kz], ws->potential[kz]);

            if (kz != number_of_layers - 1)
            {
                wdf = wdf2;
                wcnd = wcnd2;
                dpsidz = dpsidz2;
            }
        }
    }

    return 0;
}
