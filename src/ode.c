#include "cycles.h"
#define TEST

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
        soil_struct    *soil = &cycles->grid[kx].soil;
        phystate_struct *phys = &cycles->grid[kx].phys;
        wstate_struct  *ws = &cycles->grid[kx].ws;
        wflux_struct   *wf = &cycles->grid[kx].wf;

        for (kz = 0; kz < number_of_layers; kz++)
        {
            ws->potential[kz] = y[INDEX(kx, kz)];
            ws->smc[kz] = SoilWaterContent(soil->porosity[kz], soil->air_entry_potential[kz], soil->b[kz], ws->potential[kz]);
            phys->retention_capacity[kz] = RetentionCapacity(soil->porosity[kz], soil->air_entry_potential[kz], soil->b[kz], ws->potential[kz]);

            wf->lateral[kz] = 0.0;
        }
    }

    // Calculate lateral fluxes between columns (from left to right)
    for (kx = 0; kx < number_of_columns - 1; kx++)
    {
        int             kz;
        double          wcnd, wcnd2;
        double          lateral_flux;
        soil_struct    *soil = &cycles->grid[kx].soil;
        soil_struct    *soil2 = &cycles->grid[kx + 1].soil;
        phystate_struct *phys = &cycles->grid[kx].phys;
        phystate_struct *phys2 = &cycles->grid[kx + 1].phys;
        wstate_struct  *ws = &cycles->grid[kx].ws;
        wstate_struct  *ws2 = &cycles->grid[kx + 1].ws;
        wflux_struct   *wf = &cycles->grid[kx].wf;
        wflux_struct   *wf2 = &cycles->grid[kx + 1].wf;

        for (kz = 0; kz < number_of_layers; kz++)
        {
            wcnd = WaterConductivity(ws->smc[kz], soil->porosity[kz], soil->ksath[kz], soil->b[kz]);
            wcnd2 = WaterConductivity(ws2->smc[kz], soil2->porosity[kz], soil2->ksath[kz], soil2->b[kz]);

            lateral_flux = LateralFlux(phys->width, phys2->width, phys->zsoil[kz], phys2->zsoil[kz], wcnd, wcnd2, ws->potential[kz], ws2->potential[kz], phys->soil_depth[kz]);

            wf->lateral[kz] -= lateral_flux / phys->width;
            wf2->lateral[kz] += lateral_flux / phys2->width;
        }
    }

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

        // Top layer
        wdf = WaterDiffusivity(ws->smc[0], soil->porosity[0], soil->dsat[0], soil->b[0]);
        wcnd = WaterConductivity(ws->smc[0], soil->porosity[0], soil->ksatv[0], soil->b[0]);

        dpsidz = (ws->potential[0] - ws->potential[1]) / (0.5 * (phys->soil_depth[0] + phys->soil_depth[1]));

#ifdef TEST
        wf->infiltration = 0.0;
        wf->soil_evaporation = 0.0;
        wf->uptake[0] = 0.0;
#endif

        dy[INDEX(kx, 0)] = (-wcnd * dpsidz - wcnd + wf->infiltration - wf->soil_evaporation - wf->uptake[0] + wf->lateral[0]) / phys->soil_depth[0] /
            RetentionCapacity(soil->porosity[0], soil->air_entry_potential[0], soil->b[0], ws->potential[0]);

        // All other layers
        for (kz = 1; kz < number_of_layers; kz++)
        {
            // Retrieve the soil water diffusivity and hydraulic conductivity for this layer
            wdf2 = WaterDiffusivity(ws->smc[kz], soil->porosity[kz], soil->dsat[kz], soil->b[kz]);
            wcnd2 = WaterConductivity(ws->smc[kz], soil->porosity[kz], soil->ksatv[kz], soil->b[kz]);

            dpsidz2 = (kz == number_of_layers - 1) ?
                0.0 : (ws->potential[kz] - ws->potential[kz + 1]) / (0.5 * (phys->soil_depth[kz] + phys->soil_depth[kz + 1]));

            // For test purposes, drainage from the bottom is still allowed. Should be set to no flow at the bottom.
#ifdef TEST
            wf->uptake[kz] = 0.0;
#endif
            dy[INDEX(kx, kz)] = (-wcnd2 * dpsidz2 - ((kz == number_of_layers - 1) ? 0.0 : 1.0) * wcnd2 + wcnd * dpsidz + wcnd - wf->uptake[kz] + wf->lateral[kz]) / phys->soil_depth[kz] /
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

double LateralFlux(double width, double width2, double z, double z2, double wcnd, double wcnd2, double psi, double psi2, double area)
{
    double          kh;

    kh = (width * wcnd + width2 * wcnd2) / (width + width2);

    return kh * ((psi + z) - (psi2 + z2)) / (0.5 * width + 0.5 * width2) * area;
}
