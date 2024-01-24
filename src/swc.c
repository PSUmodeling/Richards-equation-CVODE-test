#include "cycles.h"
//#define TEST

void SWC(int cstep, cycles_struct *cycles, void *cvode_mem, N_Vector CV_Y)
{
    double              t;
    int                 kx;
    forcing_struct     *forcing = &cycles->forcing;

    t = (cstep + 1) * cycles->control.stepsize;

    for (kx = 0; kx < number_of_columns; kx++)
    {
        int             kz;
        double          total_depth = 0.0;
        wflux_struct   *wf = &cycles->grid[kx].wf;
        phystate_struct *phys = &cycles->grid[kx].phys;

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

    // After solving the ODE, update the state variables
    UpdateStateVariables(CV_Y, cycles);

    // Get the solution, and constrain it to be between 0.02 and porosity
    for (kx = 0; kx < number_of_columns; kx++)
    {
        int             kz;
        double          wplus = 0.0;
        double          runoff3 = 0.0;
        wstate_struct  *ws = &cycles->grid[kx].ws;
        soil_struct    *soil = &cycles->grid[kx].soil;
        phystate_struct *phys = &cycles->grid[kx].phys;

        for (kz = 0; kz < number_of_layers; kz++)
        {
            ws->smc[kz] += wplus / phys->soil_depth[kz];

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

int Ode(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *cycles_data)
{
    int             kx;
    int             kkz;
    cycles_struct  *cycles = (cycles_struct *)cycles_data;

    UpdateStateVariables(CV_Y, cycles);

    ZeroFluxes(CV_Ydot, cycles);

    for (kx = 0; kx < number_of_columns; kx++)
    {
        int             kz;
        soil_struct    *soil = &cycles->grid[kx].soil;
        phystate_struct *phys = &cycles->grid[kx].phys;
        wstate_struct  *ws = &cycles->grid[kx].ws;

        for (kz = 0; kz < number_of_layers; kz++)
        {
            phys->retention_capacity[kz] = RetentionCapacity(soil->porosity[kz], soil->air_entry_potential[kz], soil->b[kz], ws->potential[kz]);
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

    // Calculate later fluxes from boundary column to channel
    for (kkz = 0; kkz < number_of_layers; kkz++)
    {
        double          lateral_flux;

        lateral_flux = SubsurfaceToChannel(kkz, &cycles->channel, &cycles->grid[number_of_columns - 1]);
        cycles->grid[number_of_columns - 1].wf.lateral[kkz] -= lateral_flux / cycles->grid[number_of_columns - 1].phys.width;
        cycles->channel.wf.subsurface += lateral_flux / cycles->channel.phys.width;
    }

    cycles->channel.wf.discharge = Discharge(&cycles->channel) / cycles->channel.phys.width / cycles->channel.phys.length;
    //printf("stage = %lf, subsurface = %le, discharge = %le\n", cycles->channel.ws.stage, cycles->channel.wf.subsurface, cycles->channel.wf.discharge);

    for (kx = 0; kx < number_of_columns; kx++)
    {
        int             kz;
        double          wcnd, wcnd2;
        double          dpsidz, dpsidz2;
        soil_struct    *soil = &cycles->grid[kx].soil;
        phystate_struct *phys = &cycles->grid[kx].phys;
        wstate_struct  *ws = &cycles->grid[kx].ws;
        wflux_struct   *wf = &cycles->grid[kx].wf;

        // Top layer
        wcnd = WaterConductivity(ws->smc[0], soil->porosity[0], soil->ksatv[0], soil->b[0]);

        dpsidz = (ws->potential[0] - ws->potential[1]) / (0.5 * (phys->soil_depth[0] + phys->soil_depth[1]));

#ifdef TEST
        wf->infiltration = 0.0;
        wf->soil_evaporation = 0.0;
        wf->uptake[0] = 0.0;
#endif

        NV_Ith_S(CV_Ydot, INDEX(kx, 0)) = (-wcnd * dpsidz - wcnd + wf->infiltration - wf->soil_evaporation - wf->uptake[0] + wf->lateral[0]) / phys->soil_depth[0] /
            RetentionCapacity(soil->porosity[0], soil->air_entry_potential[0], soil->b[0], ws->potential[0]);

        // All other layers
        for (kz = 1; kz < number_of_layers; kz++)
        {
            // Retrieve the soil water diffusivity and hydraulic conductivity for this layer
            wcnd2 = WaterConductivity(ws->smc[kz], soil->porosity[kz], soil->ksatv[kz], soil->b[kz]);

            dpsidz2 = (kz == number_of_layers - 1) ?
                0.0 : (ws->potential[kz] - ws->potential[kz + 1]) / (0.5 * (phys->soil_depth[kz] + phys->soil_depth[kz + 1]));

            // For test purposes, drainage from the bottom is still allowed. Should be set to no flow at the bottom.
#ifdef TEST
            wf->uptake[kz] = 0.0;
#endif
            NV_Ith_S(CV_Ydot, INDEX(kx, kz)) = (-wcnd2 * dpsidz2 - ((kz == number_of_layers - 1) ? 0.0 : 1.0) * wcnd2 + wcnd * dpsidz + wcnd - wf->uptake[kz] + wf->lateral[kz]) / phys->soil_depth[kz] /
                RetentionCapacity(soil->porosity[kz], soil->air_entry_potential[kz], soil->b[kz], ws->potential[kz]);

            if (kz != number_of_layers - 1)
            {
                wcnd = wcnd2;
                dpsidz = dpsidz2;
            }
        }
    }

    NV_Ith_S(CV_Ydot, CHANNEL) = cycles->channel.wf.subsurface - cycles->channel.wf.discharge;

    return 0;
}

double LateralFlux(double width, double width2, double z, double z2, double wcnd, double wcnd2, double psi, double psi2, double area)
{
    double          kh;

    kh = (width * wcnd + width2 * wcnd2) / (width + width2);

    return kh * ((psi + z) - (psi2 + z2)) / (0.5 * width + 0.5 * width2) * area;
}

double SubsurfaceToChannel(int layer, const channel_struct *channel, const grid_struct *grid)
{
    double          contact_height = 0.0;
    double          h_diff;
    double          wcnd;

    h_diff =  (grid->ws.potential[layer] + grid->phys.zsoil[layer] + 0.5 * grid->phys.soil_depth[layer]) - (channel->ws.stage + channel->phys.bed_elevation);


    if (h_diff < 0.0 && channel->ws.stage <= 0.0)
    {
        return 0.0;
    }

    if (h_diff > 0.0)   // From soil to channel
    {
        if (grid->phys.zsoil[layer] >= channel->phys.bed_elevation)            // Layer completely above channel bed
        {
            contact_height = grid->phys.soil_depth[layer];
        }
        else if (grid->phys.zsoil[layer] + grid->phys.soil_depth[layer] >= channel->phys.bed_elevation)   // Layer partially above channel bed
        {
            contact_height = grid->phys.zsoil[layer] + grid->phys.soil_depth[layer] - channel->phys.bed_elevation;
        }
        else                                                                // Layer completely below channel bed
        {
            return 0.0;
        }
    }
    else                // From channel to soil
    {
        if (grid->phys.zsoil[layer] >= channel->ws.stage + channel->phys.bed_elevation)          // Layer completely above channel stage
        {
            return 0.0;
        }
        else if (grid->phys.zsoil[layer] + grid->phys.soil_depth[layer]<= channel->phys.bed_elevation)      // Layer completely below channel
        {
            return 0.0;
        }
        else if (grid->phys.zsoil[layer] >= channel->phys.bed_elevation)                       // Layer completely above channel bed
        {
            contact_height = MIN(grid->phys.soil_depth[layer], channel->ws.stage + channel->phys.bed_elevation - grid->phys.zsoil[layer]);
        }
        else if (grid->phys.zsoil[layer] + grid->phys.soil_depth[layer]>= channel->phys.bed_elevation)
        {
            contact_height = MIN(channel->ws.stage, grid->phys.soil_depth[layer] + grid->phys.soil_depth[layer] - channel->phys.bed_elevation);
        }
        else
        {
            printf("Not considered.\n");
        }
    }

    //printf("layer = %d, zsoil = %lf, bed_elevation = %lf, stage = %lf, contact_height = %lf, h_diff = %lf\n", layer, grid->phys.zsoil[layer], channel->phys.bed_elevation, channel->ws.stage, contact_height, h_diff);

    wcnd = WaterConductivity(grid->ws.smc[layer], grid->soil.porosity[layer], grid->soil.ksath[layer], grid->soil.b[layer]);

    return wcnd * contact_height * h_diff / (0.5 * channel->phys.width + 0.5 * grid->phys.width);
}

// Calculates channel discharge in m3/s using Manning's equation
double Discharge(channel_struct *channel)
{
    if (channel->ws.stage <= 0.0)
    {
        return 0.0;
    }

    //printf("stage = %lf, q = %lf\n", channel->ws.stage, sqrt(channel->phys.slope) * (channel->phys.width * channel->ws.stage) *
    //    pow(channel->phys.width * channel->ws.stage / (channel->phys.width + 2 * channel->ws.stage), 2.0 / 3.0) / channel->phys.roughness);
    return sqrt(channel->phys.slope) * (channel->phys.width * channel->ws.stage) *
        pow(channel->phys.width * channel->ws.stage / (channel->phys.width + 2 * channel->ws.stage), 2.0 / 3.0) / channel->phys.roughness;
}

void UpdateStateVariables(N_Vector CV_Y, cycles_struct *cycles)
{
    int             kx;

    for (kx = 0; kx < number_of_columns; kx++)
    {
        int             kz;
        soil_struct    *soil = &cycles->grid[kx].soil;
        wstate_struct  *ws = &cycles->grid[kx].ws;

        for (kz = 0; kz < number_of_layers; kz++)
        {
            ws->potential[kz] = NV_Ith_S(CV_Y, INDEX(kx, kz));
            ws->smc[kz] = SoilWaterContent(soil->porosity[kz], soil->air_entry_potential[kz], soil->b[kz], ws->potential[kz]);
        }
    }

    cycles->channel.ws.stage = MAX(NV_Ith_S(CV_Y, CHANNEL), 0.0);
}

void ZeroFluxes(N_Vector CV_Ydot, cycles_struct *cycles)
{
    int             kx;

    for (kx = 0; kx < number_of_columns; kx++)
    {
        int             kz;
        wflux_struct   *wf = &cycles->grid[kx].wf;

        for (kz = 0; kz < number_of_layers; kz++)
        {
            wf->lateral[kz] = 0.0;
            NV_Ith_S(CV_Ydot, INDEX(kx, kz)) = 0.0;
        }
    }

    cycles->channel.wf.discharge = 0.0;
    cycles->channel.wf.subsurface = 0.0;
    NV_Ith_S(CV_Ydot, CHANNEL) = 0.0;
}

double TotalWaterStorage(grid_struct *grid, channel_struct *channel)
{
    int             kx;
    double          total_storage = 0.0;
    double          total_width = 0.0;

    for (kx = 0; kx < number_of_columns; kx++)
    {
        int             kz;
        wstate_struct  *ws = &grid[kx].ws;

        for (kz = 0; kz < number_of_layers; kz++)
        {
            total_storage += ws->smc[kz] * grid[kx].phys.soil_depth[kz] * grid[kx].phys.width;
        }
        total_width += grid[kx].phys.width;
    }

    total_storage += channel->ws.stage * channel->phys.width;

    return total_storage / total_width;
}

double TotalWaterFlux(grid_struct *grid, channel_struct *channel)
{
    int             kx;
    double          total_flux = 0.0;
    double          total_width = 0.0;

    for (kx = 0; kx < number_of_columns; kx++)
    {
        int             kz;
        wflux_struct   *wf = &grid[kx].wf;

        total_flux += (wf->infiltration - wf->soil_evaporation) * grid[kx].phys.width;

        for (kz = 0; kz < number_of_layers; kz++)
        {
            total_flux += wf->uptake[kz] * grid[kx].phys.width;
        }
        total_width += grid[kx].phys.width;
    }

    total_flux += channel->wf.discharge * channel->phys.width;

    return total_flux / total_width;
}
