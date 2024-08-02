#include "cycles.h"

int Ode(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *cycles_data)
{
    int             kz;
    double         *y;
    double         *dy;
    double          wdf, wcnd;
    double          wdf2, wcnd2;
    double          dpsidz, dpsidz2;
    double          rc[NSOIL];
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
        ws->potential[kz] = y[kz];
        ws->smc[kz] = SoilWaterContent(soil->porosity[kz], soil->air_entry_pot[kz], soil->b[kz], ws->potential[kz]);
        rc[kz] = RetentionCapacity(soil->porosity[kz], soil->air_entry_pot[kz], soil->b[kz], ws->potential[kz]);
    }

    // Top layer
    WDfCnd(ws->smc[0], soil->porosity[0], soil->dsat[0], soil->ksat[0], soil->b[0], &wdf, &wcnd);

    dpsidz = (ws->potential[0] - ws->potential[1]) / (0.5 * (phys->soil_depth[0] + phys->soil_depth[1]));

    dy[0] = (-wcnd * dpsidz - wcnd + wf->infil - wf->soil_evap - wf->uptake[0]) / phys->soil_depth[0] /
        RetentionCapacity(soil->porosity[0], soil->air_entry_pot[0], soil->b[0], ws->potential[0]);

    // All other layers
    for (kz = 1; kz < NSOIL; kz++)
    {
        double          slopx;

        // Slope of bottom layer is introduced
        slopx = (kz == NSOIL - 1) ? soil->slope : 1.0;

        // Retrieve the soil water diffusivity and hydraulic conductivity for this layer
        WDfCnd(ws->smc[kz], soil->porosity[kz], soil->dsat[kz], soil->ksat[kz], soil->b[kz], &wdf2, &wcnd2);

        dpsidz2 = (kz == NSOIL - 1) ?
            0.0 : (ws->potential[kz] - ws->potential[kz + 1]) / (0.5 * (phys->soil_depth[kz] + phys->soil_depth[kz + 1]));

        dy[kz] = (-wcnd2 * dpsidz2 - slopx * wcnd2 + wcnd * dpsidz + wcnd - wf->uptake[kz]) / phys->soil_depth[kz] /
            RetentionCapacity(soil->porosity[kz], soil->air_entry_pot[kz], soil->b[kz], ws->potential[kz]);

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

double CapillaryDrive(double bexp, double air_entry_potential, double swc_b, double swc)
{
    return bexp / (bexp + 3) * air_entry_potential * (pow(swc, bexp + 3) - pow(swc_b, bexp + 3));
}

double DryDepth(double dt, const soil_struct *soil, const phystate_struct *phys, double swc_b)
{
    double          tau;
    double          g;

    tau = dt * soil->ksat[0] / (soil->porosity[0] - swc_b);
    g = CapillaryDrive(soil->b[0], soil->air_entry_pot[0], swc_b, soil->porosity[0]);

    return MIN(0.5 * (tau + sqrt(tau * tau + 4 * g * tau)), phys->soil_depth[0]);
}

void WettingFrontPosition(const phystate_struct *phys,  int *layer, double wetting_front_thickness[])
{
    int             kz;
    double          ztop;

    *layer = NSOIL - 1;
    for (kz = 0; kz < NSOIL; kz++)
    {
        if (phys->wetting_front_depth <= -phys->zsoil[kz])
        {
            *layer = kz;
            break;
        }
    }

    for (kz = 0; kz < NSOIL; kz++)
    {
        wetting_front_thickness[kz] = 0.0;
    }

    for (kz = 0; kz <= *layer; kz++)
    {
        ztop = (kz == 0) ? 0.0 : phys->zsoil[kz - 1];
        wetting_front_thickness[kz] = MIN(phys->wetting_front_depth + ztop, phys->soil_depth[kz]);
        wetting_front_thickness[kz] = MAX(wetting_front_thickness[kz], 0.0);
    }
}

double SoilMoistureAtWettingFront(int layer, double wetting_front_thickness, double soil_depth, double swc_b, double swc, double smcmax)
{
    return MIN(MAX((swc * soil_depth - swc_b * (soil_depth - wetting_front_thickness)) / wetting_front_thickness, 0.02), smcmax);
}

double InfiltrationCapacity(const soil_struct *soil, const phystate_struct *phys, int n, const double wetting_front_thickness[])
{
    int             kz;
    double          denom = 0.0;

    for (kz = 0; kz <= n; kz++)
    {
        denom += wetting_front_thickness[kz] / soil->ksat[kz];
    }

    return soil->ksat[n] * CapillaryDrive(soil->b[n], soil->air_entry_pot[n], phys->smc_below[n], soil->porosity[n]) / phys->wetting_front_depth + phys->wetting_front_depth / denom;
}

double CompositeConductivity(int n, double wetting_front_depth, const double wetting_front_thickness[], const double swc[], const soil_struct *soil)
{
    int             kz;
    double          denominator = 0.0;
    double          wdf, wcnd;

    for (kz = 0; kz <= n; kz++)
    {
        WDfCnd(swc[kz], soil->porosity[kz], soil->dsat[kz], soil->ksat[kz], soil->b[kz], &wdf, &wcnd);

        denominator += wetting_front_thickness[kz] / wcnd;
    }

    return wetting_front_depth / denominator;
}

double WettingFrontAdvanceRate(int n, double wetting_front_depth, const double wetting_front_thickness[], double swc_b, const double swc_wf[], const soil_struct *soil)
{
    double          rate = 0.0;
    rate = (swc_wf[n] == swc_b) ?
        0.0 : 1.0 / (swc_wf[n] - swc_b) * (soil->ksat[n] * CapillaryDrive(soil->b[n], soil->air_entry_pot[n], swc_b, swc_wf[n]) / wetting_front_depth + CompositeConductivity(n, wetting_front_depth, wetting_front_thickness, swc_wf, soil));

    printf("z %f, swc_wf - swc_b %f, Ksn %e, Gn %e, Kc %e\n", wetting_front_depth, swc_wf[n] - swc_b, soil->ksat[n], CapillaryDrive(soil->b[n], soil->air_entry_pot[n], swc_b, swc_wf[n]), CompositeConductivity(n, wetting_front_depth, wetting_front_thickness, swc_wf, soil));
    printf("%e\n", rate);

    return rate;
}
