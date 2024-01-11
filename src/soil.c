#include "cycles.h"

void WDfCnd(double smc, double smcmax, double dwsat, double dksat, double bexp, double *wdf, double *wcnd)
{
    double          factr2;
    double          expon;

    factr2 = MIN(1.0, MAX(0.0, smc / smcmax));
    expon = bexp + 2.0;
    *wdf = dwsat * pow(factr2, expon);

    expon = 2.0 * bexp + 3.0;
    *wcnd = dksat * pow(factr2, expon);

    return;
}

double SoilWaterPotential(double sat_wc, double air_entry_pot, double campbell_b, double wc)
{
    double          factr2;
    double          expon;

    factr2 = MIN(1.0, MAX(0.0, wc / sat_wc));
    expon = -campbell_b;

    return -air_entry_pot * pow(factr2, expon);
}

double SoilWaterContent(double sat_wc, double air_entry_pot, double campbell_b, double pot)
{
    double          factr2;
    double          expon;

    if (pot >= -air_entry_pot) return sat_wc;

    expon = -1.0 / campbell_b;
    factr2 = -pot / air_entry_pot;

    return sat_wc * pow(factr2, expon);
}

double RetentionCapacity(double sat_wc, double air_entry_pot, double campbell_b, double pot)
{
    double          factr2;
    double          expon;

    if (pot >= -air_entry_pot) return 1.0E-5;

    expon = -1.0 / campbell_b;
    factr2 = -pot / air_entry_pot;

    return sat_wc * expon * pow(factr2, expon) / pot;
}
