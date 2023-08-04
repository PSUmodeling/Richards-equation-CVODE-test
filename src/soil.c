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
