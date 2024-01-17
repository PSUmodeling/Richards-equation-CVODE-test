#include "cycles.h"

void ReadHydrologicalForcing(forcing_struct *forcing)
{
    char            cmdstr[MAXSTRING];
    int             lno = 0;
    int             k;
    FILE           *fp;

    fp = cycles_fopen("HYDRO.TXT", "r");

    // Read weather file
    FindLine(fp, "BOF", &lno, "HYDRO.TXT");

    // Read in site information and count number of weather records
    NextLine(fp, cmdstr, &lno);

    for (k = 0; k < NUMBER_OF_STEPS; k++)
    {
        double          precip, runoff;

        NextLine(fp, cmdstr, &lno);
        sscanf(cmdstr, "%*d %*d %lf %*lf %*lf %lf %*lf %*lf %*lf %*lf %lf %lf",
            &precip, &runoff, &forcing->soil_evaporation[k], &forcing->uptake[k]);
        forcing->infiltration[k] = (precip - runoff) / 1000.0 / 1800.0;   // from mm to m/s
        forcing->soil_evaporation[k] /=  1000.0 * 1800.0;
        forcing->uptake[k] /= 1000.0 * 1800.0;
    }

    fclose(fp);
}
