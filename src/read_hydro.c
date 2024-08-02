#include "cycles.h"

void ReadHydro(forc_struct *forc)
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

    for (k = 0; k < NSTEPS; k++)
    {
        double          precip, runoff;

        NextLine(fp, cmdstr, &lno);
        sscanf(cmdstr, "%*d %*d %lf %*lf %*lf %lf %*lf %*lf %*lf %*lf %lf %lf",
            &precip, &runoff, &forc->soil_evap[k], &forc->uptake[k]);
        forc->precipitation[k] = precip / 1000.0 / 1800.0;   // from mm to m/s
        forc->infil[k] = (precip - runoff) / 1000.0 / 1800.0;   // from mm to m/s
        forc->soil_evap[k] /=  1000.0 * 1800.0;
        forc->uptake[k] /= 1000.0 * 1800.0;
    }

    fclose(fp);
}
