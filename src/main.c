#include "cycles.h"

int             verbose_mode;
int             debug_mode;

int main(int argc, char *argv[])
{
    N_Vector        CV_Y;
    void           *cvode_mem;
    SUNLinearSolver sun_ls;
    cycles_struct   cycles;
    int             kstep;
    FILE           *fp;

#if defined(unix) || defined(__unix__) || defined(__unix)
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

    // Create CVode state variable array
    CV_Y = N_VNew_Serial(NSOIL);
    if (CV_Y == NULL)
    {
        printf("Error creating CVode state variable vector.\n");
        exit(EXIT_FAILURE);
    }

    // Initialize model structure and CVode state variables
    Initialize(&cycles, CV_Y, &cvode_mem);

    // Create linear solver object and set optional CVode parameters
    SetCVodeParam(cvode_mem, &sun_ls, CV_Y, &cycles);

    // Read hydrological forcing data
    ReadHydro(&cycles.forc);

    // Open output file and write header
    fp = fopen("cycles.txt", "w");
    fprintf(fp, "step,smc0,smc1,smc2,smc3\n");

    // Run model
    for (kstep = 0; kstep < NSTEPS; kstep++)
    {
        int             kz;

        // Run soil moisture time step
        SWC(kstep, &cycles, cvode_mem, CV_Y);

        // Write output
        fprintf(fp, "%d", kstep);
        for (kz = 0; kz < NSOIL; kz++)
        {
            fprintf(fp, ",%.3lf", cycles.ws.smc[kz]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);

    // Free memory
    N_VDestroy(CV_Y);

    // Free integrator memory
    CVodeFree(&cvode_mem);
    SUNLinSolFree(sun_ls);

    return EXIT_SUCCESS;
}
