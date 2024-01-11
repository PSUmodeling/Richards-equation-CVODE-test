#include "cycles.h"

int             verbose_mode;
int             debug_mode;
int             number_of_columns;
int             number_of_layers;

int main(int argc, char *argv[])
{
    N_Vector        CV_Y;
    void           *cvode_mem;
    SUNLinearSolver sun_ls;
    cycles_struct   cycles;
    int             kstep;
    int             kx;
    FILE           *smc_fp;
    FILE           *wp_fp;

#if defined(unix) || defined(__unix__) || defined(__unix)
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

    // Read domain input
    ReadDomain("input/test/domain.txt", &cycles.control);

    // Create grids
    cycles.grid = (grid_struct *)malloc(number_of_columns * sizeof(grid_struct));
    for (kx = 0; kx < number_of_columns; kx++)
    {
        cycles.grid[kx].soil.porosity = (double *)malloc(number_of_layers * sizeof(double));
        cycles.grid[kx].soil.air_entry_potential = (double *)malloc(number_of_layers * sizeof(double));
        cycles.grid[kx].soil.b = (double *)malloc(number_of_layers * sizeof(double));
        cycles.grid[kx].soil.ksath = (double *)malloc(number_of_layers * sizeof(double));
        cycles.grid[kx].soil.ksatv = (double *)malloc(number_of_layers * sizeof(double));
        cycles.grid[kx].soil.dsat = (double *)malloc(number_of_layers * sizeof(double));

        cycles.grid[kx].phys.zsoil = (double *)malloc(number_of_layers * sizeof(double));
        cycles.grid[kx].phys.soil_depth = (double *)malloc(number_of_layers * sizeof(double));
        cycles.grid[kx].phys.retention_capacity = (double *)malloc(number_of_layers * sizeof(double));

        cycles.grid[kx].ws.smc = (double *)malloc(number_of_layers * sizeof(double));
        cycles.grid[kx].ws.potential = (double *)malloc(number_of_layers * sizeof(double));

        cycles.grid[kx].wf.uptake = (double *)malloc(number_of_layers * sizeof(double));
        cycles.grid[kx].wf.lateral = (double *)malloc(number_of_layers * sizeof(double));
    }

    // Read other input data
    ReadTopography("input/test/topography.txt", &cycles);
    ReadSoil("input/test/b.txt", B, &cycles);
    ReadSoil("input/test/ksath.txt", KSATH, &cycles);
    ReadSoil("input/test/ksatv.txt", KSATV, &cycles);
    ReadSoil("input/test/porosity.txt", POROSITY, &cycles);
    ReadSoil("input/test/psi_sat.txt", AIR_ENTRY_POTENTIAL, &cycles);
    ReadSoil("input/test/initial_condition.txt", INITIAL_CONDITION, &cycles);

    // Create CVode state variable array
    CV_Y = N_VNew_Serial(number_of_columns * number_of_layers);
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
    ReadHydro(&cycles.forcing);

    // Open output file and write header
    smc_fp = fopen("smc.txt", "w");
    wp_fp = fopen("wp.txt", "w");

    // Run model
    for (kstep = 0; kstep < NSTEPS; kstep++)
    {
        // Run soil moisture time step
        SWC(kstep, &cycles, cvode_mem, CV_Y);

        // Write output
        fprintf(smc_fp, "%d", kstep);
        fprintf(wp_fp, "%d", kstep);
        for (kx = 0; kx < number_of_columns; kx++)
        {
            int             kz;
            wstate_struct  *ws = &cycles.grid[kx].ws;

            for (kz = 0; kz < number_of_layers; kz++)
            {
                fprintf(smc_fp, ",%.3lf", ws->smc[kz]);
                fprintf(wp_fp, ",%.3lf", ws->potential[kz]);
            }
        }
        fprintf(smc_fp, "\n");
        fprintf(wp_fp, "\n");
    }

    //fclose(fp);

    // Free memory
    N_VDestroy(CV_Y);

    // Free integrator memory
    CVodeFree(&cvode_mem);
    SUNLinSolFree(sun_ls);

    return EXIT_SUCCESS;
}
