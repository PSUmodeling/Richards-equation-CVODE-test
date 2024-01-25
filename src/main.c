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

#if defined(unix) || defined(__unix__) || defined(__unix)
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

    ReadDomain("input/test/domain.txt", &cycles.control);

    GenerateGrids(&cycles);

    ReadTopography("input/test/topography.txt", &cycles);

    ReadSoil("input/test/b.txt", B, &cycles);
    ReadSoil("input/test/ksath.txt", KSATH, &cycles);
    ReadSoil("input/test/ksatv.txt", KSATV, &cycles);
    ReadSoil("input/test/porosity.txt", POROSITY, &cycles);
    ReadSoil("input/test/psi_sat.txt", AIR_ENTRY_POTENTIAL, &cycles);
    ReadSoil("input/test/initial_condition.txt", INITIAL_CONDITION, &cycles);

    // Create CVode state variable array
    CV_Y = N_VNew_Serial(number_of_columns * number_of_layers + 1);
    if (CV_Y == NULL)
    {
        printf("Error creating CVode state variable vector.\n");
        exit(EXIT_FAILURE);
    }

    // Initialize model structure and CVode state variables
    Initialize(&cycles, CV_Y, &cvode_mem);

    // Create linear solver object and set optional CVode parameters
    SetCVodeParameters(cvode_mem, &sun_ls, CV_Y, &cycles);

    ReadHydrologicalForcing(&cycles.forcing);

    OpenOutputFiles(&cycles.output);

    // Run model
    for (kstep = 0; kstep < NUMBER_OF_STEPS; kstep++)
    {
        // Run soil moisture time step
        SWC(kstep, &cycles, cvode_mem, CV_Y);

        WriteOutputFiles(kstep + 1, cycles.control.stepsize / cycles.control.solver_stepsize, &cycles.output);
    }

    N_VDestroy(CV_Y);
    SUNLinSolFree(sun_ls);
    CVodeFree(&cvode_mem);

    Cleanup(&cycles);

    return EXIT_SUCCESS;
}
