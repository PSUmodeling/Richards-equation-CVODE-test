#include "cycles.h"

void GenerateGrids(cycles_struct *cycles)
{
    int             kx;

    cycles->grid = (grid_struct *)malloc(number_of_columns * sizeof(grid_struct));

    for (kx = 0; kx < number_of_columns; kx++)
    {
        cycles->grid[kx].soil.porosity = (double *)malloc(number_of_layers * sizeof(double));
        cycles->grid[kx].soil.air_entry_potential = (double *)malloc(number_of_layers * sizeof(double));
        cycles->grid[kx].soil.b = (double *)malloc(number_of_layers * sizeof(double));
        cycles->grid[kx].soil.ksath = (double *)malloc(number_of_layers * sizeof(double));
        cycles->grid[kx].soil.ksatv = (double *)malloc(number_of_layers * sizeof(double));

        cycles->grid[kx].phys.zsoil = (double *)malloc(number_of_layers * sizeof(double));
        cycles->grid[kx].phys.soil_depth = (double *)malloc(number_of_layers * sizeof(double));
        cycles->grid[kx].phys.retention_capacity = (double *)malloc(number_of_layers * sizeof(double));

        cycles->grid[kx].ws.smc = (double *)malloc(number_of_layers * sizeof(double));
        cycles->grid[kx].ws.potential = (double *)malloc(number_of_layers * sizeof(double));

        cycles->grid[kx].wf.uptake = (double *)malloc(number_of_layers * sizeof(double));
        cycles->grid[kx].wf.lateral = (double *)malloc(number_of_layers * sizeof(double));
    }
}

void Initialize(cycles_struct *cycles, N_Vector CV_Y, void **cvode_mem)
{
    int             kx;
    control_struct *control = &cycles->control;

    control->stepsize = STEPSIZE;
    control->solver_stepsize = SOLVER_STEPSIZE;

    for (kx = 0; kx < number_of_columns; kx++)
    {
        int             kz;
        soil_struct    *soil;
        wstate_struct  *ws;

        soil = &cycles->grid[kx].soil;
        ws = &cycles->grid[kx].ws;

        for (kz = 0; kz < number_of_layers; kz++)
        {
            // Initialize CVode variables for soil columns
            ws->potential[kz] = SoilWaterPotential(soil->porosity[kz], soil->air_entry_potential[kz], soil->b[kz], ws->smc[kz]);
            NV_Ith_S(CV_Y, INDEX(kx, kz)) = ws->potential[kz];
        }
    }

    // Initialize CVode variables for channel
    cycles->channel.ws.stage = 0.0;
    NV_Ith_S(CV_Y, CHANNEL) = 0.0;

    // Allocate memory for solver
    *cvode_mem = CVodeCreate(CV_BDF);
    if (*cvode_mem == NULL)
    {
        printf("Error in allocating memory for solver.\n");
        exit(EXIT_FAILURE);
    }

    return;
}
