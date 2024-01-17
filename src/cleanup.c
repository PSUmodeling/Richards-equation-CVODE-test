#include "cycles.h"

void Cleanup(cycles_struct *cycles)
{
    int             kx;

    for (kx = 0; kx < number_of_columns; kx++)
    {
        free(cycles->grid[kx].soil.porosity);
        free(cycles->grid[kx].soil.air_entry_potential);
        free(cycles->grid[kx].soil.b);
        free(cycles->grid[kx].soil.ksath);
        free(cycles->grid[kx].soil.ksatv);

        free(cycles->grid[kx].phys.zsoil);
        free(cycles->grid[kx].phys.soil_depth);
        free(cycles->grid[kx].phys.retention_capacity);

        free(cycles->grid[kx].ws.smc);
        free(cycles->grid[kx].ws.potential);

        free(cycles->grid[kx].wf.uptake);
        free(cycles->grid[kx].wf.lateral);
    }

    free(cycles->grid);

    CloseOutputFiles(&cycles->output);
}
