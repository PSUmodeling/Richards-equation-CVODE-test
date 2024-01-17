#include "cycles.h"

void ReadDomain(const char file_name[], control_struct *control)
{
    char            cmdstr[MAXSTRING];
    int             lno = 0;
    FILE           *fp;

    fp = cycles_fopen(file_name, "r");

    // Read number of soil columns
    NextLine(fp, cmdstr, &lno);
    sscanf(cmdstr, "%d", &number_of_columns);

    // Read number of soil layers
    NextLine(fp, cmdstr, &lno);
    sscanf(cmdstr, "%d", &number_of_layers);

    // Read layer bed_elevation
    NextLine(fp, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &control->layer_depth);

    fclose(fp);
}

void ReadTopography(const char file_name[], cycles_struct *cycles)
{
    int             kx;
    int             kz;
    char            cmdstr[MAXSTRING];
    int             lno = 0;
    double         *elevation;
    double         *width;
    FILE           *fp;

    fp = cycles_fopen(file_name, "r");

    elevation = (double *)malloc(number_of_columns * sizeof(double));
    width = (double *)malloc(number_of_columns * sizeof(double));

    // Read elevation
    NextLine(fp, cmdstr, &lno);
    ReadMultipleValues(cmdstr, number_of_columns, 'd', elevation);

    for (kx = 0; kx < number_of_columns; kx++)
    {
        for (kz = 0; kz < number_of_layers; kz++)
        {
            cycles->grid[kx].phys.zsoil[kz] = (kz == 0) ?
                elevation[kx] - cycles->control.layer_depth : cycles->grid[kx].phys.zsoil[kz - 1] - cycles->control.layer_depth;
            cycles->grid[kx].phys.soil_depth[kz] = cycles->control.layer_depth;
        }
    }

    // Read width
    NextLine(fp, cmdstr, &lno);
    ReadMultipleValues(cmdstr, number_of_columns, 'd', width);

    for (kx = 0; kx < number_of_columns; kx++)
    {
        for (kz = 0; kz < number_of_layers; kz++)
        {
            cycles->grid[kx].phys.width = width[kx];
        }
    }

    NextLine(fp, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cycles->channel.phys.length);

    NextLine(fp, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cycles->channel.phys.width);

    NextLine(fp, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cycles->channel.phys.bed_elevation);

    NextLine(fp, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cycles->channel.phys.slope);

    NextLine(fp, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cycles->channel.phys.roughness);

    free(elevation);
    free(width);
    fclose(fp);
}

void ReadSoil(const char file_name[], int parameter, cycles_struct *cycles)
{
    int             kz;
    int             kx;
    char            cmdstr[MAXSTRING];
    int             lno = 0;
    double         *value;
    FILE           *fp;

    fp = cycles_fopen(file_name, "r");

    value = (double *)malloc(number_of_columns * sizeof(double));

    for (kz = 0; kz < number_of_layers; kz++)
    {
        NextLine(fp, cmdstr, &lno);
        ReadMultipleValues(cmdstr, number_of_columns, 'd', value);

        for (kx = 0; kx < number_of_columns; kx++)
        {
            switch (parameter)
            {
                case B:
                    cycles->grid[kx].soil.b[kz] = value[kx];
                    break;
                case POROSITY:
                    cycles->grid[kx].soil.porosity[kz] = value[kx];
                    break;
                case AIR_ENTRY_POTENTIAL:
                    cycles->grid[kx].soil.air_entry_potential[kz] = value[kx];
                    break;
                case KSATH:
                    cycles->grid[kx].soil.ksath[kz] = value[kx];
                    break;
                case KSATV:
                    cycles->grid[kx].soil.ksatv[kz] = value[kx];
                    break;
                case INITIAL_CONDITION:
                    cycles->grid[kx].ws.smc[kz] = value[kx];
                    break;
                default:
                    break;
            }
        }
    }

    free(value);
    fclose(fp);
}
