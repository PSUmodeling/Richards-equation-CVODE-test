#include "cycles.h"

void OpenOutputFiles(output_struct *output)
{
    output->smc = fopen("smc.txt", "w");
    output->wp = fopen("wp.txt", "w");
    output->channel = fopen("channel.txt", "w");
}

void CloseOutputFiles(output_struct *output)
{
    fclose(output->smc);
    fclose(output->wp);
    fclose(output->channel);
}

void WriteOutputFiles(int kstep, const grid_struct *grid, const channel_struct *channel, output_struct *output)
{
    int             kx;
    double          total_width = 0.0;

    fprintf(output->smc, "%d", kstep);
    fprintf(output->wp, "%d", kstep);
    fprintf(output->channel, "%d", kstep);

    for (kx = 0; kx < number_of_columns; kx++)
    {
        int             kz;

        for (kz = 0; kz < number_of_layers; kz++)
        {
            fprintf(output->smc, ",%.3lf", grid[kx].ws.smc[kz]);
            fprintf(output->wp, ",%.3lf", grid[kx].ws.potential[kz]);
        }

        total_width += grid[kx].phys.width;
    }
    fprintf(output->smc, "\n");
    fprintf(output->wp, "\n");

    fprintf(output->channel, ",%.3le,%.3lf\n", channel->ws.stage, channel->wf.discharge * channel->phys.width / total_width * 86400.0 * 1000.0);
}