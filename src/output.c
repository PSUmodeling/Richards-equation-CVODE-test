#include "cycles.h"

void OpenOutputFiles(output_struct *output)
{
    output->smc_fp = fopen("smc.txt", "w");
    output->wp_fp = fopen("wp.txt", "w");
    output->channel_fp = fopen("channel.txt", "w");

    output->smc = (double *)calloc(number_of_columns * number_of_layers, sizeof(double));
    output->wp = (double *)calloc(number_of_columns * number_of_layers, sizeof(double));
    output->channel = (double *)calloc(2, sizeof(double));
}

void CloseOutputFiles(output_struct *output)
{
    fclose(output->smc_fp);
    fclose(output->wp_fp);
    fclose(output->channel_fp);

    free(output->smc);
    free(output->wp);
    free(output->channel);
}

void StoreOutput(const grid_struct *grid, const channel_struct *channel, output_struct *output)
{
    int             kx;
    double          total_width = 0.0;

    for (kx = 0; kx < number_of_columns; kx++)
    {
        int             kz;

        for (kz = 0; kz < number_of_layers; kz++)
        {
            output->smc[kx * number_of_layers + kz] += grid[kx].ws.smc[kz];
            output->wp[kx * number_of_layers + kz] += grid[kx].ws.potential[kz];
        }

        total_width += grid[kx].phys.width;
    }

    output->channel[0] += channel->ws.stage;
    output->channel[1] += channel->wf.discharge * channel->phys.width / total_width * 86400.0 * 1000.0;
}

void WriteOutputFiles(int kstep, int nstep, output_struct *output)
{
    int             k;

    fprintf(output->smc_fp, "%d", kstep);
    fprintf(output->wp_fp, "%d", kstep);
    fprintf(output->channel_fp, "%d", kstep);

    for (k = 0; k < number_of_columns * number_of_layers; k++)
    {
        fprintf(output->smc_fp, ",%.3lf", output->smc[k] / (double)nstep);
        fprintf(output->wp_fp, ",%.3lf", output->wp[k] / (double)nstep);

        output->smc[k] = 0.0;
        output->wp[k] = 0.0;
    }
    fprintf(output->smc_fp, "\n");
    fprintf(output->wp_fp, "\n");

    for (k = 0; k < 2; k++)
    {
        fprintf(output->channel_fp, ",%.3le", output->channel[k] / (double)nstep);
        output->channel[k] = 0.0;
    }
    fprintf(output->channel_fp, "\n");
}
