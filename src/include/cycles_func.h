#ifndef CYCLES_FUNC_HEADER
#define CYCLES_FUNC_HEADER

#define cycles_exit(...)        _custom_exit(__FILE__, __LINE__, __FUNCTION__, debug_mode,  __VA_ARGS__)
#define cycles_printf(...)      _custom_printf(verbose_mode, __VA_ARGS__)
#define cycles_fopen            _custom_fopen

#define MAX(x, y)               (((x) > (y)) ? (x) : (y))
#define MIN(x, y)               (((x) < (y)) ? (x) : (y))

#define SUBSURFACE(kx, kz)      ((kz) + (kx) * number_of_layers)        // index of state variable for column kx and layer kz
#define CHANNEL                 (number_of_columns * number_of_layers)  // index of channel state variable

void            CheckCVodeFlag(int);
void            Cleanup(cycles_struct *);
void            CloseOutputFiles(output_struct *);
double          Discharge(channel_struct *);
void            GenerateGrids(cycles_struct *);
void            Initialize(cycles_struct *, N_Vector, void **);
double          LateralFlux(double, double, double, double, double, double, double, double, double);
int             Ode(realtype, N_Vector, N_Vector, void *);
void            OpenOutputFiles(output_struct *);
void            ReadDomain(const char [], control_struct *);
void            ReadHydrologicalForcing(forcing_struct *);
void            ReadMultipleValues(const char [], int, char, void *);
void            ReadSoil(const char [], int, cycles_struct *);
void            ReadTopography(const char [], cycles_struct *);
double          RetentionCapacity(double, double, double, double);
void            SetCVodeParameters(void *, SUNLinearSolver *, N_Vector, cycles_struct *);
double          SoilWaterContent(double, double, double, double);
double          SoilWaterPotential(double, double, double, double);
void            SolveCVode(realtype, void *, N_Vector);
void            StoreOutput(const grid_struct *grid, const channel_struct *channel, output_struct *output);
double          SubsurfaceToChannel(int, const channel_struct *, const grid_struct *);
void            SWC(int, cycles_struct *, void *, N_Vector);
double          TotalWaterFlux(grid_struct *, channel_struct *);
double          TotalWaterStorage(grid_struct *, channel_struct *);
void            UpdateStateVariables(N_Vector, cycles_struct *);
double          WaterDiffusivity(double, double, double, double);
double          WaterConductivity(double, double, double, double);
void            WriteOutputFiles(int kstep, int nstep, output_struct *output);
void            ZeroFluxes(N_Vector, cycles_struct *);

#endif
