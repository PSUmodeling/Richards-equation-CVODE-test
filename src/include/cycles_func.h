#ifndef CYCLES_FUNC_HEADER
#define CYCLES_FUNC_HEADER

#define cycles_exit(...)        _custom_exit(__FILE__, __LINE__, __FUNCTION__, debug_mode,  __VA_ARGS__)
#define cycles_printf(...)      _custom_printf(verbose_mode, __VA_ARGS__)
#define cycles_fopen            _custom_fopen

#define MAX(x, y)               (((x) > (y)) ? (x) : (y))
#define MIN(x, y)               (((x) < (y)) ? (x) : (y))

void            CheckCVodeFlag(int);
void            Initialize(cycles_struct *, N_Vector, void **);
int             Ode(realtype, N_Vector, N_Vector, void *);
void            ReadHydro(forc_struct *);
double          RetentionCapacity(double, double, double, double);
void            SetCVodeParam(void *, SUNLinearSolver *, N_Vector, cycles_struct *);
double          SoilWaterContent(double, double, double, double);
double          SoilWaterPot(double, double, double, double);
void            SolveCVode(realtype, void *, N_Vector);
void            SWC(int, cycles_struct *, void *, N_Vector);
void            WDfCnd(double, double, double, double, double, double *, double *);

double CapillaryDrive(double bexp, double air_entry_potential, double swc_b, double swc);
double DryDepth(double dt, const soil_struct *soil, const phystate_struct *phys, double swc_b);
void WettingFrontPosition(const phystate_struct *phys,  int *layer, double wetting_front_thickness[]);
double SoilMoistureAtWettingFront(int layer, double wetting_front_thickness, double soil_depth, double swc_b, double swc, double smcmax);
double InfiltrationCapacity(const soil_struct *soil, const phystate_struct *phys, int n, const double wetting_front_thickness[]);
double CompositeConductivity(int n, double wetting_front_depth, const double wetting_front_thickness[], const double swc[], const soil_struct *soil);
double WettingFrontAdvanceRate(int n, double wetting_front_depth, const double wetting_front_thickness[], double swc_b, const double swc_wf[], const soil_struct *soil);

#endif
