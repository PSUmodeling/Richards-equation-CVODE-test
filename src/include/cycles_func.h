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
void            SetCVodeParam(void *, SUNLinearSolver *, N_Vector, cycles_struct *);
void            SolveCVode(realtype, void *, N_Vector);
void            SWC(int, cycles_struct *, void *, N_Vector);
void            WDfCnd(double, double, double, double, double, double *, double *);

#endif
