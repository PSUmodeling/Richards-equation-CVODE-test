#ifndef CYCLES_CONST_HEADER
#define CYCLES_CONST_HEADER

#define NSOIL                   4           // number of soil layers
#define NSTEPS                  1440        // number of time steps

#define RELATIVE_ERROR          1.0E-3      // relative error tolerance
#define ABSOLUTE_ERROR          1.0E-4      // absolute error tolerance (m3/m3)
#define STEPSIZE                1800.0      // step size (s)
#define INITIAL_STEP            300.0       // initial step size (s)

extern int      verbose_mode;
extern int      debug_mode;

#endif
