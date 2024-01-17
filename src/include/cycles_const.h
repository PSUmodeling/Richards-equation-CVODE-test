#ifndef CYCLES_CONST_HEADER
#define CYCLES_CONST_HEADER

#define NUMBER_OF_STEPS         1440        // number of time steps

#define RELATIVE_ERROR          1.0E-3      // relative error tolerance
#define ABSOLUTE_ERROR          1.0E-4      // absolute error tolerance (m3/m3)
#define STEPSIZE                1800.0      // step size (s)
#define INITIAL_STEP            300.0       // initial step size (s)

#define B                       0
#define KSATH                   1
#define KSATV                   2
#define POROSITY                3
#define AIR_ENTRY_POTENTIAL     4
#define INITIAL_CONDITION       5

extern int      number_of_columns;      // number of horizontal columns
extern int      number_of_layers;       // number of soil layers
extern int      verbose_mode;
extern int      debug_mode;

#endif
