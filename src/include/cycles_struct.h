#ifndef CYCLES_STRUCT_HEADER
#define CYCLES_STRUCT_HEADER

typedef struct control_struct
{
    double          stepsize;               // step size (s)
    double          layer_depth;            // soil layer thickness (m)
} control_struct;

typedef struct soil_struct
{
    double         *porosity;               // soil porosity (m3 m-3)
    double         *air_entry_potential;    // water potential at air entry (m)
    double         *b;                      // coefficient of moisture tension (-)
    double         *ksath;                  // saturated horizontal conductivity (kg s m-3)
    double         *ksatv;                  // saturated vertical conductivity (kg s m-3)
    double         *dsat;                   // saturated diffusivity (m2 s-1)
} soil_struct;

typedef struct phystate_struct
{
    double         *zsoil;                  // cumulative depth from land surface to bottom of soil layer (negative) (m)
    double         *soil_depth;             // soil layer thickness (m)
    double          width;                  // soil column width (m)
    double         *retention_capacity;     // soil water retention capacity (m-1)
} phystate_struct;

typedef struct wstate_struct
{
    double         *smc;                    // soil water content (m3 m-3)
    double         *potential;              // soil water potential (m)
} wstate_struct;

typedef struct wflux_struct
{
    double         *uptake;                 // soil water uptake (mm)
    double         *lateral;                // lateral flow (mm)
    double          soil_evaporation;       // soil evaporation (mm)
    double          infiltration;           // soil infiltration (mm)
} wflux_struct;

typedef struct forcing_struct
{
    double          uptake[NSTEPS];
    double          soil_evaporation[NSTEPS];
    double          infiltration[NSTEPS];
} forcing_struct;

typedef struct river_struct
{
    double          width;                  // river width (m)
    double          depth;                  // river depth (m)
} river_struct;

typedef struct grid_struct
{
    soil_struct     soil;
    wstate_struct   ws;
    wflux_struct    wf;
    phystate_struct phys;
} grid_struct;

typedef struct cycles_struct
{
    control_struct  control;
    forcing_struct  forcing;
    river_struct    river;
    grid_struct    *grid;
} cycles_struct;

#endif
