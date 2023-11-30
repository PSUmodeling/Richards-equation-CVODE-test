#ifndef CYCLES_STRUCT_HEADER
#define CYCLES_STRUCT_HEADER

typedef struct ctrl_struct
{
    double          stepsize;               // step size (s)
} ctrl_struct;

typedef struct soil_struct
{
    double          porosity[NSOIL];        // soil porosity (m3 m-3)
    double          air_entry_pot[NSOIL];   // water potential at air entry (m)
    double          b[NSOIL];               // coefficient of moisture tension (-)
    double          ksat[NSOIL];            // saturated conductivity (kg s m-3)
    double          dsat[NSOIL];            // saturated diffusivity (m2 s-1)
    double          slope;                  // slope for subsurface drainage calculation (-)
} soil_struct;

typedef struct phystate_struct
{
    double          zsoil[NSOIL];           // cumulative depth from land surface to bottom of soil layer (negative) (m)
    double          soil_depth[NSOIL];      // soil layer thickness (m)
} phystate_struct;

typedef struct wstate_struct
{
    double          smc[NSOIL];             // soil water content (m3 m-3)
    double          potential[NSOIL];       // soil water potential (m)
} wstate_struct;

typedef struct wflux_struct
{
    double          uptake[NSOIL];          // soil water uptake (mm)
    double          soil_evap;              // soil evaporation (mm)
    double          infil;                  // soil infiltration (mm)
    double          drainage;               // drainage (mm)
} wflux_struct;

typedef struct forc_struct
{
    double          uptake[NSTEPS];
    double          soil_evap[NSTEPS];
    double          infil[NSTEPS];
} forc_struct;

typedef struct cycles_struct
{
    ctrl_struct     ctrl;
    forc_struct     forc;
    soil_struct     soil;
    wstate_struct   ws;
    wflux_struct    wf;
    phystate_struct phys;
} cycles_struct;

#endif
