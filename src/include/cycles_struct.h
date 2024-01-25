#ifndef CYCLES_STRUCT_HEADER
#define CYCLES_STRUCT_HEADER

typedef struct control_struct
{
    int             stepsize;               // step size (s)
    int             solver_stepsize;        // solver step size (s)
    double          layer_depth;            // soil layer thickness (m)
} control_struct;

typedef struct soil_struct
{
    double         *porosity;               // soil porosity (m3 m-3)
    double         *air_entry_potential;    // water potential at air entry (m)
    double         *b;                      // coefficient of moisture tension (-)
    double         *ksath;                  // saturated horizontal conductivity (kg s m-3)
    double         *ksatv;                  // saturated vertical conductivity (kg s m-3)
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
    double          uptake[NUMBER_OF_STEPS];
    double          soil_evaporation[NUMBER_OF_STEPS];
    double          infiltration[NUMBER_OF_STEPS];
} forcing_struct;

typedef struct channel_phystate_struct
{
    double          length;                 // channel length (m)
    double          width;                  // channel width (m)
    double          bed_elevation;          // channel bed elevation (m)
    double          slope;                  // channel slope (-)
    double          roughness;              // channel roughness (s m-1/3)
} channel_phystate_struct;

typedef struct channel_wstate_struct
{
    double          stage;                  // channel stage (m)
} channel_wstate_struct;

typedef struct channel_wflux_struct
{
    double          surface;                // surface flux from soil column (m/s)
    double          subsurface;             // subsurface flux from soil column (m/s)
    double          discharge;              // channel discharge (m/s)
} channel_wflux_struct;

typedef struct grid_struct
{
    soil_struct     soil;
    wstate_struct   ws;
    wflux_struct    wf;
    phystate_struct phys;
} grid_struct;

typedef struct channel_struct
{
    channel_phystate_struct phys;
    channel_wstate_struct ws;
    channel_wflux_struct wf;
} channel_struct;

typedef struct output_struct
{
    double         *smc;
    FILE*           smc_fp;
    double         *wp;
    FILE*           wp_fp;
    double         *channel;
    FILE*           channel_fp;
} output_struct;

typedef struct cycles_struct
{
    control_struct  control;
    forcing_struct  forcing;
    channel_struct  channel;
    grid_struct    *grid;
    output_struct   output;
} cycles_struct;

#endif
