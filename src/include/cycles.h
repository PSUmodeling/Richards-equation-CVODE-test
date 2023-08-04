#ifndef CYCLES_HEADER
#define CYCLES_HEADER

#define VERSION     "0.0.0"

#if defined(_MSC_VER)
# define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>

#if defined(_WIN32) || defined(_WIN64)
# include <io.h>
# include <direct.h>
#endif

#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <stdarg.h>
#include <float.h>

#if defined(unix) || defined(__unix__) || defined(__unix)
# include <unistd.h>
# include <fenv.h>
#endif

#if defined(_WIN32) || defined(_WIN64)
# define OS         "Windows"
#endif

#if defined(unix) || defined(__unix__) || defined(__unix)
# define OS         "UNIX"
#endif

#if defined(__APPLE__) || defined(__MACH__)
# define OS         "OS X"
#endif

// SUNDIAL Header Files
#include "cvode/cvode.h"    // Prototypes for CVODE fcts., consts.
#include "sunlinsol/sunlinsol_spgmr.h"  // Access to SPGMR SUNLinearSolver
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"     // Definition of macros SUNSQR and EXP
#include "sundials/sundials_dense.h"    // Prototypes for small dense fcts.

#include "custom_io.h"

#include "cycles_const.h"
#include "cycles_struct.h"
#include "cycles_func.h"

#endif
