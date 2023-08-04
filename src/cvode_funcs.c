#include "cycles.h"

void SetCVodeParam(void *cvode_mem, SUNLinearSolver *sun_ls, N_Vector CV_Y, cycles_struct *cycles)
{
    int             cv_flag;
    static int      reset;

    if (reset)
    {
        // When model spins-up and recycles forcing, use CVodeReInit to reset solver time, which does not allocates
        // memory
        cv_flag = CVodeReInit(cvode_mem, 0.0, CV_Y);
        CheckCVodeFlag(cv_flag);
    }
    else
    {
        cv_flag = CVodeInit(cvode_mem, Ode, 0.0, CV_Y);
        CheckCVodeFlag(cv_flag);
        reset = 1;

        *sun_ls = SUNLinSol_SPGMR(CV_Y, PREC_NONE, 0);

        // Attach the linear solver
        CVodeSetLinearSolver(cvode_mem, *sun_ls, NULL);

        // Set CVode tolerances
        // Note that if gas transport is added, the absolute tolerance needs to be changed
        cv_flag = CVodeSStolerances(cvode_mem, RELATIVE_ERROR, ABSOLUTE_ERROR);
        CheckCVodeFlag(cv_flag);

        // Specifies Cycles data block and attaches it to the main cvode memory block
        cv_flag = CVodeSetUserData(cvode_mem, cycles);
        CheckCVodeFlag(cv_flag);

        // Specifies the initial step size
        cv_flag = CVodeSetInitStep(cvode_mem, INITIAL_STEP);
        CheckCVodeFlag(cv_flag);

        // Indicates if the BDF stability limit detection algorithm should be used
        cv_flag = CVodeSetStabLimDet(cvode_mem, SUNTRUE);
        CheckCVodeFlag(cv_flag);

        // Specifies an upper bound on the magnitude of the step size
        cv_flag = CVodeSetMaxStep(cvode_mem, (realtype)cycles->ctrl.stepsize);
        CheckCVodeFlag(cv_flag);

#if TEMP_DISABLED
        // Specifies the maximum number of steps to be taken by the solver in its attempt to reach the next output time
        cv_flag = CVodeSetMaxNumSteps(cvode_mem, cycles->ctrl.stepsize * 10);
        CheckCVodeFlag(cv_flag);
#endif
    }
}

void CheckCVodeFlag(int cv_flag)
{
    if (cv_flag < 0)
    {
        cycles_printf(VL_ERROR, "CVODE error %s\n", CVodeGetReturnFlagName(cv_flag));
        cycles_exit(EXIT_FAILURE);
    }
}

void SolveCVode(realtype tout, void *cvode_mem, N_Vector CV_Y)
{
    int             cv_flag;
    realtype        solvert;

    // Specifies the value of the independent variable t past which the solution is not to proceed
    cv_flag = CVodeSetStopTime(cvode_mem, tout);
    CheckCVodeFlag(cv_flag);

    cv_flag = CVode(cvode_mem, tout, CV_Y, &solvert, CV_NORMAL);
    CheckCVodeFlag(cv_flag);
}
