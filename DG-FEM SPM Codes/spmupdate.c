#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spm.h"
#include "fespm.h"
#include "spmdg.h"


/* ================================================================
 * Timestep
 * ================================================================ */


ErrorCode FEDGSPMComputeTimestep(FEDGSPM s, double h[], double *dt)
{
  ErrorCode ierr;

  if (s == NULL || dt == NULL || s->fespm == NULL)
    return (1);

  switch (s->dtmode)
  {
  case DTM_CONST:
    *dt = s->dt_const;
    break;

  case DTM_EXPLICIT:
  {
    double dt_diffusion, dt_advection;

    /* FEM explicit-Euler stability limit from the generalized eigenvalue
       problem K q = lambda M q. */
    ierr = FESPMGetStableTimestep(s->fespm, &dt_diffusion);
    if (ierr != 0) return ierr;

    /* DG-P1 advection CFL limit. */
    ierr = FEDGSPMGetStableTimestep_DG_P1(s, 0.20, &dt_advection);
    if (ierr != 0) return ierr;

    *dt = (dt_diffusion < dt_advection) ? dt_diffusion : dt_advection;
    break;
  }

  default:
    printf("*** FEDGSPM ERROR: Unknown timestep mode %d\n", (int)s->dtmode);
    return (1);
  }

  if (!isfinite(*dt) || *dt <= 0.0)
  {
    printf("*** FEDGSPM ERROR: Invalid timestep: %.16e\n", *dt);
    return (1);
  }

  /* Never enlarge a stability-limited timestep to satisfy dt_min. */
  if (s->dt_max > 0.0 && *dt > s->dt_max)
    *dt = s->dt_max;

  if (s->dt_min > 0.0 && *dt < s->dt_min)
    printf("*** FEDGSPM WARNING: stable timestep %.6e is below dt_min %.6e;\n"
           "    retaining the stable timestep.\n", *dt, s->dt_min);

  if (s->time_end > s->time && *dt > s->time_end - s->time)
    *dt = s->time_end - s->time;

  return (0);
}








ErrorCode FEDGSPMStepInitialize(FEDGSPM s)
{
  if (s == NULL || s->h == NULL || s->h_ic == NULL || s->h_c == NULL)
    return (1);

  s->time = 0.0;
  s->nt = 0;

  /* FEM homogeneous Neumann conditions are natural boundary conditions;
     no finite-difference ghost-node boundary copying is required. */
  memcpy(s->h, s->h_ic, sizeof(double) * s->nnodes);
  memcpy(s->h_c, s->h_ic, sizeof(double) * s->nnodes);

  return (0);
}

ErrorCode FEDGSPMView(FEDGSPM s)
{
  printf("FEDGSPM:\n");
  printf("  domain [%1.4e,%1.4e] x [%1.4e,%1.4e] \n", s->xrange[0], s->xrange[1], s->yrange[0], s->yrange[1]);
  printf("  nx %4d\n", s->nx);
  printf("  ny %4d\n", s->ny);
  printf("  time_end %1.4e\n", s->time_end);
  printf("  maxsteps %d\n", s->nt_max);
  printf("  output every %d\n", s->nt_io);

  printf("  timestep mode %d\n", (int)s->dtmode);
  switch (s->dtmode)
  {
  case DTM_CONST:
    printf("    dt_const %1.4e\n", s->dt_const);
    break;
  case DTM_EXPLICIT:
    printf("    dt_theta %1.4e\n", s->dt_theta);
    printf("    dt(%D) %1.4e\n", s->nt, s->dt);
    break;
  }

  printf("  diffusion mode %d\n", (int)s->mode);
  switch (s->mode)
  {
  case DM_CONST:
    printf("    k0 %1.4e\n", s->k0);
    break;
  case DM_NONLINEAR_CRITICALSLOPE:
    printf("    k0 %1.4e , m %1.4e Sc %1.4e\n", s->k0, s->d_exp, s->Sc);
    break;
  case DM_NONLINEAR_POWERLAW:
    printf("    k0 %1.4e , m %1.4e\n", s->k0, s->d_exp);
    break;
  case DM_UNDEFINED:
    break;
  }
  return (0);
}




/* ================================================================
 * AdvDiff integrator: DG-P1 advection followed by FEM diffusion
 * ================================================================ */
ErrorCode FEDGSPMStepAdvDiff(FEDGSPM s)
{
  ErrorCode ierr;
  int tk;

  printf("____ DG-P1 advection + FEM-P1 diffusion integrator ____\n");
  printf("  current time/step %1.4e/%d : target time %1.4e : max steps %d\n", s->time, s->nt, s->time_end, s->nt_max);

  if (s->time >= s->time_end)
  {
    printf("  FEDGSPM: No updates performed as already at target end time\n");
    return (0);
  }

  for (tk = 0; tk < s->nt_max; tk++)
  {
    ierr = FEDGSPMComputeTimestep(s, s->h, &s->dt);
    if (ierr != 0) {
      printf("FATAL: FEDGSPMComputeTimestep failed\n");
      return ierr;
    }

    if (s->time + s->dt > s->time_end) s->dt = s->time_end - s->time;

    /* DG-P1 advection with SSP-RK2. */
    ierr = FEDGSPMUpdate_DGAdvection_P1_SSPRK2(s, s->dt, s->h);
    if (ierr != 0) {
      printf("FATAL: FEDGSPMUpdate_DGAdvection_P1_SSPRK2 failed\n");
      return ierr;
    }

    memcpy(s->h_c, s->h, sizeof(double) * s->nnodes);

    ierr = FESPMDiffusionStep(s->fespm, s, s->h_c);
    if (ierr != 0) {
      printf("FATAL: FESPMDiffusionStep failed\n");
      return ierr;
    }

    s->time += s->dt;
    s->nt++;
    printf("  step %4d: time %1.4e, dt %1.4e\n", s->nt, s->time, s->dt);

    if (s->nt % s->nt_io == 0)
    {
      ierr = FEDGSPMFieldView(s, FVT_VTRB);
      ierr = FEDGSPMCheckpoint(s);
    }
    if (s->time >= s->time_end)
    {
      printf("  FEDGSPM: Reached target end time\n");
      break;
    }
  }
  if (tk == (s->nt_max)) printf("  FEDGSPM: Reached max number of time steps permitted\n");
  return (0);
}

ErrorCode FEDGSPMStepAdvDiff_TimeDependentInput(FEDGSPM s, ErrorCode (*user_update_function)(FEDGSPM, double, int, double, void *), void *ctx)
{
  ErrorCode ierr;
  int tk;

  printf("____ DG-P1 advection + FEM-P1 diffusion (time-dependent input) integrator ____\n");
  printf("  current time/step %1.4e/%d : target time %1.4e : max steps %d\n", s->time, s->nt, s->time_end, s->nt_max);

  if (s->time >= s->time_end)
  {
    printf("  FEDGSPM: No updates performed as already at target end time\n");
    return (0);
  }

  for (tk = 0; tk < s->nt_max; tk++)
  {
    if (user_update_function) ierr = user_update_function(s, s->time, s->nt, s->dt, ctx);

    ierr = FEDGSPMComputeTimestep(s, s->h, &s->dt);
    if (ierr != 0) return ierr;

    if (s->time + s->dt > s->time_end) s->dt = s->time_end - s->time;

    ierr = FEDGSPMUpdate_DGAdvection_P1_SSPRK2(s, s->dt, s->h);
    if (ierr != 0) return ierr;

    memcpy(s->h_c, s->h, sizeof(double) * s->nnodes);
    ierr = FESPMDiffusionStep(s->fespm, s, s->h_c);
    if (ierr != 0) return ierr;

    s->time += s->dt;
    s->nt++;
    printf("  step %4d: time %1.4e, dt %1.4e\n", s->nt, s->time, s->dt);

    if (s->nt % s->nt_io == 0)
    {
      ierr = FEDGSPMFieldView(s, FVT_VTRB);
      ierr = FEDGSPMCheckpoint(s);
    }
    if (s->time >= s->time_end)
    {
      printf("  FEDGSPM: Reached target end time\n");
      break;
    }
  }
  if (tk == (s->nt_max)) printf("  FEDGSPM: Reached max number of time steps permitted\n");
  return (0);
}


