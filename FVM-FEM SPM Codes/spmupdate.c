#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spm.h"
#include "fespm.h"


/* ================================================================
 * Timestep
 * ================================================================ */

ErrorCode FEFVSPMComputeTimestep_const(FEFVSPM s, double *dt)
{
  *dt = s->dt_const;
  return (0);
}

ErrorCode FEFVSPMComputeTimestep(FEFVSPM s, double h[], double *dt)
{
  ErrorCode ierr;

  if (!s || !dt) return 1;

  switch (s->dtmode)
  {
  case DTM_CONST:
    ierr = FEFVSPMComputeTimestep_const(s, dt);
    if (ierr != 0) return ierr;
    break;

  case DTM_EXPLICIT:
  {
    double dt_diffusion, dt_advection;

    if (!s->fespm) {
      printf("FEFVSPM ERROR: FEM module is required for explicit timestep control.\\n");
      return 1;
    }

    ierr = FESPMGetStableTimestep(s->fespm, &dt_diffusion);
    if (ierr != 0) return ierr;

    ierr = FEFVSPMGetStableTimestep_FVMAdvection(s, FEFVSPM_FVM_ADVECTION_CFL, &dt_advection);
    if (ierr != 0) return ierr;

    *dt = (dt_diffusion < dt_advection) ? dt_diffusion : dt_advection;
    *dt *= s->dt_theta;
    break;
  }

  default:
    printf("FEFVSPM ERROR: unsupported timestep mode.\\n");
    return 1;
  }

  if (!(*dt > 0.0) || isnan(*dt) || isinf(*dt)) {
    printf("FEFVSPM ERROR: invalid timestep.\\n");
    return 1;
  }
  /* Do not raise an explicitly stability-limited step to dt_min. */
  if (s->dtmode == DTM_CONST && *dt < s->dt_min) *dt = s->dt_min;
  if (*dt > s->dt_max) *dt = s->dt_max;
  if (*dt > s->time_end) *dt = s->time_end;

  return 0;
}

ErrorCode FEFVSPMView(FEFVSPM s)
{
  printf("FEFVSPM:\n");
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
    printf("    dt(%d) %1.4e\n", s->nt, s->dt);
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

ErrorCode FEFVSPMStepInitialize(FEFVSPM s)
{

  s->time = 0.0;
  s->nt = 0;

  memcpy(s->h, s->h_ic, sizeof(double) * s->nnodes);
  memcpy(s->h_c, s->h_ic, sizeof(double) * s->nnodes);

  return FEFVSPMCheckpoint(s);
}

ErrorCode FEFVSPMStepAdvDiff(FEFVSPM s)
{
  ErrorCode ierr;
  int tk;

  printf("____ FEFVSPM AdvDiff integrator ____\n");
  printf("  current time/step %1.4e/%d : target time %1.4e : max steps %d\n", s->time, s->nt, s->time_end, s->nt_max);

  if (s->time >= s->time_end)
  {
    printf("  FEFVSPM: No updates performed as already at target end time\n");
    return (0);
  }

  for (tk = 0; tk < s->nt_max; tk++)
  {
    ierr = FEFVSPMComputeTimestep(s, s->h, &s->dt);
    if (ierr != 0) {
      printf("FATAL: FEFVSPMComputeTimestep failed\n");
      return ierr;
    }

    /* Enforce MUSCL FVM advection CFL */
    {
      double dt_adv_stable;
      ierr = FEFVSPMGetStableTimestep_FVMAdvection(s, FEFVSPM_FVM_ADVECTION_CFL, &dt_adv_stable);
      if (ierr != 0) {
        printf("FATAL: FEFVSPMGetStableTimestep_FVMAdvection failed\n");
        return ierr;
      }
      if (dt_adv_stable > 0.0 && s->dt > dt_adv_stable) {
        s->dt = dt_adv_stable;
      }
    }

    if (s->time + s->dt > s->time_end) s->dt = s->time_end - s->time;

    /* Second-order MUSCL FVM advection */
    ierr = FEFVSPMUpdate_FVMAdvection_2ndOrder(s, s->dt, s->h);
    if (ierr != 0) {
      printf("FATAL: FEFVSPMUpdate_FVMAdvection_2ndOrder failed\n");
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
      ierr = FEFVSPMFieldView(s, FVT_VTRB);
      ierr = FEFVSPMCheckpoint(s);
    }
    if (s->time >= s->time_end)
    {
      printf("  FEFVSPM: Reached target end time\n");
      break;
    }
  }
  if (tk == (s->nt_max)) printf("  FEFVSPM: Reached max number of time steps permitted\n");
  return (0);
}

