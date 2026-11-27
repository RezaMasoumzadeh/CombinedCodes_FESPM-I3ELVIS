#ifndef __fespm_h__
#define __fespm_h__

#include "spm.h"

/* Create the FEM diffusion module on the existing FEDGSPM structured grid.
 * Set the desired diffusivity mode and parameters before calling this.
 * Supported modes: DM_CONST (sfcmod=4), DM_NONLINEAR_POWERLAW (5),
 * and DM_NONLINEAR_CRITICALSLOPE (6).
 */
ErrorCode FESPMCreate(FEDGSPM spm, FESPM *fespm);

/* Destroy the FEM diffusion module. */
ErrorCode FESPMDestroy(FESPM *fespm);

/* Perform one FEM diffusion/uplift timestep after DG advection and DG-to-CG L2 projection.
 * PDE: dh/dt = u_v + div(kappa grad h).
 * Uses explicit diffusion: M H^(n+1) = M H* + dt (M U_v - K(h*) H*).
 * For nonlinear modes, K is rebuilt from the supplied h_old at every step.
 */
ErrorCode FESPMDiffusionStep(FESPM fespm, FEDGSPM spm, double h_old[]);

/* Estimate explicit-Euler stability limit from K q = lambda M q. */
ErrorCode FESPMGetStableTimestep(FESPM fespm, double *dt);

#endif
