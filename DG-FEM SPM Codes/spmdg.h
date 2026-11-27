#ifndef __spmdg_h__
#define __spmdg_h__

#include "spm.h"

/*
 * P1 discontinuous-Galerkin advection on the same two-triangle-per-cell
 * mesh used by FESPM diffusion.
 *
 * Solves h_t + u_H . grad_H(h) = 0 with a spatially varying, temporally
 * frozen nodal P1 velocity field. The equivalent flux formulation
 * h_t + div_H(u_H h) = h div_H(u_H) is used, including the velocity-
 * divergence correction. Time integration is SSP-RK2 and DG-P1 values
 * are transferred to the continuous nodal field by consistent L2
 * projection.
 *
 * At physical inflow boundaries the current implementation applies
 * zero-order extrapolation (exterior trace equals the interior trace).
 * An edge is split at the zero of its linearly varying normal velocity
 * when that velocity changes sign along the edge.
 */
ErrorCode FEDGSPMUpdate_DGAdvection_P1_SSPRK2(FEDGSPM s, double dt,
                                               double h[]);

/*
 * Practical elementwise advective CFL estimate for DG-P1.
 * This is an advective restriction, not a proof of stability for the
 * complete operator when div_H(u_H) is nonzero.
 * Returns DBL_MAX if the horizontal velocity is identically zero.
 */
ErrorCode FEDGSPMGetStableTimestep_DG_P1(FEDGSPM s, double cfl,
                                         double *dt);

#endif
