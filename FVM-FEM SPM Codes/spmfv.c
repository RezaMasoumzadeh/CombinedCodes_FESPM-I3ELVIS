#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "spm.h"

/* ============================================================
 * SECOND-ORDER MUSCL-HANCOCK FINITE-VOLUME ADVECTION
 * ============================================================
 *
 * Governing non-conservative equation:
 *     h_t + u_h . grad(h) = 0
 *
 * Equivalent flux-plus-source form:
 *     h_t + div(u_h h) = h div(u_h)
 *
 * The horizontal velocity may vary in space, but is frozen in time
 * during one advection substep.
 *
 * Discretization:
 *   - vertex-centred median-dual control volumes
 *   - piecewise-linear MUSCL reconstruction
 *   - monotonized-central (MC) slope limiter
 *   - Hancock half-step predictor
 *   - upwind numerical flux
 *
 * The flux divergence and the discrete divergence source use the
 * same face-normal velocities and face lengths. This pairing is
 * important for consistency and preservation of a constant h field.
 *
 * Boundary treatment:
 *   zero-order extrapolation at the physical boundary faces:
 *       h_f = h*_boundary.
 * This is not an externally prescribed inflow boundary condition.
 * ============================================================ */

static double fvm_mc_limiter(double dl, double dr)
{
    double centered;
    double magnitude;

    if (!isfinite(dl) || !isfinite(dr))
        return 0.0;

    if ((dl <= 0.0 && dr >= 0.0) ||
        (dl >= 0.0 && dr <= 0.0))
        return 0.0;

    centered = 0.5 * dl + 0.5 * dr;
    magnitude = 2.0 * fmin(fabs(dl), fabs(dr));
    magnitude = fmin(magnitude, fabs(centered));

    return copysign(magnitude, centered);
}

/* Return 1 for a valid grid/velocity/optional-height state. */
static int fvm_validate_state(FEFVSPM s, const double hc[])
{
    int i;
    int expected_nodes;

    if (s == NULL || s->u_h == NULL)
        return 0;

    if (s->nx < 2 || s->ny < 2)
        return 0;

    if (s->nx > INT_MAX / s->ny)
        return 0;

    expected_nodes = s->nx * s->ny;
    if (s->nnodes != expected_nodes)
        return 0;

    /* u_h is indexed as 2*idx and 2*idx+1. */
    if (s->nnodes > INT_MAX / 2)
        return 0;

    if (!(s->dx > 0.0) || !(s->dy > 0.0) ||
        !isfinite(s->dx) || !isfinite(s->dy))
        return 0;

    for (i = 0; i < 2 * s->nnodes; ++i) {
        if (!isfinite(s->u_h[i]))
            return 0;
    }

    if (hc != NULL) {
        for (i = 0; i < s->nnodes; ++i) {
            if (!isfinite(hc[i]))
                return 0;
        }
    }

    return 1;
}

static double fvm_control_volume(
    int i, int j, int nx, int ny, double dx, double dy)
{
    int boundary_x = (i == 0 || i == nx - 1);
    int boundary_y = (j == 0 || j == ny - 1);

    if (boundary_x && boundary_y)
        return 0.25 * dx * dy;
    if (boundary_x || boundary_y)
        return 0.50 * dx * dy;
    return dx * dy;
}

static void fvm_face_lengths(
    int i, int j, int nx, int ny, double dx, double dy,
    double *len_ew, double *len_ns)
{
    *len_ew = ((j == 0 || j == ny - 1) ? 0.5 * dy : dy);
    *len_ns = ((i == 0 || i == nx - 1) ? 0.5 * dx : dx);
}

ErrorCode FEFVSPMUpdate_FVMAdvection_2ndOrder(
    FEFVSPM s, double dt, double hc[])
{
    int nx, ny, nnodes;
    int i, j, idx;
    double dx, dy;
    double *slope_x = NULL;
    double *slope_y = NULL;
    double *hpred = NULL;
    double *hnew = NULL;

    if (!fvm_validate_state(s, hc)) {
        fprintf(stderr, "FEFVSPM ERROR: Invalid state in FVM advection.\n");
        return 1;
    }
    if (!(dt > 0.0) || !isfinite(dt)) {
        fprintf(stderr, "FEFVSPM ERROR: Invalid timestep in FVM advection.\n");
        return 1;
    }

    nx = s->nx;
    ny = s->ny;
    nnodes = s->nnodes;
    dx = s->dx;
    dy = s->dy;

    if ((size_t)nnodes > ((size_t)-1) / sizeof(double)) {
        fprintf(stderr, "FEFVSPM ERROR: Array size overflow.\n");
        return 1;
    }

    slope_x = (double *)malloc((size_t)nnodes * sizeof(double));
    slope_y = (double *)malloc((size_t)nnodes * sizeof(double));
    hpred = (double *)malloc((size_t)nnodes * sizeof(double));
    hnew = (double *)malloc((size_t)nnodes * sizeof(double));

    if (slope_x == NULL || slope_y == NULL || hpred == NULL || hnew == NULL) {
        fprintf(stderr, "FEFVSPM ERROR: Memory allocation failed.\n");
        free(slope_x);
        free(slope_y);
        free(hpred);
        free(hnew);
        return 1;
    }

    /* STEP 1: MC-limited slopes; one-sided slopes at boundaries. */
    for (j = 0; j < ny; ++j) {
        for (i = 0; i < nx; ++i) {
            double dl, dr;
            idx = i + j * nx;

            if (i > 0 && i < nx - 1) {
                dl = (hc[idx] - hc[idx - 1]) / dx;
                dr = (hc[idx + 1] - hc[idx]) / dx;
                slope_x[idx] = fvm_mc_limiter(dl, dr);
            } else if (i == 0) {
                slope_x[idx] = (hc[idx + 1] - hc[idx]) / dx;
            } else {
                slope_x[idx] = (hc[idx] - hc[idx - 1]) / dx;
            }

            if (j > 0 && j < ny - 1) {
                dl = (hc[idx] - hc[idx - nx]) / dy;
                dr = (hc[idx + nx] - hc[idx]) / dy;
                slope_y[idx] = fvm_mc_limiter(dl, dr);
            } else if (j == 0) {
                slope_y[idx] = (hc[idx + nx] - hc[idx]) / dy;
            } else {
                slope_y[idx] = (hc[idx] - hc[idx - nx]) / dy;
            }
        }
    }

    /* STEP 2: Hancock half-step predictor for h_t + u.grad(h) = 0. */
    for (j = 0; j < ny; ++j) {
        for (i = 0; i < nx; ++i) {
            double ux, uy;
            idx = i + j * nx;
            ux = s->u_h[2 * idx];
            uy = s->u_h[2 * idx + 1];
            hpred[idx] = hc[idx] - 0.5 * dt *
                (ux * slope_x[idx] + uy * slope_y[idx]);

            if (!isfinite(hpred[idx])) {
                fprintf(stderr, "FEFVSPM ERROR: Non-finite Hancock predictor.\n");
                free(slope_x); free(slope_y); free(hpred); free(hnew);
                return 1;
            }
        }
    }

    /* STEP 3: Face fluxes and matching discrete div(u) correction. */
    for (j = 0; j < ny; ++j) {
        for (i = 0; i < nx; ++i) {
            double len_ew, len_ns, vol;
            double flux_sum = 0.0;
            double div_u_sum = 0.0;
            idx = i + j * nx;

            fvm_face_lengths(i, j, nx, ny, dx, dy, &len_ew, &len_ns);
            vol = fvm_control_volume(i, j, nx, ny, dx, dy);

            /* EAST */
            if (i < nx - 1) {
                int eidx = idx + 1;
                double uf = 0.5 * (s->u_h[2 * idx] + s->u_h[2 * eidx]);
                double hL = hpred[idx] + 0.5 * dx * slope_x[idx];
                double hR = hpred[eidx] - 0.5 * dx * slope_x[eidx];
                double hf = (uf > 0.0) ? hL : ((uf < 0.0) ? hR : 0.5 * (hL + hR));
                flux_sum += uf * hf * len_ew;
                div_u_sum += uf * len_ew;
            } else {
                double uf = s->u_h[2 * idx];
                flux_sum += uf * hpred[idx] * len_ew;
                div_u_sum += uf * len_ew;
            }

            /* WEST */
            if (i > 0) {
                int widx = idx - 1;
                double uf = 0.5 * (s->u_h[2 * idx] + s->u_h[2 * widx]);
                double hL = hpred[widx] + 0.5 * dx * slope_x[widx];
                double hR = hpred[idx] - 0.5 * dx * slope_x[idx];
                double hf = (uf > 0.0) ? hL : ((uf < 0.0) ? hR : 0.5 * (hL + hR));
                flux_sum -= uf * hf * len_ew;
                div_u_sum -= uf * len_ew;
            } else {
                double uf = s->u_h[2 * idx];
                flux_sum -= uf * hpred[idx] * len_ew;
                div_u_sum -= uf * len_ew;
            }

            /* NORTH */
            if (j < ny - 1) {
                int nidx = idx + nx;
                double vf = 0.5 * (s->u_h[2 * idx + 1] + s->u_h[2 * nidx + 1]);
                double hL = hpred[idx] + 0.5 * dy * slope_y[idx];
                double hR = hpred[nidx] - 0.5 * dy * slope_y[nidx];
                double hf = (vf > 0.0) ? hL : ((vf < 0.0) ? hR : 0.5 * (hL + hR));
                flux_sum += vf * hf * len_ns;
                div_u_sum += vf * len_ns;
            } else {
                double vf = s->u_h[2 * idx + 1];
                flux_sum += vf * hpred[idx] * len_ns;
                div_u_sum += vf * len_ns;
            }

            /* SOUTH */
            if (j > 0) {
                int sidx = idx - nx;
                double vf = 0.5 * (s->u_h[2 * idx + 1] + s->u_h[2 * sidx + 1]);
                double hL = hpred[sidx] + 0.5 * dy * slope_y[sidx];
                double hR = hpred[idx] - 0.5 * dy * slope_y[idx];
                double hf = (vf > 0.0) ? hL : ((vf < 0.0) ? hR : 0.5 * (hL + hR));
                flux_sum -= vf * hf * len_ns;
                div_u_sum -= vf * len_ns;
            } else {
                double vf = s->u_h[2 * idx + 1];
                flux_sum -= vf * hpred[idx] * len_ns;
                div_u_sum -= vf * len_ns;
            }

            /* STEP 4: update the non-conservative equation. */
            hnew[idx] = hc[idx] - (dt / vol) *
                (flux_sum - hpred[idx] * div_u_sum);

            if (!isfinite(hnew[idx])) {
                fprintf(stderr, "FEFVSPM ERROR: Non-finite updated height.\n");
                free(slope_x); free(slope_y); free(hpred); free(hnew);
                return 1;
            }
        }
    }

    for (idx = 0; idx < nnodes; ++idx)
        hc[idx] = hnew[idx];

    free(slope_x);
    free(slope_y);
    free(hpred);
    free(hnew);
    return 0;
}

/*
 * Conservative local timestep estimate based on the sum of absolute
 * normal speeds over all faces of each median-dual control volume:
 *
 *   dt_i = CFL * V_i / sum_f( |u_n,f| A_f ),
 *   dt   = min_i dt_i.
 *
 * For uniform velocity on a uniform interior cell, CFL=1 corresponds
 * to the familiar coefficient 0.5 in
 *   dt <= 0.5 / (|u|/dx + |v|/dy).
 * This is a timestep estimate, not a proof of positivity or stability
 * for every variable-velocity MUSCL-Hancock problem.
 */
ErrorCode FEFVSPMGetStableTimestep_FVMAdvection(
    FEFVSPM s, double cfl, double *dt)
{
    int nx, ny, i, j, idx;
    double dx, dy;
    double dt_min = DBL_MAX;

    if (dt == NULL) {
        fprintf(stderr, "FEFVSPM ERROR: NULL timestep pointer in CFL routine.\n");
        return 1;
    }
    *dt = 0.0;

    if (!fvm_validate_state(s, NULL)) {
        fprintf(stderr, "FEFVSPM ERROR: Invalid state in CFL routine.\n");
        return 1;
    }
    if (!(cfl > 0.0) || cfl > 1.0 || !isfinite(cfl)) {
        fprintf(stderr, "FEFVSPM ERROR: CFL must satisfy 0 < CFL <= 1.\n");
        return 1;
    }

    nx = s->nx;
    ny = s->ny;
    dx = s->dx;
    dy = s->dy;

    for (j = 0; j < ny; ++j) {
        for (i = 0; i < nx; ++i) {
            double len_ew, len_ns, vol;
            double u_e, u_w, v_n, v_s, rate;
            idx = i + j * nx;

            fvm_face_lengths(i, j, nx, ny, dx, dy, &len_ew, &len_ns);
            vol = fvm_control_volume(i, j, nx, ny, dx, dy);

            if (i < nx - 1) {
                int eidx = idx + 1;
                u_e = 0.5 * (s->u_h[2 * idx] + s->u_h[2 * eidx]);
            } else {
                u_e = s->u_h[2 * idx];
            }

            if (i > 0) {
                int widx = idx - 1;
                u_w = 0.5 * (s->u_h[2 * idx] + s->u_h[2 * widx]);
            } else {
                u_w = s->u_h[2 * idx];
            }

            if (j < ny - 1) {
                int nidx = idx + nx;
                v_n = 0.5 * (s->u_h[2 * idx + 1] + s->u_h[2 * nidx + 1]);
            } else {
                v_n = s->u_h[2 * idx + 1];
            }

            if (j > 0) {
                int sidx = idx - nx;
                v_s = 0.5 * (s->u_h[2 * idx + 1] + s->u_h[2 * sidx + 1]);
            } else {
                v_s = s->u_h[2 * idx + 1];
            }

            rate = (fabs(u_e) + fabs(u_w)) * len_ew +
                   (fabs(v_n) + fabs(v_s)) * len_ns;

            if (!isfinite(rate) || rate < 0.0) {
                fprintf(stderr, "FEFVSPM ERROR: Invalid local CFL rate.\n");
                return 1;
            }

            if (rate > 0.0) {
                double dt_local = cfl * vol / rate;
                if (!(dt_local > 0.0) || !isfinite(dt_local)) {
                    fprintf(stderr, "FEFVSPM ERROR: Invalid local CFL timestep.\n");
                    return 1;
                }
                if (dt_local < dt_min)
                    dt_min = dt_local;
            }
        }
    }

    /* DBL_MAX means that all face-normal velocities are zero. */
    *dt = dt_min;
    return 0;
}
