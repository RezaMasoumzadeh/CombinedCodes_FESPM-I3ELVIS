/*
 * P1 discontinuous-Galerkin advection on the same two-triangle-per-cell
 * mesh used by FESPM diffusion.
 *
 * Governing equation:
 *
 *     h_t + u . grad(h) = 0
 *
 * The velocity is frozen in time during one surface substep but may vary
 * spatially. No divergence-free assumption is made. The equivalent
 * flux-form equation is
 *
 *     h_t + div(u h) = h div(u).
 *
 * Therefore the DG residual includes the elementwise velocity-divergence
 * contribution. For a nodal P1 velocity field, div(u) is constant in each
 * triangle and M^{-1} S = (div u) U exactly.
 *
 * The boundary inflow state currently uses zero-order extrapolation
 * (exterior trace = interior trace). This is a zero-gradient inflow
 * treatment; applications requiring prescribed inflow data should extend
 * the interface to supply those data explicitly.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "spmdg.h"

typedef struct {
    int a, b, tri, face;
} DGEdge;

typedef struct {
    int node[3];
    int neighbor[3];
    int neighbor_face[3];
    double area;
    double x[3], y[3];
    double gx[3], gy[3];
} DGTriangle;

static int dg_edge_cmp(const void *pa, const void *pb)
{
    const DGEdge *a = (const DGEdge *)pa;
    const DGEdge *b = (const DGEdge *)pb;
    if (a->a != b->a) return (a->a < b->a) ? -1 : 1;
    if (a->b != b->b) return (a->b < b->b) ? -1 : 1;
    return a->tri - b->tri;
}

static int dg_node(int i, int j, int nx)
{
    return i + j * nx;
}

static int dg_make_tri(DGTriangle *t, FEDGSPM s, int a, int b, int c)
{
    int q;
    double det;

    t->node[0] = a;
    t->node[1] = b;
    t->node[2] = c;

    for (q = 0; q < 3; ++q) {
        int n = t->node[q];
        int i = n % s->nx;
        int j = n / s->nx;

        t->x[q] = s->xrange[0] + i * s->dx;
        t->y[q] = s->yrange[0] + j * s->dy;
        t->neighbor[q] = -1;
        t->neighbor_face[q] = -1;
    }

    det = (t->x[1] - t->x[0]) * (t->y[2] - t->y[0])
        - (t->x[2] - t->x[0]) * (t->y[1] - t->y[0]);

    if (!isfinite(det) || fabs(det) <= DBL_MIN)
        return 1;

    /* Ensure counter-clockwise orientation so the face normal below
       is consistently outward. */
    if (det < 0.0) {
        int tmp = t->node[1];
        double tx, ty;

        t->node[1] = t->node[2];
        t->node[2] = tmp;

        tx = t->x[1]; t->x[1] = t->x[2]; t->x[2] = tx;
        ty = t->y[1]; t->y[1] = t->y[2]; t->y[2] = ty;
        det = -det;
    }

    t->area = 0.5 * det;

    t->gx[0] = (t->y[1] - t->y[2]) / det;
    t->gx[1] = (t->y[2] - t->y[0]) / det;
    t->gx[2] = (t->y[0] - t->y[1]) / det;

    t->gy[0] = (t->x[2] - t->x[1]) / det;
    t->gy[1] = (t->x[0] - t->x[2]) / det;
    t->gy[2] = (t->x[1] - t->x[0]) / det;

    return 0;
}

static int dg_build_mesh(FEDGSPM s, DGTriangle **tout, int *ntout)
{
    int i, j, k = 0;
    int nt, ne, p = 0;
    DGTriangle *tri = NULL;
    DGEdge *edges = NULL;

    if (!s || !tout || !ntout || s->nx < 2 || s->ny < 2 ||
        s->mx != s->nx - 1 || s->my != s->ny - 1 ||
        !s->u_h || !isfinite(s->dx) || !isfinite(s->dy) ||
        s->dx <= 0.0 || s->dy <= 0.0)
        return 1;

    if (s->my <= 0 || s->mx > (INT_MAX / 2) / s->my)
        return 1;

    nt = 2 * s->mx * s->my;
    ne = 3 * nt;

    tri = (DGTriangle *)calloc((size_t)nt, sizeof(*tri));
    edges = (DGEdge *)malloc((size_t)ne * sizeof(*edges));
    if (!tri || !edges) {
        free(tri);
        free(edges);
        return 1;
    }

    for (j = 0; j < s->my; ++j) {
        for (i = 0; i < s->mx; ++i) {
            int n0 = dg_node(i,     j,     s->nx);
            int n1 = dg_node(i + 1, j,     s->nx);
            int n2 = dg_node(i + 1, j + 1, s->nx);
            int n3 = dg_node(i,     j + 1, s->nx);

            if (dg_make_tri(&tri[k++], s, n0, n1, n2) ||
                dg_make_tri(&tri[k++], s, n0, n2, n3)) {
                free(edges);
                free(tri);
                return 1;
            }
        }
    }

    for (k = 0; k < nt; ++k) {
        for (i = 0; i < 3; ++i) {
            int a = tri[k].node[i];
            int b = tri[k].node[(i + 1) % 3];

            edges[p].a = (a < b) ? a : b;
            edges[p].b = (a < b) ? b : a;
            edges[p].tri = k;
            edges[p].face = i;
            ++p;
        }
    }

    qsort(edges, (size_t)ne, sizeof(*edges), dg_edge_cmp);

    for (p = 0; p < ne;) {
        int q = p + 1;

        while (q < ne &&
               edges[q].a == edges[p].a &&
               edges[q].b == edges[p].b)
            ++q;

        if (q - p == 2) {
            DGEdge e0 = edges[p];
            DGEdge e1 = edges[p + 1];

            tri[e0.tri].neighbor[e0.face] = e1.tri;
            tri[e0.tri].neighbor_face[e0.face] = e1.face;
            tri[e1.tri].neighbor[e1.face] = e0.tri;
            tri[e1.tri].neighbor_face[e1.face] = e0.face;
        } else if (q - p > 2) {
            fprintf(stderr, "DG ERROR: non-manifold edge detected\n");
            free(edges);
            free(tri);
            return 1;
        }

        p = q;
    }

    free(edges);
    *tout = tri;
    *ntout = nt;
    return 0;
}

static double dg_interp_edge(const double u[3], int f, double t)
{
    return (1.0 - t) * u[f] + t * u[(f + 1) % 3];
}

static void dg_velocity_at(const DGTriangle *t, const double *uh,
                           const double lam[3], double *ux, double *uy)
{
    int a;

    *ux = 0.0;
    *uy = 0.0;

    for (a = 0; a < 3; ++a) {
        int n = t->node[a];
        *ux += lam[a] * uh[2 * n];
        *uy += lam[a] * uh[2 * n + 1];
    }
}

/* Elementwise divergence of the nodal P1 velocity interpolant. */
static double dg_velocity_divergence(const DGTriangle *t, const double *uh)
{
    int a;
    double divu = 0.0;

    for (a = 0; a < 3; ++a) {
        int n = t->node[a];
        divu += uh[2 * n]     * t->gx[a]
              + uh[2 * n + 1] * t->gy[a];
    }

    return divu;
}

/*
 * Compute the DG spatial residual R such that
 *
 *     M_K dU_K/dt = V_K - F_K + S_K.
 *
 * The volume and divergence terms use a degree-two triangle rule.
 * The edge flux uses the two-point Gauss rule, exact through degree 3
 * on each segment. For nodal P1 velocity, u.n is linear along an edge;
 * an edge is split at its exact zero whenever u.n changes sign.
 */
static void dg_rhs(FEDGSPM s, const DGTriangle *tri, int nt,
                   const double *U, double *rhs)
{
    static const double lq[3][3] = {
        {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
        {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
        {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}
    };
    const double tq[2] = {
        0.5 - 0.5 / sqrt(3.0),
        0.5 + 0.5 / sqrt(3.0)
    };

    int e, a, q, f;

    memset(rhs, 0, (size_t)3 * (size_t)nt * sizeof(double));

    for (e = 0; e < nt; ++e) {
        const DGTriangle *t = &tri[e];
        const double *ue = &U[3 * e];
        double r[3] = {0.0, 0.0, 0.0};
        double divu = dg_velocity_divergence(t, s->u_h);

        /* V_i = integral_K h (u . grad(phi_i)) dK. */
        for (q = 0; q < 3; ++q) {
            double ux, uy, hq;

            dg_velocity_at(t, s->u_h, lq[q], &ux, &uy);
            hq = lq[q][0] * ue[0]
               + lq[q][1] * ue[1]
               + lq[q][2] * ue[2];

            for (a = 0; a < 3; ++a) {
                r[a] += (t->area / 3.0) * hq
                      * (ux * t->gx[a] + uy * t->gy[a]);
            }
        }

        /* F_i = integral_boundary h* (u . n) phi_i ds. */
        for (f = 0; f < 3; ++f) {
            int v0 = f;
            int v1 = (f + 1) % 3;
            int nb = t->neighbor[f];
            int nbf = t->neighbor_face[f];
            double dx = t->x[v1] - t->x[v0];
            double dy = t->y[v1] - t->y[v0];
            double len = sqrt(dx * dx + dy * dy);
            double nx = dy / len;
            double ny = -dx / len;

            {
                /*
                 * For a nodal P1 velocity, u.n is linear along the edge.
                 * Split the edge at its exact zero, if the normal velocity
                 * changes sign, so each integration segment has one
                 * well-defined upwind direction.
                 */
                int nseg = 1, seg;
                double zbreak[3] = {0.0, 1.0, 1.0};
                int node0 = t->node[v0];
                int node1 = t->node[v1];
                double un0 = s->u_h[2 * node0] * nx
                           + s->u_h[2 * node0 + 1] * ny;
                double un1 = s->u_h[2 * node1] * nx
                           + s->u_h[2 * node1 + 1] * ny;

                if (un0 * un1 < 0.0) {
                    double zzero = un0 / (un0 - un1);
                    if (zzero > 0.0 && zzero < 1.0) {
                        zbreak[1] = zzero;
                        zbreak[2] = 1.0;
                        nseg = 2;
                    }
                }

                for (seg = 0; seg < nseg; ++seg) {
                    double za = zbreak[seg];
                    double zb = zbreak[seg + 1];

                    for (q = 0; q < 2; ++q) {
                        double z = za + (zb - za) * tq[q];
                        double lam[3] = {0.0, 0.0, 0.0};
                        double ux, uy, un, hminus, hplus, hhat;
                        double weight = 0.5 * len * (zb - za);

                        lam[v0] = 1.0 - z;
                        lam[v1] = z;

                        dg_velocity_at(t, s->u_h, lam, &ux, &uy);
                        un = ux * nx + uy * ny;

                        hminus = dg_interp_edge(ue, f, z);

                        if (nb >= 0) {
                            const double *unb = &U[3 * nb];
                            int same_orientation =
                                (tri[nb].node[nbf] == t->node[v0]);
                            double zn = same_orientation ? z : 1.0 - z;

                            hplus = dg_interp_edge(unb, nbf, zn);
                        } else {
                            /*
                             * Physical boundary: zero-order extrapolation
                             * at inflow (exterior trace = interior trace).
                             */
                            hplus = hminus;
                        }

                        hhat = (un >= 0.0) ? hminus : hplus;

                        for (a = 0; a < 3; ++a)
                            r[a] -= weight * lam[a] * un * hhat;
                    }
                }
            }
        }

        /*
         * Apply M_K^{-1} to V-F. For
         * M_K = A/12 * [[2,1,1],[1,2,1],[1,1,2]],
         * M_K^{-1} = 1/A * [[9,-3,-3],[-3,9,-3],[-3,-3,9]].
         */
        rhs[3 * e + 0] =
            (3.0 / t->area) * (3.0 * r[0] - r[1] - r[2]);
        rhs[3 * e + 1] =
            (3.0 / t->area) * (-r[0] + 3.0 * r[1] - r[2]);
        rhs[3 * e + 2] =
            (3.0 / t->area) * (-r[0] - r[1] + 3.0 * r[2]);

        /*
         * S_K = div(u)_K M_K U_K. Hence M_K^{-1} S_K
         * = div(u)_K U_K exactly. This is the correction required
         * by the original non-conservative equation.
         */
        rhs[3 * e + 0] += divu * ue[0];
        rhs[3 * e + 1] += divu * ue[1];
        rhs[3 * e + 2] += divu * ue[2];
    }
}

/* Matrix-free multiplication by the consistent continuous P1 mass matrix. */
static void dg_cg_mass_apply(const DGTriangle *tri, int nt,
                             const double *x, double *y, int nnodes)
{
    int e;

    memset(y, 0, (size_t)nnodes * sizeof(double));

    for (e = 0; e < nt; ++e) {
        const DGTriangle *t = &tri[e];
        int n0 = t->node[0], n1 = t->node[1], n2 = t->node[2];
        double c = t->area / 12.0;
        double x0 = x[n0], x1 = x[n1], x2 = x[n2];

        y[n0] += c * (2.0 * x0 + x1 + x2);
        y[n1] += c * (x0 + 2.0 * x1 + x2);
        y[n2] += c * (x0 + x1 + 2.0 * x2);
    }
}

/*
 * Consistent L2 projection:
 *
 *     M_CG h_CG = b_DG,
 *     (b_DG)_i = sum_K integral_K phi_i h_DG dK.
 *
 * The local DG-to-CG right-hand-side integral is evaluated exactly.
 * PCG solves the global consistent mass-matrix system.
 */
static int dg_project_to_nodes_L2(FEDGSPM s, const DGTriangle *tri, int nt,
                                  const double *U, double *h)
{
    int e, a, it, n = s->nnodes;
    double *b = NULL, *diag = NULL, *r = NULL, *z = NULL;
    double *p = NULL, *Ap = NULL;
    double rz_old, bnorm2, tol, denom;
    int maxit;

    b = (double *)calloc((size_t)n, sizeof(double));
    diag = (double *)calloc((size_t)n, sizeof(double));
    r = (double *)calloc((size_t)n, sizeof(double));
    z = (double *)calloc((size_t)n, sizeof(double));
    p = (double *)calloc((size_t)n, sizeof(double));
    Ap = (double *)calloc((size_t)n, sizeof(double));

    if (!b || !diag || !r || !z || !p || !Ap) {
        free(b); free(diag); free(r); free(z); free(p); free(Ap);
        return 1;
    }

    for (e = 0; e < nt; ++e) {
        const DGTriangle *t = &tri[e];
        int n0 = t->node[0], n1 = t->node[1], n2 = t->node[2];
        const double *ue = &U[3 * e];
        double c = t->area / 12.0;

        b[n0] += c * (2.0 * ue[0] + ue[1] + ue[2]);
        b[n1] += c * (ue[0] + 2.0 * ue[1] + ue[2]);
        b[n2] += c * (ue[0] + ue[1] + 2.0 * ue[2]);

        diag[n0] += t->area / 6.0;
        diag[n1] += t->area / 6.0;
        diag[n2] += t->area / 6.0;
    }

    memset(h, 0, (size_t)n * sizeof(double));
    memcpy(r, b, (size_t)n * sizeof(double));

    bnorm2 = 0.0;
    for (a = 0; a < n; ++a) {
        bnorm2 += b[a] * b[a];

        if (diag[a] <= 0.0 || !isfinite(diag[a])) {
            free(b); free(diag); free(r); free(z); free(p); free(Ap);
            return 1;
        }

        z[a] = r[a] / diag[a];
        p[a] = z[a];
    }

    if (bnorm2 <= 1.0e-60) {
        free(b); free(diag); free(r); free(z); free(p); free(Ap);
        return 0;
    }

    rz_old = 0.0;
    for (a = 0; a < n; ++a)
        rz_old += r[a] * z[a];

    tol = 1.0e-12 * sqrt(bnorm2) + 1.0e-14;
    maxit = (n > 100) ? 10 * n : 1000;

    for (it = 0; it < maxit; ++it) {
        double rnorm2 = 0.0, alpha, beta, rz_new;

        dg_cg_mass_apply(tri, nt, p, Ap, n);

        denom = 0.0;
        for (a = 0; a < n; ++a)
            denom += p[a] * Ap[a];

        if (denom <= 0.0 || !isfinite(denom)) {
            free(b); free(diag); free(r); free(z); free(p); free(Ap);
            return 1;
        }

        alpha = rz_old / denom;

        for (a = 0; a < n; ++a) {
            h[a] += alpha * p[a];
            r[a] -= alpha * Ap[a];
            rnorm2 += r[a] * r[a];
        }

        if (sqrt(rnorm2) <= tol)
            break;

        for (a = 0; a < n; ++a)
            z[a] = r[a] / diag[a];

        rz_new = 0.0;
        for (a = 0; a < n; ++a)
            rz_new += r[a] * z[a];

        beta = rz_new / rz_old;

        for (a = 0; a < n; ++a)
            p[a] = z[a] + beta * p[a];

        rz_old = rz_new;
    }

    if (it == maxit) {
        fprintf(stderr,
                "DG L2 projection: PCG did not converge in %d iterations\n",
                maxit);
        free(b); free(diag); free(r); free(z); free(p); free(Ap);
        return 1;
    }

    free(b); free(diag); free(r); free(z); free(p); free(Ap);
    return 0;
}

ErrorCode FEDGSPMUpdate_DGAdvection_P1_SSPRK2(FEDGSPM s, double dt,
                                               double h[])
{
    DGTriangle *tri = NULL;
    int nt = 0, e, a;
    double *u0 = NULL, *u1 = NULL, *u2 = NULL, *r = NULL;

    if (!s || !h || !s->u_h || !isfinite(dt) || dt <= 0.0)
        return 1;

    if (dg_build_mesh(s, &tri, &nt))
        return 1;

    u0 = (double *)calloc((size_t)3 * (size_t)nt, sizeof(double));
    u1 = (double *)calloc((size_t)3 * (size_t)nt, sizeof(double));
    u2 = (double *)calloc((size_t)3 * (size_t)nt, sizeof(double));
    r  = (double *)calloc((size_t)3 * (size_t)nt, sizeof(double));

    if (!u0 || !u1 || !u2 || !r) {
        free(tri); free(u0); free(u1); free(u2); free(r);
        return 1;
    }

    /* Continuous nodal values are injected into each local DG-P1 element. */
    for (e = 0; e < nt; ++e)
        for (a = 0; a < 3; ++a)
            u0[3 * e + a] = h[tri[e].node[a]];

    /* SSP-RK2; the same frozen-in-time spatial velocity is used twice. */
    dg_rhs(s, tri, nt, u0, r);
    for (e = 0; e < 3 * nt; ++e)
        u1[e] = u0[e] + dt * r[e];

    dg_rhs(s, tri, nt, u1, r);
    for (e = 0; e < 3 * nt; ++e)
        u2[e] = 0.5 * u0[e] + 0.5 * (u1[e] + dt * r[e]);

    if (dg_project_to_nodes_L2(s, tri, nt, u2, h)) {
        free(tri); free(u0); free(u1); free(u2); free(r);
        return 1;
    }

    free(tri); free(u0); free(u1); free(u2); free(r);
    return 0;
}

ErrorCode FEDGSPMGetStableTimestep_DG_P1(FEDGSPM s, double cfl,
                                         double *dt)
{
    DGTriangle *tri = NULL;
    int nt = 0, e, f;
    double dmin = DBL_MAX;
    int has_advection = 0;

    if (!s || !dt || !isfinite(cfl) || cfl <= 0.0)
        return 1;

    if (dg_build_mesh(s, &tri, &nt))
        return 1;

    for (e = 0; e < nt; ++e) {
        double denom = 0.0;

        for (f = 0; f < 3; ++f) {
            int a = f, b = (f + 1) % 3;
            double dx = tri[e].x[b] - tri[e].x[a];
            double dy = tri[e].y[b] - tri[e].y[a];
            double len = sqrt(dx * dx + dy * dy);
            double nx = dy / len;
            double ny = -dx / len;
            double ux0 = s->u_h[2 * tri[e].node[a]];
            double uy0 = s->u_h[2 * tri[e].node[a] + 1];
            double ux1 = s->u_h[2 * tri[e].node[b]];
            double uy1 = s->u_h[2 * tri[e].node[b] + 1];
            double bn0 = fabs(ux0 * nx + uy0 * ny);
            double bn1 = fabs(ux1 * nx + uy1 * ny);
            double bn = (bn0 > bn1) ? bn0 : bn1;

            denom += bn * len;
        }

        if (denom > 1.0e-20) {
            double de = cfl * (2.0 * tri[e].area) / (3.0 * denom);
            if (de < dmin)
                dmin = de;
            has_advection = 1;
        }
    }

    free(tri);

    /* With identically zero velocity, there is no advective restriction. */
    *dt = has_advection ? dmin : DBL_MAX;
    return 0;
}
