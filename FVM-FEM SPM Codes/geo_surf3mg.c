/* COUPLING I3ELVIS thermomechanical code with FEFVSPM, commit  381c06d */
/* incremental change to have true change rates due to advection, uplift, and evolution shown separate */

#ifdef MAINFILE
#undef MAINFILE
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#include "geo_head3mg.h"

/* LAPACK dense linear solver used by planeq(). */
extern void dgesv_(int *n, int *nrhs, double *a, int *lda,
                   int *ipiv, double *b, int *ldb, int *info);

/* ------------------------------------------------------------------- */
/* Calculate surface positions from contact points external code */
int surface()
{

    ErrorCode ierr;
    /* Internal counters and service */
    int n;
    int status = 1; /* catch whether SPM returns successfully */
    int sfc_num_max = 7500;

    /* calculate velocity field for SPM nodes*/
    // veloc_sfc();

    /* set times : reset start to zero  (only what SPM sees)*/
    sfc_time = 0;
    spm->time = 0.0;
    sfc_tfinal = timestep / 3.15576e+7;
    // spm->time_end = timestep;
    sfc_tstp = timestep / sfc_num_tstp;
    // spm->dt = timestep/sfc_num_tstp;

    ierr = FEFVSPMSetEndTime(spm, timestep);
    ierr = FEFVSPMSetTimestepBounds(spm, 1e-6, sfc_tstp);
    ierr = FEFVSPMSetTimestep(spm, sfc_tstp);
    ierr = FEFVSPMSetConstantTimestep(spm, sfc_tstp);

    if (sfcmod >= 4 && sfcmod <= 6 && spm->fespm != NULL)
    {
        double dt_fem_stable;
        ierr = FESPMGetStableTimestep(spm->fespm, &dt_fem_stable);
        if (ierr != 0)
        {
            printf("FATAL: FESPMGetStableTimestep failed in surface\n");
            exit(EXIT_FAILURE);
        }
        if (dt_fem_stable > 0.0 && spm->dt > dt_fem_stable)
        {
            printf("FESPM: reducing surface timestep from %.6e to %.6e\n",
                   spm->dt, dt_fem_stable);
            spm->dt = dt_fem_stable;
            ierr = FEFVSPMSetConstantTimestep(spm, spm->dt);
        }
    }



    // if(spm->dt<1e-6 || isnan(spm->dt))spm->dt=1e-6;
    fflush(stdout);

    // ----
    // Surface internal loop, if desired
    for (n = 0; (n < sfc_num_max && spm->time < spm->time_end && spm->dt > 0); n++)
    {

        if (spm->time_end - spm->time < spm->dt && (spm->time_end - spm->time) / spm->dt > 1e-9)
        {
            spm->dt = spm->time_end - spm->time;
            ierr = FEFVSPMSetConstantTimestep(spm, (spm->time_end - spm->time));
        }
        else if (n >= sfc_num_max)
        {
            printf("MAXIMUM NUMBER OF SURFACE TIME STEPS EXCEEDED. EXITING\n");
            exit(EXIT_FAILURE);
        }

        /* FEM explicit-diffusion stability limit */
        if (sfcmod >= 4 && sfcmod <= 6 && spm->fespm != NULL)
        {
            double dt_fem_stable;
            ierr = FESPMGetStableTimestep(spm->fespm, &dt_fem_stable);
            if (ierr != 0)
            {
                printf("FATAL: FESPMGetStableTimestep failed in surface loop\n");
                exit(EXIT_FAILURE);
            }
            if (dt_fem_stable > 0.0 && spm->dt > dt_fem_stable)
            {
                printf("FESPM: reducing surface timestep from %.6e to %.6e\n",
                       spm->dt, dt_fem_stable);
                spm->dt = dt_fem_stable;
                ierr = FEFVSPMSetConstantTimestep(spm, spm->dt);
            }
        }



        /* calculate velocity field for SPM nodes*/
        veloc_sfc();
        /* Enforce the FVM CFL using the velocity field just computed. */
        if (sfcmod >= 1 && sfcmod <= 6)
        {
            double dt_adv_stable;
            ierr = FEFVSPMGetStableTimestep_FVMAdvection(spm, FEFVSPM_FVM_ADVECTION_CFL, &dt_adv_stable);
            if (ierr != 0)
            {
                printf("FATAL: FEFVSPMGetStableTimestep_FVMAdvection failed in surface loop\n");
                exit(EXIT_FAILURE);
            }
            if (dt_adv_stable > 0.0 && spm->dt > dt_adv_stable)
            {
                printf("FVM: reducing surface timestep from %.6e to %.6e\n",
                       spm->dt, dt_adv_stable);
                spm->dt = dt_adv_stable;
                ierr = FEFVSPMSetConstantTimestep(spm, spm->dt);
            }
        }
        // uplift and advect
        move_sfc();
        /* this is the culprit - write everything evolution here */
        evolve_sfc();

        spm->time += spm->dt;
        if ((spm->time_end - spm->time) / spm->dt < 1e-9)
            spm->time = spm->time_end;
        printf("time %8.3f a\n", spm->time / 3.15576e+7);

        // if(n==1)write_sfc(1);
    }
    // ----

    if (n >= sfc_num_max - 1)
    {
        printf("Diffusion not resolved within %i timesteps (spm_dt: 4%3f, thm_dt: %4.3f) - EXITING\n", sfc_num_max, spm->dt / 3.15576e+7, sfc_tfinal);
        exit(-666);
    };

    printf("%i surface cycles ...\n", n);
    // write_sfc(1);

    return_sfc();

    status = 0;
    return status;
}
/* END calculate surface positions from external code*/
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* Loop through all surface nodes, assign them a velocity */
void veloc_sfc()
{
    long int mm1, mm2, mm3, mm4, ctr;
    long int wn1[500];
    double vxyz1[10];
    double x, z;

    ErrorCode ierr;

    for (mm2 = 0; mm2 < nsfcx; mm2++) // loop all sfc nodes
        for (mm3 = 0; mm3 < nsfcz; mm3++)
        {
            mm1 = mm2 * nsfcz + mm3;
            mm4 = mm3 * nsfcx + mm2;

            x = mm2 * xsstp;
            z = mm3 * zsstp;

            sfc_u[mm4] = 0;
            sfc_v[mm4] = 0;
            sfc_w[mm4] = 0;

            if (mm2 == 0)
            {
                vxyzcalc1(xstpx, waterlev - spm->h[mm4] + ystpy * 0.75, z, vxyz1, wn1);
            }
            else if (mm2 == nsfcx - 1)
            {
                vxyzcalc1(xsize - xstpx, waterlev - spm->h[mm4] + ystpy * 0.75, z, vxyz1, wn1);
            }
            else if (mm3 == 0)
            {
                vxyzcalc1(x, waterlev - spm->h[mm4] + ystpy * 0.75, zstpz, vxyz1, wn1);
            }
            else if (mm3 == nsfcz - 1)
            {
                vxyzcalc1(x, waterlev - spm->h[mm4] + ystpy * 0.75, zsize - zstpz, vxyz1, wn1);
            }
            else
            {
                vxyzcalc1(x, waterlev - spm->h[mm4] + ystpy * 0.75, z, vxyz1, wn1);
            }

            if (isinf(vxyz1[0]) || isnan(vxyz1[0]) || isinf(vxyz1[1]) || isnan(vxyz1[1]) || isinf(vxyz1[2]) || isnan(vxyz1[2]))
            {
                printf("invalid velocity \n");
                exit(-1);
            }

            sfc_u[mm4] = vxyz1[0];
            sfc_v[mm4] = vxyz1[2];
            sfc_w[mm4] = vxyz1[1];

        } // end for all surface nodes

    /* set uplift */
    {
        double *uvert;

        uvert = malloc(sizeof(double) * spm->nnodes);
        NodalIteratorBegin(spm) double xn[2];
        int nidx;

        nidx = NodeIdx(spm, ni, nj);
        FEFVSPMGetNodeCoordinate(spm, ni, nj, xn);
        uvert[nidx] = -sfc_w[nidx];
        NodalIteratorEnd(spm)

            ierr = FEFVSPMSetVerticalVelocity(spm, uvert);
        free(uvert);
    }
    /* set  horiz vel */
    {
        double *uhoriz;

        uhoriz = malloc(sizeof(double) * 2 * spm->nnodes);
        NodalIteratorBegin(spm) double xn[2];
        int nidx;

        nidx = NodeIdx(spm, ni, nj);
        FEFVSPMGetNodeCoordinate(spm, ni, nj, xn);

        uhoriz[2 * nidx] = sfc_u[nidx];
        uhoriz[2 * nidx + 1] = sfc_v[nidx];

        NodalIteratorEnd(spm)

            ierr = FEFVSPMSetHorizontalVelocity(spm, uhoriz);
        free(uhoriz);
    }
}
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* Loop through all fine, moving surface nodes, advect them */
void move_sfc()
{
    long int mm1, mm2, mm3, mm4, ctr;
    ErrorCode ierr;

    for (mm2 = 0; mm2 < nsfcx; mm2++) // loop all sfc nodes
        for (mm3 = 0; mm3 < nsfcz; mm3++)
        {
            mm1 = mm2 * nsfcz + mm3;
            mm4 = mm3 * nsfcx + mm2;
            sfc_old[mm1] = sfc[mm1] = spm->h[mm4];
            sfc_upl[mm1] = spm->dt * spm->u_v[mm4];
        }

    if (sfcmod < 4)
    {
        /* Preserve the original SL ordering for sfcmod 1/2/3:
           apply uplift first, then advect the uplifted surface with FVM. */
        for (mm2 = 0; mm2 < nsfcx; mm2++)
            for (mm3 = 0; mm3 < nsfcz; mm3++)
            {
                mm1 = mm2 * nsfcz + mm3;
                mm4 = mm3 * nsfcx + mm2;
                spm->h[mm4] += spm->dt * spm->u_v[mm4];
                sfc[mm1] += spm->dt * spm->u_v[mm4];
            }

        /* Second-order MUSCL-Hancock FVM advection. */
        ierr = FEFVSPMUpdate_FVMAdvection_2ndOrder(spm, spm->dt, spm->h);
        if (ierr != 0)
        {
            printf("FATAL: FEFVSPMUpdate_FVMAdvection_2ndOrder (h) failed in move_sfc\n");
            exit(EXIT_FAILURE);
        }

        /* Record the advection increment relative to the uplifted surface. */
        for (mm2 = 0; mm2 < nsfcx; mm2++)
            for (mm3 = 0; mm3 < nsfcz; mm3++)
            {
                mm1 = mm2 * nsfcz + mm3;
                mm4 = mm3 * nsfcx + mm2;
                sfc_adv[mm1] = spm->h[mm4] - sfc[mm1];
                sfc[mm1] = spm->h[mm4];
            }
    }
    else
    {
        /* sfcmod 4/5/6: FVM advection followed by FEM diffusion/uplift. */
        spm->dtmode = DTM_CONST;
        ierr = FEFVSPMSetConstantTimestep(spm, spm->dt);

        ierr = FEFVSPMUpdate_FVMAdvection_2ndOrder(spm, spm->dt, spm->h);
        if (ierr != 0)
        {
            printf("FATAL: FEFVSPMUpdate_FVMAdvection_2ndOrder (h) failed in move_sfc\n");
            exit(EXIT_FAILURE);
        }

        /* Advect spatially varying diffusivity parameters where required. */
        if (sfcmod > 4) {
            ierr = FEFVSPMUpdate_FVMAdvection_2ndOrder(spm, spm->dt, spm->k0_n);
            if (ierr != 0) {
                printf("FATAL: FEFVSPMUpdate_FVMAdvection_2ndOrder (k0_n) failed in move_sfc\n");
                exit(EXIT_FAILURE);
            }
        }

        for (mm2 = 0; mm2 < nsfcx; mm2++)
            for (mm3 = 0; mm3 < nsfcz; mm3++)
            {
                mm1 = mm2 * nsfcz + mm3;
                mm4 = mm3 * nsfcx + mm2;
                sfc_adv[mm1] = spm->h[mm4] - sfc[mm1];
            }
    }

}
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* Evolve surface grid according to anything */
void evolve_sfc()
{

    long int mm1, mm2, mm3, mm4;
    int j;
    double x, z, exc, dh;
    ErrorCode ierr;

    switch (sfcmod)
    {
    case 1:
        for (mm2 = 0; mm2 < nsfcx; mm2++)
            for (mm3 = 0; mm3 < nsfcz; mm3++)
            {
                mm1 = mm2 * nsfcz + mm3;
                x = mm2 * xsstp;
                z = mm3 * zsstp;
                if (waterlev - sfc[mm1] < eroslev)
                    sfc[mm1] = waterlev - eroslev;
                if (waterlev - sfc[mm1] > sedilev)
                    sfc[mm1] = waterlev - sedilev;
            }
        break;
    case 2:
        for (mm2 = 0; mm2 < nsfcx; mm2++)
            for (mm3 = 0; mm3 < nsfcz; mm3++)
            {
                mm1 = mm2 * nsfcz + mm3;
                mm4 = mm3 * nsfcx + mm2;
                x = mm2 * xsstp;
                z = mm3 * zsstp;
                if (waterlev - sfc[mm1] < eroslev)
                {
                    sfc[mm1] = MAXV(waterlev - eroslev, sfc[mm1] - eroscon * sfc_tstp);
                }
                if (waterlev - sfc[mm1] > sedilev)
                {
                    sfc[mm1] = MINV(waterlev - sedilev, sfc[mm1] + sedicon * sfc_tstp);
                }
                dh = sfc[mm1] - spm->h[mm4];
                mm4 = mm1;
            }
        break;
    case 3:
        for (mm2 = 0; mm2 < nsfcx; mm2++)
            for (mm3 = 0; mm3 < nsfcz; mm3++)
            {
                mm1 = mm2 * nsfcz + mm3;
                mm4 = mm3 * nsfcx + mm2;
                x = mm2 * xsstp;
                z = mm3 * zsstp;
                if (waterlev - sfc[mm1] < eroslev)
                {
                    exc = eroslev - waterlev + sfc[mm1];
                    sfc[mm1] = MAXV(waterlev - eroslev, sfc[mm1] - eroscon * sfc_tstp * exc);
                }
                if (waterlev - sfc[mm1] > sedilev)
                {
                    exc = waterlev - sedilev - sfc[mm1];
                    sfc[mm1] = MINV(waterlev - sedilev, sfc[mm1] + sedicon * sfc_tstp * exc);
                }
            }
        break;
    case 4:
        ierr = FEFVSPMSetDiffusionMode(spm, DM_CONST);
        diffuse_sfc();
        break;
    case 5:
        ierr = FEFVSPMSetDiffusionMode(spm, DM_NONLINEAR_POWERLAW);
        diffuse_sfc();
        break;
    case 6:
        ierr = FEFVSPMSetDiffusionMode(spm, DM_NONLINEAR_CRITICALSLOPE);
        diffuse_sfc();
        break;
    }

    if (sfcmod < 4)
    {
        for (mm2 = 0; mm2 < nsfcx; mm2++)
            for (mm3 = 0; mm3 < nsfcz; mm3++)
            {
                mm1 = mm2 * nsfcz + mm3;
                mm4 = mm3 * nsfcx + mm2;
                spm->h[mm4] = sfc[mm1];
                spm->h_c[mm4] = sfc_old[mm1];

                if (sfcmod < 4)
                    sfc_evl[mm1] = sfc[mm1] - (sfc_old[mm1] + sfc_adv[mm1] + sfc_upl[mm1]);
            }
    }
}
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* FEM diffusion for sfcmod 4/5/6 */
void diffuse_sfc()
{
    ErrorCode ierr;
    int mm1, mm2, mm3, mm4;

    /* FVM advection has already been applied in move_sfc().
     * h_c stores the post-advection field h*.
     * FEM applies diffusion and vertical velocity:
     * M H^(n+1) = M H* + dt (M U_V - K H*).
     * Homogeneous Neumann conditions are natural FEM boundary conditions. */
    memcpy(spm->h_c, spm->h, sizeof(double) * spm->nnodes);

    ierr = FESPMDiffusionStep(spm->fespm, spm, spm->h_c);
    if (ierr != 0)
    {
        printf("FATAL: FESPMDiffusionStep failed in diffuse_sfc\n");
        exit(EXIT_FAILURE);
    }

    for (mm2 = 0; mm2 < nsfcx; mm2++)
        for (mm3 = 0; mm3 < nsfcz; mm3++)
        {
            mm1 = mm2 * nsfcz + mm3;
            mm4 = mm3 * nsfcx + mm2;
            sfc_evl[mm1] = spm->h[mm4] - spm->h_c[mm4] - spm->dt * spm->u_v[mm4];
        }

    memcpy(spm->h_c, spm->h, sizeof(double) * spm->nnodes);
}
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* Transfer elevation from FEFVSPM to ELVIS */
void return_sfc()
{
    long int mm1, mm2, mm3, mm4;

    for (mm2 = 0; mm2 < nsfcx; mm2++)
        for (mm3 = 0; mm3 < nsfcz; mm3++)
        {
            mm4 = mm3 * nsfcx + mm2;
            mm1 = mm2 * nsfcz + mm3;
            sfc[mm1] = spm->h_c[mm4];
            kvar[mm4] = spm->k0_n[mm4];
        }
}
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* Write fields of interest */
void write_sfc(int f0)
{
    long int mm1, mm2, mm3, mm4, mm5, mm6, mm7;
    int i, stll;

    char prefx[50];
    char filesc[50] = {};
    char fls4[] = "sfc_";
    char postf[] = ".vtk";

    int si = 0;
    char cr = 'a';
    while ((cr != '\0') && (si < 50))
    {
        prefx[si] = fl1out[si];
        cr = fl1out[si];
        si++;
    }
    sprintf(filesc, "%s", fls4);
    stll = strlen(prefx);
    for (i = 0; i < stll - 4; i++)
    {
        filesc[i + 4] = prefx[i];
    }
    strncat(filesc, postf, 50);

    fl = fopen(filesc, "wt");

    fprintf(fl, "# vtk DataFile Version 3.0\nFEFVSPM surface\nASCII\nDATASET UNSTRUCTURED_GRID\n");
    fprintf(fl, "POINTS %i float\n", nsfc);
    for (mm2 = 0; mm2 < nsfcx; mm2++)
    {
        for (mm3 = 0; mm3 < nsfcz; mm3++)
        {
            mm1 = mm2 * nsfcz + mm3;
            mm4 = mm3 * nsfcx + mm2;
            fprintf(fl, "%6.3f %6.3f %6.3f\n", spm->dx * mm2, waterlev - spm->h_c[mm4], spm->dy * mm3);
        }
    }

    fprintf(fl, "CELLS %i %i\n", (nsfcx - 1) * (nsfcz - 1), (nsfcx - 1) * (nsfcz - 1) * 5);
    for (mm2 = 0; mm2 < nsfcx - 1; mm2++)
    {
        for (mm3 = 0; mm3 < nsfcz - 1; mm3++)
        {
            mm4 = mm2 * nsfcz + mm3;
            mm6 = mm4 + 1;
            mm5 = mm4 + nsfcz;
            mm7 = mm4 + nsfcz + 1;
            fprintf(fl, "4 %i %i %i %i \n", mm4, mm5, mm6, mm7);
        }
    }

    fprintf(fl, "CELL_TYPES %i \n", (nsfcx - 1) * (nsfcz - 1));
    for (mm2 = 0; mm2 < (nsfcx - 1) * (nsfcz - 1); mm2++)
    {
        fprintf(fl, "8 \n");
    }

    fprintf(fl, "POINT_DATA %i\nSCALARS elevation float 1\nLOOKUP_TABLE default\n", nsfc);
    for (mm2 = 0; mm2 < nsfcx; mm2++)
    {
        for (mm3 = 0; mm3 < nsfcz; mm3++)
        {
            mm4 = mm3 * nsfcx + mm2;
            fprintf(fl, "%6.6f\n", spm->h_c[mm4]);
        }
    }

    fprintf(fl, "SCALARS horz_vel_x float 1\nLOOKUP_TABLE default\n");
    for (mm2 = 0; mm2 < nsfcx; mm2++)
    {
        for (mm3 = 0; mm3 < nsfcz; mm3++)
        {
            mm4 = mm3 * 2 * nsfcx + mm2 * 2;
            fprintf(fl, "%6.3e\n", spm->u_h[mm4]);
        }
    }

    fprintf(fl, "SCALARS horz_vel_z float 1\nLOOKUP_TABLE default\n");
    for (mm2 = 0; mm2 < nsfcx; mm2++)
    {
        for (mm3 = 0; mm3 < nsfcz; mm3++)
        {
            mm4 = mm3 * 2 * nsfcx + mm2 * 2 + 1;
            fprintf(fl, "%6.3e\n", spm->u_h[mm4]);
        }
    }

    fprintf(fl, "SCALARS dhdt_advection float 1\nLOOKUP_TABLE default\n");
    for (mm2 = 0; mm2 < nsfcx; mm2++)
    {
        for (mm3 = 0; mm3 < nsfcz; mm3++)
        {
            mm4 = mm3 * nsfcx + mm2;
            mm1 = mm2 * nsfcz + mm3;
            if (spm->dt > 0) fprintf(fl, "%6.3e\n", (sfc_adv[mm1]) / spm->dt);
            else fprintf(fl, "%6.3e\n", 0.0);
        }
    }
    fprintf(fl, "SCALARS dhdt_uplift float 1\nLOOKUP_TABLE default\n");
    for (mm2 = 0; mm2 < nsfcx; mm2++)
    {
        for (mm3 = 0; mm3 < nsfcz; mm3++)
        {
            mm4 = mm3 * nsfcx + mm2;
            mm1 = mm2 * nsfcz + mm3;
            if (spm->dt > 0) fprintf(fl, "%6.3e\n", (sfc_upl[mm1]) / spm->dt);
            else fprintf(fl, "%6.3e\n", 0.0);
        }
    }
    fprintf(fl, "SCALARS dhdt_ers_sed float 1\nLOOKUP_TABLE default\n");
    for (mm2 = 0; mm2 < nsfcx; mm2++)
    {
        for (mm3 = 0; mm3 < nsfcz; mm3++)
        {
            mm4 = mm3 * nsfcx + mm2;
            mm1 = mm2 * nsfcz + mm3;
            if (spm->dt > 0) fprintf(fl, "%6.3e\n", (sfc_evl[mm1]) / spm->dt);
            else fprintf(fl, "%6.3e\n", 0.0);
        }
    }

    fprintf(fl, "SCALARS base_diffusivity float 1\nLOOKUP_TABLE default\n");
    for (mm2 = 0; mm2 < nsfcx; mm2++)
    {
        for (mm3 = 0; mm3 < nsfcz; mm3++)
        {
            mm4 = mm3 * nsfcx + mm2;
            fprintf(fl, "%6.3e\n", spm->k0_n[mm4]);
        }
    }

    fprintf(fl, "\nCELL_DATA %i\nSCALARS eff_cell_diffusivity float 1\nLOOKUP_TABLE default\n", spm->ncells);
    for (mm2 = 0; mm2 < nsfcx - 1; mm2++)
    {
        for (mm3 = 0; mm3 < nsfcz - 1; mm3++)
        {
            mm4 = mm3 * (nsfcx - 1) + mm2;
            fprintf(fl, "%6.3e\n", spm->kappa[mm4]);
        }
    }

    fclose(fl);
}
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* Initialize level surface on grid domain */
int init_surface()
{
    int status = 1;
    int km_midi, km_begi, km_endi;
    double km_fact;
    int mm1, mm2, mm3, mm4;
    long int i, j, ij;
    char bcon[4] = "NNNN";

    ErrorCode ierr;

    nsfcx = (int)ceil((xnumx - 1) * SFC2G + 1);
    nsfcz = (int)ceil((znumz - 1) * SFC2G + 1);
    nsfc = nsfcx * nsfcz;
    xsstp = xsize / (nsfcx - 1);
    zsstp = zsize / (nsfcz - 1);

    ierr = FEFVSPMCreate(&spm);

    ierr = FEFVSPMSetGeometry(spm, 0, xsize, 0, zsize);
    ierr = FEFVSPMSetResolution(spm, nsfcx, nsfcz);

    ierr = FEFVSPMSetMaxNumberSteps(spm, 500);
    ierr = FEFVSPMSetEndTime(spm, 1e-3);
    ierr = FEFVSPMSetTimestepMode(spm, DTM_CONST);
    ierr = FEFVSPMSetTimestepIO(spm, 10);

    ierr = FEFVSPMSetBoundaryConditions(spm, bcon);

    spm->k0 = sfc_diff_k;
    spm->d_exp = sfc_diff_exp;
    spm->Sc = tan(sfc_diff_threshold * 3.141592653 / 180.0);

    if (sfcmod == 4)
        ierr = FEFVSPMSetConstantDiffusivity(spm, sfc_diff_k);
    if (sfcmod == 5)
        ierr = FEFVSPMSetNonlinearDiffusivityPowerLaw(spm, sfc_diff_k, sfc_diff_exp);
    if (sfcmod == 6)
        ierr = FEFVSPMSetNonlinearDiffusivityCriticalSlope(spm, sfc_diff_k, tan(sfc_diff_threshold * 3.141592653 / 180.0), sfc_diff_exp);

    if (sfcmod >= 4 && sfcmod <= 6)
    {
        ierr = FESPMCreate(spm, &spm->fespm);
        if (ierr != 0)
        {
            printf("FATAL: FESPMCreate failed in init_surface\n");
            exit(EXIT_FAILURE);
        }
    }

    km_fact = 1.0;
    km_midi = round(spm->ny / 2);
    km_begi = km_midi - round(0.05 * zsize / zsstp);
    km_endi = km_midi + round(0.05 * zsize / zsstp);

    NodalIteratorBegin(spm) if (nj < km_begi)
    {
        spm->k0_n[NodeIdx(spm, ni, nj)] = sfc_diff_k;
    }
    else if (nj < km_endi + 1)
    {
        spm->k0_n[NodeIdx(spm, ni, nj)] = sfc_diff_k + sfc_diff_k * (km_fact - 1) * (nj - km_begi) / (km_endi - km_begi);
    }
    else
    {
        spm->k0_n[NodeIdx(spm, ni, nj)] = sfc_diff_k * km_fact;
    }
    NodalIteratorEnd(spm)

        NodalIteratorBegin(spm)
            spm->Sc_n[NodeIdx(spm, ni, nj)] = tan(sfc_diff_threshold * 3.141592653 / 180.0);
    spm->d_exp_n[NodeIdx(spm, ni, nj)] = sfc_diff_exp;
    NodalIteratorEnd(spm)

        printf("surface nodes: %i x %i, total: %i\n", nsfcx, nsfcz, nsfc);

    for (j = 0; j < nsfcz; j++)
    {
        for (i = 0; i < nsfcx; i++)
        {
            ij = i * nsfcz + j;
            if ((i == 0) || (i == (nsfcx - 1)) || (j == 0) || (j == (nsfcz - 1)))
            {
                sfc[ij] = 0;
            }
            else
            {
                sfc[ij] += 0.5 * (2 * (double)(rand()) / RAND_MAX);
            }
            sfc_u[ij] = 0;
            sfc_v[ij] = 0;
            sfc_w[ij] = 0;
        }
    }

    {
        double *hic;
        hic = malloc(sizeof(double) * spm->nnodes);
        for (mm2 = 0; mm2 < nsfcx; mm2++)
            for (mm3 = 0; mm3 < nsfcz; mm3++)
            {
                mm4 = mm3 * nsfcx + mm2;
                mm1 = mm2 * nsfcz + mm3;
                hic[mm4] = sfc[mm1];
            }
        ierr = FEFVSPMSetInitialHeight(spm, hic);
        free(hic);
    }

    {
        double *uvert;
        double *uhoriz;
        uvert = malloc(sizeof(double) * spm->nnodes);
        uhoriz = malloc(sizeof(double) * 2 * spm->nnodes);
        NodalIteratorBegin(spm) double xn[2];
        int nidx;
        nidx = NodeIdx(spm, ni, nj);
        FEFVSPMGetNodeCoordinate(spm, ni, nj, xn);
        uvert[nidx] = 0.0;
        uhoriz[2 * nidx] = 0.0;
        uhoriz[2 * nidx + 1] = 0.0;
        NodalIteratorEnd(spm)
            ierr = FEFVSPMSetVerticalVelocity(spm, uvert);
        free(uvert);
        ierr = FEFVSPMSetHorizontalVelocity(spm, uhoriz);
        free(uhoriz);
    }

    ierr = FEFVSPMStepInitialize(spm);

    for (j = 0; j < spm->nnodes; j++)
    {
        spm->h_c[j] = spm->h[j];
    }

    status = 0;
    return status;
}
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* Re-init surface with existing triangulation */
int reinit_surface()
{
    long int mm1, mm2, mm3, mm4;
    int status = 1;
    char bcon[4] = "NNNN";
    ErrorCode ierr;

    printf("reinit_surface:  nsurface nodes: %i\n", nsfc);

    nsfcx = (int)ceil((xnumx - 1) * SFC2G + 1);
    nsfcz = (int)ceil((znumz - 1) * SFC2G + 1);
    nsfc = nsfcx * nsfcz;
    xsstp = xsize / (nsfcx - 1);
    zsstp = zsize / (nsfcz - 1);

    fflush(stdout);

    ierr = FEFVSPMCreate(&spm);

    spm->k0 = sfc_diff_k;
    spm->d_exp = sfc_diff_exp;
    spm->Sc = sfc_diff_threshold;

    ierr = FEFVSPMSetGeometry(spm, 0, xsize, 0, zsize);
    ierr = FEFVSPMSetResolution(spm, nsfcx, nsfcz);

    if (sfcmod == 4)
        ierr = FEFVSPMSetConstantDiffusivity(spm, sfc_diff_k);
    if (sfcmod == 5)
        ierr = FEFVSPMSetNonlinearDiffusivityPowerLaw(spm, sfc_diff_k, sfc_diff_exp);
    if (sfcmod == 6)
        ierr = FEFVSPMSetNonlinearDiffusivityCriticalSlope(
            spm, sfc_diff_k,
            tan(sfc_diff_threshold * 3.141592653 / 180.0),
            sfc_diff_exp);

    if (sfcmod >= 4 && sfcmod <= 6)
    {
        ierr = FESPMCreate(spm, &spm->fespm);
        if (ierr != 0)
        {
            printf("FATAL: FESPMCreate failed in reinit_surface\n");
            exit(EXIT_FAILURE);
        }
    }

    ierr = FEFVSPMSetMaxNumberSteps(spm, 500);
    ierr = FEFVSPMSetEndTime(spm, 80.04e-3);
    ierr = FEFVSPMSetTimestepMode(spm, DTM_CONST);
    ierr = FEFVSPMSetTimestepIO(spm, 10);

    ierr = FEFVSPMSetBoundaryConditions(spm, bcon);

    NodalIteratorBegin(spm)
        spm->Sc_n[NodeIdx(spm, ni, nj)] = tan(sfc_diff_threshold * 3.141592653 / 180.0);
    spm->d_exp_n[NodeIdx(spm, ni, nj)] = sfc_diff_exp;
    NodalIteratorEnd(spm)

    for (mm2 = 0; mm2 < nsfcx; mm2++)
        for (mm3 = 0; mm3 < nsfcz; mm3++)
        {
            mm1 = mm2 * nsfcz + mm3;
            mm4 = mm3 * nsfcx + mm2;
            spm->h[mm4] = sfc[mm1];
            spm->h_c[mm4] = sfc[mm1];
        }

    for (mm2 = 0; mm2 < nsfcx; mm2++)
        for (mm3 = 0; mm3 < nsfcz; mm3++)
        {
            mm1 = mm2 * (nsfcz) + mm3;
            mm4 = mm3 * (nsfcx) + mm2;
            spm->k0_n[mm4] = kvar[mm4];
        }

    status = 0;
    return status;
}
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* Equation of plane */
void planeq(double x[], double y[], double z[], double coef[])
{
    int n = 3;
    int nrhs = 1;
    double G[9];
    int lda = 3;
    int ldb = 3;
    int ipiv[3];
    double b[3];
    int info;

    G[0] = 1; G[1] = 1; G[2] = 1;
    G[3] = x[0]; G[4] = x[1]; G[5] = x[2];
    G[6] = y[0]; G[7] = y[1]; G[8] = y[2];

    b[0] = z[0]; b[1] = z[1]; b[2] = z[2];

    dgesv_(&n, &nrhs, G, &lda, ipiv, b, &ldb, &info);

    if (info != 0)
    {
        printf("inversion state: %i\n", info);
    }

    coef[0] = b[0];
    coef[1] = b[1];
    coef[2] = b[2];

    return;
}

/* ------------------------------------------------------------------- */