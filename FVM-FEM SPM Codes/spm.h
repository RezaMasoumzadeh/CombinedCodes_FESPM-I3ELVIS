#ifndef __spm_h__
#define __spm_h__

#include <stdio.h>
#include <stdlib.h>

typedef int ErrorCode;

typedef enum {  DM_CONST, DM_NONLINEAR_CRITICALSLOPE, DM_NONLINEAR_POWERLAW,
                DM_UNDEFINED } FEFVSPMDiffusionMode;
typedef enum { DTM_CONST, DTM_EXPLICIT } FEFVSPMTimestepMode;

typedef enum { DIR_JMAX=0, DIR_IMAX, DIR_JMIN, DIR_IMIN } FEFVSPMDirection;

typedef enum { FVT_ASCII, FVT_VTR, FVT_VTRB, FVT_UNDEFINED } FEFVSPMFieldViewType;

#define FEFVSPM_IGNORE -1
#define FEFVSPM_DECIDE -2

/* Initial CFL coefficient for the MUSCL-Hancock FVM advection bound. */
#define FEFVSPM_FVM_ADVECTION_CFL 1.0

typedef struct _p_FEFVSPM *FEFVSPM;

typedef struct _p_FESPM *FESPM;

struct _p_FEFVSPM {
    double xrange[2];
    double yrange[2];
    double time_end,dt,dt_const,time;
    int nx,ny,mx,my,nnodes,ncells;
    double dx,dy,dt_theta;
    int nt,nt_io,nt_max;
    FEFVSPMDiffusionMode mode;
    FEFVSPMTimestepMode dtmode;
    double *h_ic,*h,*h_c,*u_v,*u_h,*k0_n,*d_exp_n,*Sc_n; /* node */
    double *kappa; /* cell */
    double k0,d_exp,Sc;
    char bctype[4];
    double dt_min,dt_max;

        /*
     * FEM diffusion module.
     */
    FESPM fespm;
};

#define FEFVSPMCheckIsIgnore(a) ( (int)((a)) == FEFVSPM_IGNORE ? 1 : 0 )

#include "spmgriditerator.h"

#define ERRQ(str) { printf("*** FEFVSPM ERROR: -------------------------- Error Message --------------------------\n"\
                           "*** FEFVSPM ERROR: %s\n"\
                           "*** FEFVSPM ERROR: %s() line %d in %s\n"\
                           "*** FEFVSPM ERROR: -------------------------- End of Error Message --------------------------\n",\
                           (str),__FUNCTION__,__LINE__,__FILE__); exit(EXIT_FAILURE); }

#define ERRQ2(str1,str2) { printf("*** FEFVSPM ERROR: -------------------------- Error Message --------------------------\n"\
                                  "*** FEFVSPM ERROR: %s %s\n"\
                                  "*** FEFVSPM ERROR: %s() line %d in %s\n"\
                                  "*** FEFVSPM ERROR: -------------------------- End of Error Message --------------------------\n",\
                                  (str1),(str2),__FUNCTION__,__LINE__,__FILE__); exit(EXIT_FAILURE); }


ErrorCode FEFVSPMCreate(FEFVSPM *s);
ErrorCode FEFVSPMDestroy(FEFVSPM *s);

ErrorCode FEFVSPMSetGeometry(FEFVSPM s,double x0,double x1,double y0,double y1);
ErrorCode FEFVSPMSetResolution(FEFVSPM s,int nx,int ny);

ErrorCode FEFVSPMSetTimestepBounds(FEFVSPM spm,double dt_min,double dt_max);
ErrorCode FEFVSPMSetTimestep(FEFVSPM s,double dt);
ErrorCode FEFVSPMSetTimestepTheta(FEFVSPM s,double theta);
ErrorCode FEFVSPMSetMaxNumberSteps(FEFVSPM s,int nt_max);
ErrorCode FEFVSPMSetEndTime(FEFVSPM s,double tend);
ErrorCode FEFVSPMSetTimestepMode(FEFVSPM s,FEFVSPMTimestepMode mode);
ErrorCode FEFVSPMSetTimestepIO(FEFVSPM s,int nt);
ErrorCode FEFVSPMSetConstantTimestep(FEFVSPM s,double dt);
ErrorCode FEFVSPMSetExplicitTimestep(FEFVSPM s,double theta);

ErrorCode FEFVSPMSetDiffusionMode(FEFVSPM s,FEFVSPMDiffusionMode mode);
ErrorCode FEFVSPMSetConstantDiffusivity(FEFVSPM s,double k0);
ErrorCode FEFVSPMSetNonlinearDiffusivityCriticalSlope(FEFVSPM s,double k0,double Sc,double d_exp);
ErrorCode FEFVSPMSetNonlinearDiffusivityPowerLaw(FEFVSPM s,double k0,double d_exp);

ErrorCode FEFVSPMSetInitialHeight(FEFVSPM s,double h0[]);
ErrorCode FEFVSPMSetVerticalVelocity(FEFVSPM s,double v[]);
ErrorCode FEFVSPMSetHorizontalVelocity(FEFVSPM s,double uw[]);

ErrorCode FEFVSPMSetBoundaryConditions(FEFVSPM s,const char bctype[]);

ErrorCode FEFVSPMGetNodeCoordinate(FEFVSPM s,int i,int j,double coor[]);
ErrorCode FEFVSPMGetCentroidCoordinate(FEFVSPM s,int i,int j,double coor[]);
ErrorCode FEFVSPMGetVerticalVelocity(FEFVSPM s,int *m,int *n,double **data);
ErrorCode FEFVSPMGetHorizontalVelocity(FEFVSPM s,int *m,int *n,double **data);
ErrorCode FEFVSPMGetKappa(FEFVSPM s,int *m,int *n,double **data);
ErrorCode FEFVSPMGetHeight(FEFVSPM s,int *m,int *n,double **data);

ErrorCode FEFVSPMFieldView_ascii(FEFVSPM s);
ErrorCode FEFVSPMFieldView_VTRB(FEFVSPM s,const char filename[]);
ErrorCode FEFVSPMFieldView_VTR(FEFVSPM s,const char filename[]);
ErrorCode FEFVSPMView_Binary(FEFVSPM s,const char name[]);
ErrorCode FEFVSPMCheckpoint(FEFVSPM s);
ErrorCode FEFVSPMFieldView(FEFVSPM s,FEFVSPMFieldViewType type);
ErrorCode FEFVSPMLoadFromCheckpoint(FEFVSPM s,const char name[]);


ErrorCode FEFVSPMComputeTimestep(FEFVSPM s,double h[],double *dt);

ErrorCode FEFVSPMView(FEFVSPM s);
ErrorCode FEFVSPMStepInitialize(FEFVSPM s);

ErrorCode FEFVSPMStepAdvDiff(FEFVSPM s);

ErrorCode FEFVSPMNodeDataGetRange(FEFVSPM s,double x[],int start,int stride,double range[]);
ErrorCode FEFVSPMCellDataGetRange(FEFVSPM s,double x[],int start,int stride,double range[]);
ErrorCode FEFVSPMNodeDataSet(FEFVSPM s,double x[],double val);
ErrorCode FEFVSPMCellDataSet(FEFVSPM s,double x[],double val);
ErrorCode FEFVSPMNodeDataScale(FEFVSPM s,double x[],double val);
ErrorCode FEFVSPMCellDataScale(FEFVSPM s,double x[],double val);
ErrorCode FEFVSPMNodeDataShift(FEFVSPM s,double x[],double val);
ErrorCode FEFVSPMCellDataShift(FEFVSPM s,double x[],double val);

ErrorCode FEFVSPMCreateNodalField(FEFVSPM s,int ncomponents,double *nodadata[]);
ErrorCode FEFVSPMCreateCellField(FEFVSPM s,int ncomponents,double *celldata[]);
ErrorCode FEFVSPMDestroyNodalField(FEFVSPM s,double **data);
ErrorCode FEFVSPMDestroyCellField(FEFVSPM s,double **data);
ErrorCode FEFVSPMCellGetNodeValues(FEFVSPM s,int ci,int cj,int ncomponents,double nodalfield[],double vertexvalues[]);
ErrorCode FEFVSPMProjectNodalField(FEFVSPM s,double nodedata[]);

/* ============================================================
 * Second-order MUSCL-Hancock finite-volume advection.
 *
 * Solves the non-conservative equation
 *
 *     dh/dt + u_h . grad(h) = 0
 *
 * over [t, t+dt].  The velocity may vary in space but is frozen in
 * time during this substep.  The implementation uses face fluxes
 * together with the matching source h div(u_h), so div(u_h) need not
 * vanish.
 *
 * Median-dual control volumes and face lengths at boundaries.
 * Explicit boundary face fluxes on all four sides.
 * ============================================================ */
ErrorCode FEFVSPMUpdate_FVMAdvection_2ndOrder(FEFVSPM s, double dt, double hc[]);
ErrorCode FEFVSPMGetStableTimestep_FVMAdvection(FEFVSPM s, double cfl, double *dt);

#endif