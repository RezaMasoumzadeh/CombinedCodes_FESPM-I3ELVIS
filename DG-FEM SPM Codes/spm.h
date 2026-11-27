#ifndef __spm_h__
#define __spm_h__

#include <stdio.h>
#include <stdlib.h>

typedef int ErrorCode;

typedef enum {  DM_CONST, DM_NONLINEAR_CRITICALSLOPE, DM_NONLINEAR_POWERLAW,
                DM_UNDEFINED } FEDGSPMDiffusionMode;
typedef enum { DTM_CONST, DTM_EXPLICIT } FEDGSPMTimestepMode;

typedef enum { DIR_JMAX=0, DIR_IMAX, DIR_JMIN, DIR_IMIN } FEDGSPMDirection;

typedef enum { FVT_ASCII, FVT_VTR, FVT_VTRB, FVT_UNDEFINED } FEDGSPMFieldViewType;

#define FEDGSPM_IGNORE -1
#define FEDGSPM_DECIDE -2

typedef struct _p_FEDGSPM *FEDGSPM;

typedef struct _p_FESPM *FESPM;

struct _p_FEDGSPM {
    double xrange[2];
    double yrange[2];
    double time_end,dt,dt_const,time;
    int nx,ny,mx,my,nnodes,ncells;
    double dx,dy,dt_theta;
    int nt,nt_io,nt_max;
    FEDGSPMDiffusionMode mode;
    FEDGSPMTimestepMode dtmode;
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

#define FEDGSPMCheckIsIgnore(a) ( (int)((a)) == FEDGSPM_IGNORE ? 1 : 0 )

#include "spmgriditerator.h"

#define ERRQ(str) { printf("*** FEDGSPM ERROR: -------------------------- Error Message --------------------------\n"\
                           "*** FEDGSPM ERROR: %s\n"\
                           "*** FEDGSPM ERROR: %s() line %d in %s\n"\
                           "*** FEDGSPM ERROR: -------------------------- End of Error Message --------------------------\n",\
                           (str),__FUNCTION__,__LINE__,__FILE__); exit(EXIT_FAILURE); }

#define ERRQ2(str1,str2) { printf("*** FEDGSPM ERROR: -------------------------- Error Message --------------------------\n"\
                                  "*** FEDGSPM ERROR: %s %s\n"\
                                  "*** FEDGSPM ERROR: %s() line %d in %s\n"\
                                  "*** FEDGSPM ERROR: -------------------------- End of Error Message --------------------------\n",\
                                  (str1),(str2),__FUNCTION__,__LINE__,__FILE__); exit(EXIT_FAILURE); }


ErrorCode FEDGSPMCreate(FEDGSPM *s);
ErrorCode FEDGSPMDestroy(FEDGSPM *s);

ErrorCode FEDGSPMSetGeometry(FEDGSPM s,double x0,double x1,double y0,double y1);
ErrorCode FEDGSPMSetResolution(FEDGSPM s,int nx,int ny);

ErrorCode FEDGSPMSetTimestepBounds(FEDGSPM spm,double dt_min,double dt_max);
ErrorCode FEDGSPMSetTimestep(FEDGSPM s,double dt);
ErrorCode FEDGSPMSetTimestepTheta(FEDGSPM s,double theta);
ErrorCode FEDGSPMSetMaxNumberSteps(FEDGSPM s,int nt_max);
ErrorCode FEDGSPMSetEndTime(FEDGSPM s,double tend);
ErrorCode FEDGSPMSetTimestepMode(FEDGSPM s,FEDGSPMTimestepMode mode);
ErrorCode FEDGSPMSetTimestepIO(FEDGSPM s,int nt);
ErrorCode FEDGSPMSetConstantTimestep(FEDGSPM s,double dt);
ErrorCode FEDGSPMSetExplicitTimestep(FEDGSPM s,double theta);

ErrorCode FEDGSPMSetDiffusionMode(FEDGSPM s,FEDGSPMDiffusionMode mode);
ErrorCode FEDGSPMSetConstantDiffusivity(FEDGSPM s,double k0);
ErrorCode FEDGSPMSetNonlinearDiffusivityCriticalSlope(FEDGSPM s,double k0,double Sc,double d_exp);
ErrorCode FEDGSPMSetNonlinearDiffusivityPowerLaw(FEDGSPM s,double k0,double d_exp);

ErrorCode FEDGSPMSetInitialHeight(FEDGSPM s,double h0[]);
ErrorCode FEDGSPMSetVerticalVelocity(FEDGSPM s,double v[]);
ErrorCode FEDGSPMSetHorizontalVelocity(FEDGSPM s,double uw[]);

ErrorCode FEDGSPMSetBoundaryConditions(FEDGSPM s,const char bctype[]);

ErrorCode FEDGSPMGetNodeCoordinate(FEDGSPM s,int i,int j,double coor[]);
ErrorCode FEDGSPMGetCentroidCoordinate(FEDGSPM s,int i,int j,double coor[]);
ErrorCode FEDGSPMGetVerticalVelocity(FEDGSPM s,int *m,int *n,double **data);
ErrorCode FEDGSPMGetHorizontalVelocity(FEDGSPM s,int *m,int *n,double **data);
ErrorCode FEDGSPMGetKappa(FEDGSPM s,int *m,int *n,double **data);
ErrorCode FEDGSPMGetHeight(FEDGSPM s,int *m,int *n,double **data);

ErrorCode FEDGSPMFieldView_ascii(FEDGSPM s);
ErrorCode FEDGSPMFieldView_VTRB(FEDGSPM s,const char filename[]);
ErrorCode FEDGSPMFieldView_VTR(FEDGSPM s,const char filename[]);
ErrorCode FEDGSPMView_Binary(FEDGSPM s,const char name[]);
ErrorCode FEDGSPMCheckpoint(FEDGSPM s);
ErrorCode FEDGSPMFieldView(FEDGSPM s,FEDGSPMFieldViewType type);
ErrorCode FEDGSPMLoadFromCheckpoint(FEDGSPM s,const char name[]);


ErrorCode FEDGSPMComputeTimestep(FEDGSPM s,double h[],double *dt);
ErrorCode FEDGSPMStepInitialize(FEDGSPM s);

ErrorCode FEDGSPMView(FEDGSPM s);

ErrorCode FEDGSPMStepAdvDiff(FEDGSPM s);
ErrorCode FEDGSPMStepAdvDiff_TimeDependentInput(FEDGSPM s,ErrorCode (*user_update_function)(FEDGSPM,double,int,double,void*),void* ctx);

ErrorCode FEDGSPMNodeDataGetRange(FEDGSPM s,double x[],int start,int stride,double range[]);
ErrorCode FEDGSPMCellDataGetRange(FEDGSPM s,double x[],int start,int stride,double range[]);
ErrorCode FEDGSPMNodeDataSet(FEDGSPM s,double x[],double val);
ErrorCode FEDGSPMCellDataSet(FEDGSPM s,double x[],double val);
ErrorCode FEDGSPMNodeDataScale(FEDGSPM s,double x[],double val);
ErrorCode FEDGSPMCellDataScale(FEDGSPM s,double x[],double val);
ErrorCode FEDGSPMNodeDataShift(FEDGSPM s,double x[],double val);
ErrorCode FEDGSPMCellDataShift(FEDGSPM s,double x[],double val);

ErrorCode FEDGSPMCreateNodalField(FEDGSPM s,int ncomponents,double *nodadata[]);
ErrorCode FEDGSPMCreateCellField(FEDGSPM s,int ncomponents,double *celldata[]);
ErrorCode FEDGSPMDestroyNodalField(FEDGSPM s,double **data);
ErrorCode FEDGSPMDestroyCellField(FEDGSPM s,double **data);
ErrorCode FEDGSPMCellGetNodeValues(FEDGSPM s,int ci,int cj,int ncomponents,double nodalfield[],double vertexvalues[]);
ErrorCode FEDGSPMProjectNodalField(FEDGSPM s,double nodedata[]);

/* P1 DG advection with upwind flux and the velocity-divergence
   correction required by the original non-conservative equation. */
ErrorCode FEDGSPMUpdate_DGAdvection_P1_SSPRK2(FEDGSPM s, double dt, double hc[]);

#endif
