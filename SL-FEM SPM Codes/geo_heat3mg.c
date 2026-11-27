/* Refactored heat module; extracted from the documented candidate. */
#ifdef MAINFILE
#undef MAINFILE
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#ifndef FEMH_MODEL_HEADER
#define FEMH_MODEL_HEADER "geo_head3mg.h"
#endif
#include FEMH_MODEL_HEADER
#ifndef FEMH_PCG_RTOL
#define FEMH_PCG_RTOL 1.e-10
#endif
#ifndef FEMH_PCG_ATOL
#define FEMH_PCG_ATOL 1.e-8
#endif
#ifndef FEMH_PCG_MAXITER
#define FEMH_PCG_MAXITER 2500
#endif

/* Non-owning aliases: the host retains allocation and ownership. */
typedef struct {
    const float *density, *heat_capacity, *conductivity, *heating, *expansion;
    const double *pressure, *coordinate[3];
    const float *stress[6], *strain[6], *marker_position[3];
    float *temperature, *previous, *source, *outer_temperature, *marker_temperature;
    double *weight_sum, *correction_sum;
    int (*marker_type)(long);
    long node_count, marker_count, shape[3];
    double extent[3];
} HeatHost;

/* Production markt[] is char; the test header uses int. Read values without
 * casting either storage to an incompatible pointer type. */

typedef struct {
    HeatHost host;
    long n,offset[27],lift_count;
    long *lift_rows;
    int exterior_fixed;
    unsigned long zero_lifts,local_lifts,full_lifts;
    unsigned char *fixed;
    double *storage;
    double *mass, *load, *stiffness;
    double *u, *stage, *rhs, *residual, *scaled, *direction;
    double *product, *boundary, *diagonal;
} HeatCall;

static double hr_boundary_time = -DBL_MAX;

/* PCG is expressed against a linear operator, not the Q1 matrix layout. */
typedef struct {
    const HeatCall *assembly;
    double coefficient;
    void (*apply)(const HeatCall *,double,const double *,double *);
    double (*apply_dot)(const HeatCall *,double,const double *,double *);
} HeatOperator;

typedef struct {
    long count;
    const unsigned char *constrained;
    const double *rhs, *diagonal;
    double *answer, *residual, *scaled, *direction, *image;
    double rho, error_squared, tolerance_squared;
} HeatIteration;

typedef struct {
    long center,corner[8];
    double fraction[3],derivative_scale[3],velocity[100];
} HeatSample;

typedef double (*HeatMechanism)(const HeatHost *,const HeatSample *);

typedef struct {
    double old_temperature,new_temperature,length[3],weight[8];
    long cell[3],node[8];
} MarkerSample;

typedef struct {
    double shape[8][8],diffusion[8][36][3];
    unsigned char first[36],second[36],slot[36];
} HeatQuadrature;

/* Internal component contracts; implementations follow the public entry points. */
static void hr_transfer_markers(const HeatHost *host);
static void hr_marker_forward(const HeatHost *host,long marker);
static MarkerSample hr_marker_sample(const HeatHost *host,long marker);
static double hr_marker_diffusion(const HeatHost *host,long marker,const MarkerSample *sample);
static void hr_marker_scatter(const HeatHost *host,const MarkerSample *sample,double change);
static int hr_marker_inside(const HeatHost *host,long marker);
static double hr_freeze_sources(const HeatCall *c);
static double hr_source(const HeatHost *host,long x,long y,long z);
static HeatSample hr_sample(const HeatHost *host,long x,long y,long z);
static double hr_gravity_heating(const HeatHost *host,const HeatSample *sample);
static double hr_pressure_heating(const HeatHost *host,const HeatSample *sample);
static double hr_shear_heating(const HeatHost *host,const HeatSample *sample);
static double hr_corner_weight(const HeatSample *sample,int vertex,int omitted);
static void hr_advance(HeatCall *c, double h, double start);
static int hr_stage(HeatCall *c, double a, double time, double *guess);
static int hr_solve(HeatCall *c, double a, double *x);
static void hr_iteration_step(HeatIteration *it,const HeatOperator *op);
static void hr_iteration_refresh(HeatIteration *it,const HeatOperator *op);
static void hr_iteration_direction(HeatIteration *it,double beta);
static void hr_apply(const HeatCall *c, double a, const double *x, double *out);
static double hr_apply_dot(const HeatCall *,double,const double *,double *);
static double hr_product(const HeatCall *,double,const double *,double *,int);
static inline double hr_row_value(const HeatCall *,double,const double *,long,long,long,long);
static void hr_lift_shell(HeatCall *);
static void hr_reference(HeatQuadrature *);
static inline void hr_assemble_cell(HeatCall *,const HeatQuadrature *,long,long,long);
static void hr_assemble(HeatCall *c);
static void hr_shape(int point, const double length[3],
                     double shape[8], double gradient[8][3]);
static HeatCall hr_open(long count);
static void hr_close(HeatCall *c);
static HeatHost hr_bind_host(void);
static void hr_geometry(void);
static void hr_prepare(HeatCall *c);
static double hr_prescribed(long node, double time);
static int hr_is_fixed(long node);
static long hr_id(long x, long y, long z);
static int hr_slot(int dx, int dy, int dz);
static void *hr_storage(size_t count, size_t width);
static void hr_stop(const char *message);
static int hr_host_marker_type(long marker);
void titerate(int cycle)
{
    HeatCall c;
    HeatHost host=hr_bind_host();
    double duration=timestep, elapsed=0, limit, rate;
    long i;
    int invalid=0;
    if (duration==0) return;
    if (!(duration>0) || !isfinite(duration) || !isfinite(timesum) ||
        !isfinite(timesum+duration)) hr_stop("invalid thermal time interval");
    if (heatdif<0 || !isfinite(heatdif) || marknum<0)
        hr_stop("invalid marker configuration");
#ifdef MAXMRK
    if (marknum>MAXMRK) hr_stop("marker count exceeds host storage");
#endif
    hr_geometry(); c=hr_open(nodenum); hr_prepare(&c);
    #pragma omp parallel for schedule(static)
    for (i=0; i<c.n; ++i) {
        c.u[i]=c.fixed[i] ? hr_prescribed(i,timesum) : host.temperature[i];
        host.previous[i]=host.outer_temperature[i]=(float)c.u[i];
    }
    rate=hr_freeze_sources(&c); limit=duration;
    if (rate>0) {
        if (!(maxtkstep>0) || !isfinite(maxtkstep)) hr_stop("invalid temperature step limit");
        if (maxtkstep/rate<limit) limit=maxtkstep/rate;
    }
    hr_assemble(&c);
    while (elapsed<duration) {
        double remaining=duration-elapsed, h=remaining<limit ? remaining : limit;
        double start=timesum+elapsed;
        if (start<timebeg && timebeg-start<h) h=timebeg-start;
        if (start<timeend && timeend-start<h) h=timeend-start;
        if (!(h>0) || !isfinite(h) || elapsed+h==elapsed)
            hr_stop("thermal substep underflow");
        #pragma omp parallel for schedule(static)
        for (i=0; i<c.n; ++i) host.previous[i]=(float)c.u[i];
        hr_advance(&c,h,start);
        elapsed=h==remaining ? duration : elapsed+h;
    }
    #pragma omp parallel for reduction(|:invalid) schedule(static)
    for (i=0; i<c.n; ++i) {
        if (!isfinite(c.u[i]) || fabs(c.u[i])>FLT_MAX) invalid=1;
    }
    if (invalid) hr_stop("nonfinite/unrepresentable final temperature");
    #pragma omp parallel for schedule(static)
    for (i=0; i<c.n; ++i) host.temperature[i]=(float)c.u[i];
    hr_boundary_time=timesum+duration;
    hr_close(&c);
    hr_transfer_markers(&host);
    if (printmod) printf("heat candidate cycle=%d completed\n",cycle+1);
}

void tkrecalc(void)
{
    HeatHost host=hr_bind_host();
    long i;
    double time=timesum>hr_boundary_time ? timesum : hr_boundary_time;
    hr_geometry();
    for (i=0; i<nodenum; ++i) if (hr_is_fixed(i)) host.temperature[i]=(float)hr_prescribed(i,time);
}

double tbonderr(long node, int mode, double *eps, long *indices, double *values)
{
    (void)node; (void)mode; (void)eps; (void)indices; (void)values;
    hr_stop("legacy boundary row routine is disabled");
    return 0;
}

int gausmatt(int mode, long last, long first)
{
    (void)mode; (void)last; (void)first;
    hr_stop("legacy thermal relaxation is disabled");
    return -1;
}
static void hr_transfer_markers(const HeatHost *host)
{
    long index;
    #pragma omp parallel for schedule(static)
    for(index=0;index<host->node_count;++index) host->weight_sum[index]=host->correction_sum[index]=0;
    #pragma omp parallel for schedule(static)
    for(index=0;index<host->marker_count;++index)
        if(host->marker_type(index)<100 && hr_marker_inside(host,index)) hr_marker_forward(host,index);
    if(heatdif) {
        #pragma omp parallel for schedule(static)
        for(index=0;index<host->node_count;++index)
            if(host->weight_sum[index]!=0) host->correction_sum[index]/=host->weight_sum[index];
        #pragma omp parallel for schedule(static)
        for(index=0;index<host->marker_count;++index) {
            if(host->marker_type(index)<50 && hr_marker_inside(host,index)) {
                MarkerSample sample=hr_marker_sample(host,index);
                int neighbor;
                double correction=0;
                for(neighbor=0;neighbor<8;++neighbor)
                    correction-=sample.weight[neighbor]*host->correction_sum[sample.node[neighbor]];
                host->marker_temperature[index]+=(float)correction;
            }
        }
    }
}

static void hr_marker_forward(const HeatHost *host,long marker)
{
    MarkerSample sample=hr_marker_sample(host,marker);
    if(host->marker_type(marker)>=50) host->marker_temperature[marker]=(float)sample.new_temperature;
    else {
        if(host->marker_temperature[marker]<=0) {
            host->marker_temperature[marker]=(float)sample.old_temperature;
            if(host->marker_position[1][marker]>zdeep && host->marker_temperature[marker]<tdeep)
                host->marker_temperature[marker]=tdeep;
        }
        if(heatdif) {
            double correction=hr_marker_diffusion(host,marker,&sample);
            host->marker_temperature[marker]+=(float)correction;
            hr_marker_scatter(host,&sample,correction);
        }
        host->marker_temperature[marker]+=(float)(sample.new_temperature-sample.old_temperature);
    }
}

static MarkerSample hr_marker_sample(const HeatHost *host,long marker)
{
    MarkerSample sample;
    double temperatures[100],weights[100];
    long nodes[100];
    int index;
    tkcalc1(host->marker_position[0][marker],host->marker_position[1][marker],
            host->marker_position[2][marker],temperatures,nodes,weights);
    sample.old_temperature=temperatures[3];sample.new_temperature=temperatures[2];
    for(index=0;index<3;++index){sample.cell[index]=nodes[index];sample.length[index]=weights[index];}
    for(index=0;index<8;++index){sample.node[index]=nodes[3+index];sample.weight[index]=weights[3+index];}
    return sample;
}

static double hr_marker_diffusion(const HeatHost *host,long marker,const MarkerSample *sample)
{
    double conductivity=0,capacity=0,maximum=-DBL_MAX,strength=0;
    int index,boundary=0;
    for(index=0;index<8;++index) {
        long node=sample->node[index];
        conductivity+=sample->weight[index]*host->conductivity[node];
        capacity+=(host->density[node]*host->heat_capacity[node])*sample->weight[index];
        maximum=fmax(maximum,host->temperature[node]);
    }
    for(index=0;index<3;++index) {
        boundary|=sample->cell[index]<=0 || sample->cell[index]>=host->shape[index]-2;
        strength+=1/(sample->length[index]*sample->length[index]);
    }
    if(boundary && host->marker_temperature[marker]+sample->new_temperature-sample->old_temperature>maximum)
        host->marker_temperature[marker]=(float)(maximum-sample->new_temperature+sample->old_temperature);
    strength=fmax(-150,-2*heatdif*conductivity/capacity*timestep*strength);
    return (sample->old_temperature-(double)host->marker_temperature[marker])*(1-exp(strength));
}

static void hr_marker_scatter(const HeatHost *host,const MarkerSample *sample,double change)
{
    int index;
    for(index=0;index<8;++index) {
        long node=sample->node[index];
        double weight=sample->weight[index],contribution=weight*change;
        #pragma omp atomic update
        host->weight_sum[node]+=weight;
        #pragma omp atomic update
        host->correction_sum[node]+=contribution;
    }
}

static int hr_marker_inside(const HeatHost *host,long marker)
{
    int dimension;
    for(dimension=0;dimension<3;++dimension) {
        double coordinate=host->marker_position[dimension][marker];
        if(!(coordinate>0 && coordinate<host->extent[dimension])) return 0;
    }
    return 1;
}
double heaterr(long x, long y, long z, int mode, int level,
               double *eps, long *indices, double *values)
{
    HeatHost host=hr_bind_host();
    long i;
    double source;
    (void)eps; (void)indices; (void)values;
    if (mode!=0 || level!=0) hr_stop("legacy heaterr mode/level is disabled");
    if (x<1 || x>=xnumx-1 || y<1 || y>=ynumy-1 || z<1 || z>=znumz-1)
        hr_stop("physical source requires an internal finest-grid node");
    i=hr_id(x,y,z); source=hr_source(&host,x,y,z);
    host.source[i]=(float)source;
    return source/((double)host.density[i]*host.heat_capacity[i]);
}

static double hr_freeze_sources(const HeatCall *c)
{
    long x,y,z,i;
    double rate=0;
    if(!adiabyn && !frictyn) {
        #pragma omp parallel for reduction(max:rate) schedule(static)
        for(i=0;i<c->n;++i) {
            c->host.source[i]=c->host.heating[i];
            if(!isfinite(c->host.source[i]))hr_stop("nonfinite physical source");
            {double value=c->fixed[i]?0:fabs((double)c->host.source[i]/((double)c->host.density[i]*c->host.heat_capacity[i]));if(value>rate)rate=value;}
        }
        return rate;
    }
    /* Interior first; exterior extension reads only completed interior data. */
    for (z=1; z<znumz-1; ++z) for (x=1; x<xnumx-1; ++x)
    for (y=1; y<ynumy-1; ++y) c->host.source[hr_id(x,y,z)]=(float)hr_source(&c->host,x,y,z);
    #pragma omp parallel for reduction(max:rate) schedule(static)
    for (i=0; i<c->n; ++i) {
        long iy=i%ynumy, ix=(i/ynumy)%xnumx, iz=i/(xnumx*ynumy);
        if (ix==0 || ix==xnumx-1 || iy==0 || iy==ynumy-1 || iz==0 || iz==znumz-1) {
            long near=hr_id(ix==0?1:(ix==xnumx-1?xnumx-2:ix),
                            iy==0?1:(iy==ynumy-1?ynumy-2:iy),
                            iz==0?1:(iz==znumz-1?znumz-2:iz));
            c->host.source[i]=(float)((double)c->host.heating[i]+(double)c->host.source[near]-c->host.heating[near]);
        }
        if (!isfinite(c->host.source[i])) hr_stop("nonfinite physical source");
        {double value=c->fixed[i]?0:fabs((double)c->host.source[i]/((double)c->host.density[i]*c->host.heat_capacity[i]));if(value>rate)rate=value;}
    }
    return rate;
}

static double hr_source(const HeatHost *host,long x,long y,long z)
{
    const HeatMechanism mechanisms[3]={hr_gravity_heating,hr_pressure_heating,hr_shear_heating};
    const int active[3]={adiabyn==1,adiabyn==2,frictyn!=0};
    long node=hr_id(x,y,z);
    double power=host->heating[node];
    int mechanism;
    if((adiabyn || frictyn) && host->density[node]>2000 && host->coordinate[1][y]<ysize-stp100) {
        HeatSample sample=hr_sample(host,x,y,z);
        for(mechanism=0;mechanism<3;++mechanism)
            if(active[mechanism]) power+=mechanisms[mechanism](host,&sample);
    }
    return power;
}

static HeatSample hr_sample(const HeatHost *host,long x,long y,long z)
{
    HeatSample sample={0};
    const long location[3]={x,y,z};
    int dimension,vertex;
    sample.center=hr_id(x,y,z);
    for(dimension=0;dimension<3;++dimension) {
        const double *axis=host->coordinate[dimension];
        long index=location[dimension];
        double span=axis[index+1]-axis[index-1];
        sample.fraction[dimension]=(axis[index]-axis[index-1])/span;
        sample.derivative_scale[dimension]=2/span;
    }
    for(vertex=0;vertex<8;++vertex)
        sample.corner[vertex]=hr_id(x+((vertex>>1)&1),y+(vertex&1),z+((vertex>>2)&1));
    if(adiabyn) {
        long scratch[100];
        vxyzcalc1(host->coordinate[0][x],host->coordinate[1][y],host->coordinate[2][z],sample.velocity,scratch);
    }
    return sample;
}

static double hr_gravity_heating(const HeatHost *host,const HeatSample *sample)
{
    const double gravity[3]={GXKOEF,GYKOEF,GZKOEF};
    double power=0;
    int dimension;
    for(dimension=0;dimension<3;++dimension) power+=gravity[dimension]*sample->velocity[dimension];
    return (double)host->expansion[sample->center]*host->previous[sample->center]
                  *host->density[sample->center]*power;
}

static double hr_pressure_heating(const HeatHost *host,const HeatSample *sample)
{
    const int bit[3]={1,0,2};
    double transport=0;
    int dimension,vertex;
    for(dimension=0;dimension<3;++dimension) {
        double gradient=0;
        for(vertex=0;vertex<8;++vertex)
            gradient+=(((vertex>>bit[dimension])&1)?1:-1)
                        *hr_corner_weight(sample,vertex,dimension)*host->pressure[sample->corner[vertex]];
        transport+=sample->velocity[dimension]*sample->derivative_scale[dimension]*gradient;
    }
    return (double)host->expansion[sample->center]*host->previous[sample->center]*transport;
}

static double hr_shear_heating(const HeatHost *host,const HeatSample *sample)
{
    const long negative[3]={ynumy,1,xnumx*ynumy};
    const int shear_component[3]={5,4,3};
    double power=0;
    int vertex,dimension;
    for(vertex=0;vertex<8;++vertex) {
        long node=sample->corner[vertex];
        float contraction=0;
        for(dimension=0;dimension<3;++dimension)
            contraction+=host->stress[dimension][node]*host->strain[dimension][node];
        power+=hr_corner_weight(sample,vertex,-1)*contraction;
    }
    for(dimension=0;dimension<3;++dimension) {
        int component=shear_component[dimension];
        long previous=sample->center-negative[dimension];
        double lower=host->stress[component][previous]*host->strain[component][previous];
        double upper=host->stress[component][sample->center]*host->strain[component][sample->center];
        power+=2*(lower+(upper-lower)*sample->fraction[dimension]);
    }
    return power;
}

static double hr_corner_weight(const HeatSample *sample,int vertex,int omitted)
{
    const int bit[3]={1,0,2};
    double value=1;
    int dimension;
    for(dimension=0;dimension<3;++dimension) if(dimension!=omitted)
        value*=((vertex>>bit[dimension])&1) ? sample->fraction[dimension] : 1-sample->fraction[dimension];
    return value;
}
static void hr_advance(HeatCall *c, double h, double start)
{
    const double gamma=1.0-1.0/sqrt(2.0), ratio=(1.0-gamma)/gamma;
    const double a=gamma*h;
    long i;
    int invalid=0, first, second;
    #pragma omp parallel for reduction(|:invalid) schedule(static)
    for (i=0; i<c->n; ++i) {
        c->diagonal[i]=c->mass[i]+a*c->stiffness[27*i+13];
        if (!(c->diagonal[i]>0) || !isfinite(c->diagonal[i])) invalid=1;
        c->rhs[i]=c->mass[i]*c->u[i]+a*c->load[i];
        c->stage[i]=c->u[i];
    }
    if (invalid) hr_stop("invalid stage diagonal");
    first=hr_stage(c,a,start+gamma*h,c->stage);
    #pragma omp parallel for schedule(static)
    for (i=0; i<c->n; ++i) {
        double change=c->stage[i]-c->u[i];
        c->rhs[i]=c->mass[i]*(c->u[i]+ratio*change)+a*c->load[i];
        c->u[i]=c->stage[i]+ratio*change;
    }
    second=hr_stage(c,a,start+h,c->u);
    if (printmod) printf("heat candidate dt=%e PCG=%d,%d\n",h,first,second);
}

static int hr_stage(HeatCall *c, double a, double time, double *guess)
{
    long i;
    int iterations,nonzero=0;
    #pragma omp parallel for reduction(|:nonzero) schedule(static)
    for(i=0;i<c->n;++i) {
        double value=c->fixed[i] ? hr_prescribed(i,time) : 0;
        c->boundary[i]=value;
        nonzero|=value!=0;
        if(c->fixed[i]) guess[i]=c->rhs[i]=0;
    }
    if(!nonzero) ++c->zero_lifts;
    else if(c->exterior_fixed) {
        ++c->local_lifts;
        #pragma omp parallel for schedule(static)
        for(i=0;i<c->lift_count;++i) {
            long row=c->lift_rows[i],y=row%ynumy;
            long x=(row/ynumy)%xnumx,z=row/(xnumx*ynumy);
            c->rhs[row]-=hr_row_value(c,a,c->boundary,row,x,y,z);
        }
    } else {
        ++c->full_lifts;
        hr_apply(c,a,c->boundary,c->product);
        #pragma omp parallel for schedule(static)
        for(i=0;i<c->n;++i) if(!c->fixed[i]) c->rhs[i]-=c->product[i];
    }
    iterations=hr_solve(c,a,guess);
    #pragma omp parallel for schedule(static)
    for(i=0;i<c->n;++i) if(c->fixed[i]) guess[i]=c->boundary[i];
    return iterations;
}
static int hr_solve(HeatCall *c, double a, double *x)
{
    HeatOperator op={c,a,hr_apply,hr_apply_dot};
    HeatIteration it={c->n,c->fixed,c->rhs,c->diagonal,x,c->residual,
                      c->scaled,c->direction,c->product,0,0,0};
    long row,freedom=0;
    double rhs_norm=0;
    int completed=0,maximum=cyc2max>0 ? cyc2max : FEMH_PCG_MAXITER;
    if (!(FEMH_PCG_RTOL>0) || !(FEMH_PCG_ATOL>=0)) hr_stop("invalid PCG tolerances");
    #pragma omp parallel for reduction(+:rhs_norm,freedom) schedule(static)
    for(row=0;row<it.count;++row) {
        if(it.constrained[row]) it.answer[row]=0;
        else {
            double scaled_rhs=it.rhs[row]/it.diagonal[row];
            rhs_norm+=scaled_rhs*scaled_rhs; ++freedom;
        }
    }
    it.tolerance_squared=FEMH_PCG_RTOL*FEMH_PCG_RTOL*rhs_norm
                        +FEMH_PCG_ATOL*FEMH_PCG_ATOL*(double)freedom;
    if(!isfinite(it.tolerance_squared)) hr_stop("nonfinite PCG tolerance");
    hr_iteration_refresh(&it,&op);
    hr_iteration_direction(&it,0);
    for(;;) {
        if(it.error_squared<=it.tolerance_squared || (completed && completed%50==0)) {
            hr_iteration_refresh(&it,&op);
            if(it.error_squared<=it.tolerance_squared) return completed;
            hr_iteration_direction(&it,0);
        }
        if(completed==maximum) hr_stop("PCG iteration limit reached");
        hr_iteration_step(&it,&op);
        ++completed;
    }
}

static void hr_iteration_step(HeatIteration *it,const HeatOperator *op)
{
    long row;
    double previous_rho=it->rho,alpha,next_rho=0,squared=0;
    double curvature=op->apply_dot(op->assembly,op->coefficient,it->direction,it->image);
    if(!(curvature>0) || !(previous_rho>0) || !isfinite(curvature))
        hr_stop("PCG nonpositive curvature or breakdown");
    alpha=previous_rho/curvature;
    #pragma omp parallel for reduction(+:next_rho,squared) schedule(static)
    for(row=0;row<it->count;++row) {
        double residual=0,scaled;
        if(it->constrained[row]) it->answer[row]=0;
        else {
            it->answer[row]+=alpha*it->direction[row];
            residual=it->residual[row]-alpha*it->image[row];
        }
        scaled=residual/it->diagonal[row];
        it->residual[row]=residual;it->scaled[row]=scaled;
        next_rho+=residual*scaled;squared+=scaled*scaled;
    }
    if(!isfinite(next_rho) || !isfinite(squared)) hr_stop("nonfinite PCG recurrence");
    it->rho=next_rho;it->error_squared=squared;
    hr_iteration_direction(it,next_rho/previous_rho);
}

static void hr_iteration_refresh(HeatIteration *it,const HeatOperator *op)
{
    long row;
    double rho=0,squared=0;
    op->apply(op->assembly,op->coefficient,it->answer,it->image);
    #pragma omp parallel for reduction(+:rho,squared) schedule(static)
    for(row=0;row<it->count;++row) {
        double residual=it->constrained[row] ? 0 : it->rhs[row]-it->image[row];
        double scaled=residual/it->diagonal[row];
        it->residual[row]=residual;it->scaled[row]=scaled;
        rho+=residual*scaled;squared+=scaled*scaled;
    }
    if(!isfinite(rho) || !isfinite(squared)) hr_stop("nonfinite PCG residual");
    it->rho=rho;it->error_squared=squared;
}

static void hr_iteration_direction(HeatIteration *it,double beta)
{
    long row;
    #pragma omp parallel for schedule(static)
    for(row=0;row<it->count;++row)
        it->direction[row]=beta==0 ? it->scaled[row] : it->scaled[row]+beta*it->direction[row];
}



static void hr_apply(const HeatCall *c,double a,const double *x,double *out)
{
    (void)hr_product(c,a,x,out,0);
}

static double hr_apply_dot(const HeatCall *c,double a,const double *x,double *out)
{
    return hr_product(c,a,x,out,1);
}

static double hr_product(const HeatCall *c,double a,const double *input,double *out,int need_dot)
{
    long x,z;
    double dot=0;
    #pragma omp parallel for collapse(2) reduction(+:dot) schedule(static)
    for(z=0;z<znumz;++z) for(x=0;x<xnumx;++x) {
        long y,base=(z*xnumx+x)*ynumy;
        for(y=0;y<ynumy;++y) {
            long row=base+y;
            double value=hr_row_value(c,a,input,row,x,y,z);
            out[row]=value;
            if(need_dot) dot+=input[row]*value;
        }
    }
    return dot;
}

static inline double hr_row_value(const HeatCall *c,double a,const double *input,
                                  long row,long x,long y,long z)
{
    const double *coefficients=c->stiffness+27*row;
    double value=c->mass[row]*input[row];
    int slot;
    if(x>0 && x<xnumx-1 && y>0 && y<ynumy-1 && z>0 && z<znumz-1) {
        for(slot=0;slot<27;++slot)
            value+=a*coefficients[slot]*input[row+c->offset[slot]];
    } else {
        for(slot=0;slot<27;++slot) {
            long nx=x+(slot/3)%3-1,ny=y+slot%3-1,nz=z+slot/9-1;
            if(nx>=0 && nx<xnumx && ny>=0 && ny<ynumy && nz>=0 && nz<znumz)
                value+=a*coefficients[slot]*input[row+c->offset[slot]];
        }
    }
    return value;
}
static void hr_assemble(HeatCall *c)
{
    HeatQuadrature reference;
    hr_reference(&reference);
    #pragma omp parallel
    {
        int color;
        for(color=0;color<8;++color) {
            long x,y,z;
            #pragma omp for collapse(3) schedule(static)
            for(z=(color>>2)&1;z<znumz-1;z+=2)
            for(x=(color>>1)&1;x<xnumx-1;x+=2)
            for(y=color&1;y<ynumy-1;y+=2)hr_assemble_cell(c,&reference,x,y,z);
            /* omp for's implicit barrier separates the eight cell colors. */
        }
    }
}

static void hr_reference(HeatQuadrature *reference)
{
    const double unit[3]={1,1,1};
    int point,first,second,pair=0,dimension;
    for(first=0;first<8;++first)for(second=first;second<8;++second) {
        int dx=((second>>1)&1)-((first>>1)&1),dy=(second&1)-(first&1);
        int dz=((second>>2)&1)-((first>>2)&1);
        reference->first[pair]=first;reference->second[pair]=second;
        reference->slot[pair]=(unsigned char)hr_slot(dx,dy,dz);
        ++pair;
    }
    for(point=0;point<8;++point) {
        double gradient[8][3];
        hr_shape(point,unit,reference->shape[point],gradient);
        for(pair=0;pair<36;++pair)for(dimension=0;dimension<3;++dimension)
            reference->diffusion[point][pair][dimension]=
                gradient[reference->first[pair]][dimension]*gradient[reference->second[pair]][dimension];
    }
}

static void hr_shape(int point, const double length[3],
                     double shape[8], double gradient[8][3])
{
    const int bit[3] = {1,0,2};
    const double g = 1.0/sqrt(3.0);
    int a, d;
    double position[3];
    for (d=0; d<3; ++d) position[d]=((point>>bit[d])&1) ? g : -g;
    for (a=0; a<8; ++a) {
        double sign[3], factor[3];
        for (d=0; d<3; ++d) {
            sign[d]=((a>>bit[d])&1) ? 1.0 : -1.0;
            factor[d]=1.0+sign[d]*position[d];
        }
        shape[a]=factor[0]*factor[1]*factor[2]/8.0;
        for (d=0; d<3; ++d)
            gradient[a][d]=sign[d]*factor[(d+1)%3]*factor[(d+2)%3]
                            /(4.0*length[d]);
    }
}



static inline void hr_assemble_cell(HeatCall *c,const HeatQuadrature *reference,long x,long y,long z)
{
    long node[8];
    double local_mass[8]={0},local_load[8]={0},local_k[36]={0};
    double length[3]={c->host.coordinate[0][x+1]-c->host.coordinate[0][x],
                      c->host.coordinate[1][y+1]-c->host.coordinate[1][y],
                      c->host.coordinate[2][z+1]-c->host.coordinate[2][z]};
    double weight=length[0]*length[1]*length[2]/8,metric[3];
    int vertex,point,pair,dimension;
    for(dimension=0;dimension<3;++dimension)metric[dimension]=1/(length[dimension]*length[dimension]);
    for(vertex=0;vertex<8;++vertex)
        node[vertex]=hr_id(x+((vertex>>1)&1),y+(vertex&1),z+((vertex>>2)&1));
    for(point=0;point<8;++point) {
        const double *shape=reference->shape[point];
        double capacity=0,conductivity=0,source=0,scale[3];
        for(vertex=0;vertex<8;++vertex) {
            long row=node[vertex];
            capacity+=shape[vertex]*((double)c->host.density[row]*c->host.heat_capacity[row]);
            conductivity+=shape[vertex]*c->host.conductivity[row];
            source+=shape[vertex]*c->host.source[row];
        }
        for(vertex=0;vertex<8;++vertex) {
            local_mass[vertex]+=weight*capacity*shape[vertex];
            local_load[vertex]+=weight*source*shape[vertex];
        }
        for(dimension=0;dimension<3;++dimension)scale[dimension]=weight*conductivity*metric[dimension];
        for(pair=0;pair<36;++pair) {
            const double *product=reference->diffusion[point][pair];
            local_k[pair]+=product[0]*scale[0]+product[1]*scale[1]+product[2]*scale[2];
        }
    }
    for(vertex=0;vertex<8;++vertex) {
        c->mass[node[vertex]]+=local_mass[vertex];
        c->load[node[vertex]]+=local_load[vertex];
    }
    for(pair=0;pair<36;++pair) {
        int first=reference->first[pair],second=reference->second[pair];
        c->stiffness[27*node[first]+reference->slot[pair]]+=local_k[pair];
        if(first!=second)c->stiffness[27*node[second]+(26-reference->slot[pair])]+=local_k[pair];
    }
}



static HeatCall hr_open(long count)
{
    HeatCall c = {0};
    double **vectors[]={&c.mass,&c.load,&c.u,&c.stage,&c.rhs,&c.residual,
                        &c.scaled,&c.direction,&c.product,&c.boundary,&c.diagonal};
    size_t segment;
    c.n=count; c.host=hr_bind_host();
    for(segment=0;segment<27;++segment) {
        long slot=(long)segment;
        c.offset[segment]=(slot/9-1)*xnumx*ynumy+((slot/3)%3-1)*ynumy+slot%3-1;
    }
    c.fixed=hr_storage((size_t)count,sizeof(*c.fixed));
    c.storage=hr_storage((size_t)count,38*sizeof(double));
    c.stiffness=c.storage;
    for(segment=0;segment<sizeof(vectors)/sizeof(vectors[0]);++segment)
        *vectors[segment]=c.storage+(27+segment)*(size_t)count;
    return c;
}

static void hr_close(HeatCall *c)
{
    free(c->fixed); free(c->storage); free(c->lift_rows);
    memset(c,0,sizeof(*c));
}

static HeatHost hr_bind_host(void)
{
    HeatHost host = {
        .density=ro, .heat_capacity=cp, .conductivity=kt,
        .heating=ht, .expansion=et, .pressure=pr,
        .coordinate={gx,gy,gz},
        .stress={sxx,syy,szz,sxy,sxz,syz},
        .strain={exx,eyy,ezz,exy,exz,eyz},
        .marker_position={markx,marky,markz},
        .temperature=tk, .previous=tk0, .source=tk1, .outer_temperature=tk2,
        .marker_temperature=markk, .weight_sum=sol0, .correction_sum=sol1,
        .marker_type=hr_host_marker_type,
        .node_count=nodenum, .marker_count=marknum,
        .shape={xnumx,ynumy,znumz}, .extent={xsize,ysize,zsize}
    };
    return host;
}

static void hr_geometry(void)
{
    long j;
    if (xnumx < 3 || ynumy < 3 || znumz < 3 ||
        ynumy > MAXNNN || xnumx > MAXNNN/ynumy ||
        znumz > MAXNNN/(xnumx*ynumy)) hr_stop("invalid finest grid size");
    if (nodenum != xnumx*ynumy*znumz || nodenum > MAXPPP ||
        nodenum8 < 0 || nodenum > MAXPAR || nodenum8 > MAXPAR-nodenum)
        hr_stop("inconsistent host grid storage");
    if (!isfinite(gx[0]) || !isfinite(gy[0]) || !isfinite(gz[0]))
        hr_stop("nonfinite coordinate origin");
    for (j=1; j<xnumx; ++j)
        if (!isfinite(gx[j]) || gx[j]<=gx[j-1]) hr_stop("invalid x coordinates");
    for (j=1; j<ynumy; ++j)
        if (!isfinite(gy[j]) || gy[j]<=gy[j-1]) hr_stop("invalid y coordinates");
    for (j=1; j<znumz; ++j)
        if (!isfinite(gz[j]) || gz[j]<=gz[j-1]) hr_stop("invalid z coordinates");
}

static void hr_prepare(HeatCall *c)
{
    long j;
    int invalid = 0;
    #pragma omp parallel for reduction(|:invalid) schedule(static)
    for (j=0; j<c->n; ++j) {
        if (!(c->host.density[j]>0 && c->host.heat_capacity[j]>0 && c->host.conductivity[j]>0) ||
            !isfinite(c->host.density[j]) || !isfinite(c->host.heat_capacity[j]) || !isfinite(c->host.conductivity[j]) ||
            !isfinite(c->host.heating[j]) || !isfinite(c->host.temperature[j])) invalid = 1;
    }
    if (invalid) hr_stop("nonpositive/nonfinite material or temperature/source");
    /* Boundary topology is frozen for this thermal call. */
    free(c->lift_rows);c->lift_rows=NULL;c->lift_count=0;
    c->exterior_fixed=1;
    {
        int any_fixed=0;
        for(j=0;j<c->n;++j) {
            c->fixed[j]=(unsigned char)hr_is_fixed(j);
            if(c->fixed[j]) {
                long y=j%ynumy,x=(j/ynumy)%xnumx,z=j/(xnumx*ynumy);
                any_fixed=1;
                if(x>0 && x<xnumx-1 && y>0 && y<ynumy-1 && z>0 && z<znumz-1)
                    c->exterior_fixed=0;
            }
        }
        if(any_fixed && c->exterior_fixed) {
            long inner[3]={xnumx>4?xnumx-4:0,ynumy>4?ynumy-4:0,znumz>4?znumz-4:0};
            long capacity=c->n-inner[0]*inner[1]*inner[2];
            c->lift_rows=hr_storage((size_t)capacity,sizeof(*c->lift_rows));
            hr_lift_shell(c);
        }
    }
}

static double hr_prescribed(long node, double time)
{
    double value = bondv1[bondm[nodenum8 + node]][0];
    if (node % ynumy == ynumy-1 && time > timebeg) {
        if (!(timeend > timebeg)) hr_stop("invalid boundary ramp interval");
        value = time >= timeend ? tempend :
            value + (tempend-value) * ((time-timebeg)/(timeend-timebeg));
    }
    if (!isfinite(value)) hr_stop("nonfinite prescribed temperature");
    return value;
}

static int hr_is_fixed(long node)
{
    long tag = bondm[nodenum8 + node], other;
    long x, y, z, xx, yy, zz;
    double constant, coefficient;
    if (!tag) return 0;
    if (tag < 0 || tag >= MAXBON) hr_stop("invalid boundary index");
    constant = bondv1[tag][0]; coefficient = bondv1[tag][1];
    if (!isfinite(constant) || !isfinite(coefficient))
        hr_stop("nonfinite boundary condition");
    if (!bondn1[tag] || coefficient == 0) return 1;
    other = (long)bondn1[tag] - 1;
    if (other < 0 || other >= nodenum || constant != 0 || coefficient != 1)
        hr_stop("unsupported thermal affine constraint");
    y = node % ynumy; x = (node / ynumy) % xnumx;
    z = node / (xnumx * ynumy);
    yy = other % ynumy; xx = (other / ynumy) % xnumx;
    zz = other / (xnumx * ynumy);
    if (y == yy && z == zz &&
        ((x == 0 && xx == 1) || (x == xnumx-1 && xx == xnumx-2))) return 0;
    if (x == xx && z == zz &&
        ((y == 0 && yy == 1) || (y == ynumy-1 && yy == ynumy-2))) return 0;
    if (x == xx && y == yy &&
        ((z == 0 && zz == 1) || (z == znumz-1 && zz == znumz-2))) return 0;
    hr_stop("zero-gradient relation is not an exterior normal neighbor");
    return 0;
}

static long hr_id(long x, long y, long z)
{
    return (z * xnumx + x) * ynumy + y;
}

static int hr_slot(int dx, int dy, int dz)
{
    return 9 * (dz + 1) + 3 * (dx + 1) + dy + 1;
}

static void *hr_storage(size_t count, size_t width)
{
    void *memory;
    if (!count || !width || count > SIZE_MAX / width)
        hr_stop("invalid workspace size");
    memory = calloc(count, width);
    if (!memory) hr_stop("workspace allocation failed");
    return memory;
}

static void hr_stop(const char *message)
{
    fprintf(stderr, "heat candidate: %s\n", message);
    exit(EXIT_FAILURE);
}

static int hr_host_marker_type(long marker)
{
    return markt[marker];
}

static void hr_lift_shell(HeatCall *c)
{
    long x,z;
    for(z=0;z<znumz;++z)for(x=0;x<xnumx;++x) {
        int full=z<2 || z>=znumz-2 || x<2 || x>=xnumx-2;
        int part;
        for(part=0;part<(full?1:2);++part) {
            long begin=part?(ynumy-2>2?ynumy-2:2):0;
            long end=part?ynumy:(full?ynumy:2),y;
            for(y=begin;y<end;++y) {
                long row=hr_id(x,y,z);
                if(!c->fixed[row])c->lift_rows[c->lift_count++]=row;
            }
        }
    }
}