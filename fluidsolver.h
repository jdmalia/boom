#ifndef FLUIDSOLVER_H
#define FLUIDSOLVER_H

#include <QObject>
#include "Eigen"
#include "Sparse"

using namespace Eigen;

#define DIM 100 // columns & rows
#define SIZE ((DIM+2)*(DIM+2))
#define IX(i,j) ((i)*(DIM+2)+(j))
#define SWAP(x0, x) {float* tmp=x0; x0=x; x=tmp;}

#define MAX_PARTICLES 100000
#define CG_ITERATIONS 40
#define CG_TOLERANCE 1e-8
#define CFL_LIMIT 1

enum Numeric_Method {CONJ_GRAD, GAUSS_SEIDEL};
enum Velocity_Method {FORWARD_EULER, MIDPOINT};
enum Performance {LOW, NORMAL, HIGH};
enum Blast_Shape {CIRCLE, DIAMOND, SQUARE};

class FluidSolver : QObject
{
    Q_OBJECT

public:
    FluidSolver();
    ~FluidSolver();

    void step();
    float max_vel();
    void vel_field(int);
    void scalar_field(int);

    float* u;       // horz component of velocity
    float* u_prev;
    float* u_bfecc_1;
    float* u_bfecc_2;
    float* u_acc;
    float* v;       // vert component of velocity
    float* v_prev;
    float* v_bfecc_1;
    float* v_bfecc_2;
    float* v_acc;
    float* d;       // densities
    float* d_prev;
    float* d_bfecc_1;
    float* d_bfecc_2;
    float* scalars;
    float* t;       // temperatures
    float* t_prev;
    float* t_bfecc_1;
    float* t_bfecc_2;
    float* divergence;

    float total_divergence;

    bool isDone;

    float d0; //threshold density

    float t_amb; //ambient temp
    float t_max;
    float cr;

    float blast_radius;
    float blast_intensity;
    Blast_Shape blast_shape;

    float* div_constants;
    float div_falloff;

    int num_particles;
    Performance performance;
    float *particle_pos_X, *particle_pos_Y, *particle_temp;
    float ah, r, cm;

    float dt;
    float h_space;

    bool buoyancy_enabled;
    bool bad_diffusion;
    bool bfecc;
    bool simEnabled;
    bool oneStep;

    bool temp_diff_enable;
    bool temp_adv_enable;
    bool temp_env_enable;

    bool multi_threading;

    Numeric_Method method;
    Velocity_Method vel_method;

    int cfl_violations;

    bool sum_divergence;
    int iterations;

public slots:
    void set_dt(float);
    void set_diff(float);
    void set_visc(float);
    void set_blast_rad(float);
    void set_blast_intensity(float);
    void diffusion_test();
    void advection_test();
    void buoyancy_test();
    void cig_test();
    void reset();
    void buoyancyEnable(bool);
    void set_buoy(float);

private:

    float diff;
    float heat_diff;
    float visc;
    float buoy;

    SparseMatrix<float> conj_grad_A_diff;
    SparseMatrix<float> conj_grad_A_visc;
    SparseMatrix<float> conj_grad_A_temp;
    SparseMatrix<float> conj_grad_A_proj;
    VectorXf conj_grad_x;
    VectorXf conj_grad_b;

    float cg_tolerance;

    void allocateMemory();
    void deallocateMemory();
    void initializeConjGradMat(SparseMatrix<float> * conj_grad_A, float k, float dt, float h);
    void initializeConjGradProj();

    void arrToVec(VectorXf* vec, float* arr, int size);
    void vecToArr(float* arr, VectorXf* vec, int size);

    void add_source(int n, float * x, float * s, float dt);
    void diffuse(int n, int b, float * x, float * x0, float diff, float dt);
    void diffuse_cg(int n, int b, float * x, float * x0, SparseMatrix<float> * conj_grad_A);
    void diffuse_bad(int n, int b, float * x, float * x0, float diff, float dt);
    void advect(int n, int b, float * d, float * d0, float * u, float * v, float dt);
    void advect_back(int n, int b, float * d, float * d0, float * u, float * v, float dt);
    void advect_bfecc(int n, int b, float * d, float * d0, float * d_bfecc_1, float* d_bfecc_2, float * u, float * v, float dt);
    void project ( int n, float * u, float * v, float * p, float * div );
    void project_cg( int n, float * u, float * v, float * p, float * div );
    void dens_step(int n, float * x, float * x0, float * u, float * v, float diff, float dt);
    void vel_step ( int n, float * u, float * v, float * u0, float * v0, float visc, float dt );
    void temp_step();
    void environment_loss(int n, float * x, float * x0, float x_amb, float x_max, float cr, float dt);
    void external_forces();
    void set_bnd(int n, int b, float * x);
    float myAbs(float);
    void cfl_check();

    void update_particles();
    void update_particles_MT();
    static void *update_particles_ST(void* arg);
};

#endif // FLUIDSOLVER_H
