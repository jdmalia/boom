#include <cstdlib>
#include <cstdio>

#include "fluidsolver.h"
#include "Dense"

#include <iostream>
#include <fstream>

#define ENERGY_FILE "/Users/jason/Desktop/cs4699_data/energies.csv"
#define DIVERGENCE_FILE "/Users/jason/Desktop/cs4699_data/divergence_sums%02d.csv"

#define DIV_ITERATIONS 500

#include "logging.h"

#include <pthread.h>
#define NUM_THREADS 1

struct thread_data {
    int tid;
    FluidSolver* solver;
};

using namespace std;

bool sum_velocities = false;
float energy[1000];
int file_num = 0;

FluidSolver::FluidSolver() {

    isDone = true;

    simEnabled = true;

    allocateMemory();
    dt = 0.1;
    diff = 0.00001;
    heat_diff = 0.00001;
    visc = 0.0001;
    buoy = 0.03;

    d0 = 0.0;

    t_amb = 0.0;
    t_max = 1.0;
    cr = .01;   // completely arbitrary

    temp_diff_enable = true;
    temp_adv_enable = true;
    temp_env_enable = true;

    num_particles = 0;

    bad_diffusion = false;

    h_space = 1.0/(DIM+2);

    ah = 1;
    r = 1;
    cm = 1;
    z = 1;
    Bh = 200;
    temp_threshold = 0.9;
    exploding_particles = true;
    explosiveParticles = true;

    initializeConjGradMat(&conj_grad_A_diff, diff, dt, h_space);
    initializeConjGradMat(&conj_grad_A_temp, diff, dt, h_space);
    initializeConjGradMat(&conj_grad_A_visc, visc, dt, h_space);
    initializeConjGradProj();

    conj_grad_x = VectorXf(DIM*DIM);
    conj_grad_b = VectorXf(DIM*DIM);

    cg_tolerance = 0.1;

    method = GAUSS_SEIDEL;
    bfecc = true;
    oneStep = false;

    performance = NORMAL;

    vel_method = FORWARD_EULER;

    div_falloff = 0.005;

    blast_radius = 8;
    blast_intensity = .01;
    blast_shape = CIRCLE;

    buoyancy_enabled = false;

    multi_threading = false;

    cfl_violations = 0;

    total_divergence = 0;
    sum_divergence = false;
    iterations = 0;
}

FluidSolver::~FluidSolver() {
    deallocateMemory();
}

void FluidSolver::allocateMemory() {

    u = (float *)calloc(SIZE, sizeof(float));
    if(NULL==u) {
        // Handle error
        exit(-1);
    }
    u_prev = (float *)calloc(SIZE, sizeof(float));
    if(NULL==u_prev) {
        // Handle error
        exit(-1);
    }
    u_bfecc_1 = (float *)calloc(SIZE, sizeof(float));
    if(NULL==u_bfecc_1) {
        // Handle error
        exit(-1);
    }
    u_bfecc_2 = (float *)calloc(SIZE, sizeof(float));
    if(NULL==u_bfecc_2) {
        // Handle error
        exit(-1);
    }
    v = (float *)calloc(SIZE, sizeof(float));
    if(NULL==v) {
        // Handle error
        exit(-1);
    }
    v_prev = (float *)calloc(SIZE, sizeof(float));
    if(NULL==v_prev) {
        // Handle error
        exit(-1);
    }
    v_bfecc_1 = (float *)calloc(SIZE, sizeof(float));
    if(NULL==v_bfecc_1) {
        // Handle error
        exit(-1);
    }
    v_bfecc_2 = (float *)calloc(SIZE, sizeof(float));
    if(NULL==v_bfecc_2) {
        // Handle error
        exit(-1);
    }
    t = (float *)calloc(SIZE, sizeof(float));
    if(NULL==t) {
        // Handle error
        exit(-1);
    }
    t_prev = (float *)calloc(SIZE, sizeof(float));
    if(NULL==t_prev) {
        // Handle error
        exit(-1);
    }
    t_bfecc_1 = (float *)calloc(SIZE, sizeof(float));
    if(NULL==t_bfecc_1) {
        // Handle error
        exit(-1);
    }
    t_bfecc_2 = (float *)calloc(SIZE, sizeof(float));
    if(NULL==t_bfecc_2) {
        // Handle error
        exit(-1);
    }
    d = (float *)calloc(SIZE, sizeof(float));
    if(NULL==d) {
        // Handle error
        exit(-1);
    }
    d_prev = (float *)calloc(SIZE, sizeof(float));
    if(NULL==d_prev) {
        // Handle error
        exit(-1);
    }
    d_bfecc_1 = (float *)calloc(SIZE, sizeof(float));
    if(NULL==d_bfecc_1) {
        // Handle error
        exit(-1);
    }
    d_bfecc_2 = (float *)calloc(SIZE, sizeof(float));
    if(NULL==d_bfecc_2) {
        // Handle error
        exit(-1);
    }
    particle_pos_X = (float *)calloc(MAX_PARTICLES, sizeof(float));
    if(NULL==particle_pos_X) {
        // Handle error
        exit(-1);
    }
    particle_pos_Y = (float *)calloc(MAX_PARTICLES, sizeof(float));
    if(NULL==particle_pos_Y) {
        // Handle error
        exit(-1);
    }
    particle_temp = (float *)calloc(MAX_PARTICLES, sizeof(float));
    if(NULL==particle_temp) {
        // Handle error
        exit(-1);
    }
    particle_life = (float *)calloc(MAX_PARTICLES, sizeof(float));
    if(NULL==particle_life) {
        // Handle error
        exit(-1);
    }
    scalars = (float *)calloc(SIZE, sizeof(float));
    if(NULL==scalars) {
        // Handle error
        exit(-1);
    }
    div_constants = (float *)calloc(SIZE, sizeof(float));
    if(NULL==div_constants) {
        // Handle error
        exit(-1);
    }
    u_acc = (float *)calloc(SIZE, sizeof(float));
    if(NULL==u_acc) {
        // Handle error
        exit(-1);
    }
    v_acc = (float *)calloc(SIZE, sizeof(float));
    if(NULL==v_acc) {
        // Handle error
        exit(-1);
    }
    divergence = (float *)calloc(SIZE, sizeof(float));
    if(NULL==divergence) {
        // Handle error
        exit(-1);
    }
}

void FluidSolver::deallocateMemory() {

    if(NULL != u) free(u);
    if(NULL != u_prev) free(u_prev);
    if(NULL != u_bfecc_1) free(u_bfecc_1);
    if(NULL != u_bfecc_2) free(u_bfecc_2);
    if(NULL != v) free(v);
    if(NULL != v_prev) free(v_prev);
    if(NULL != v_bfecc_1) free(v_bfecc_1);
    if(NULL != v_bfecc_2) free(v_bfecc_2);
    if(NULL != d) free(d);
    if(NULL != d_prev) free(d_prev);
    if(NULL != d_bfecc_1) free(d_bfecc_1);
    if(NULL != d_bfecc_2) free(d_bfecc_2);
    if(NULL != t) free(t);
    if(NULL != t_prev) free(t_prev);
    if(NULL != t_bfecc_1) free(t_bfecc_1);
    if(NULL != t_bfecc_2) free(t_bfecc_2);
    if(NULL != particle_pos_X) free(particle_pos_X);
    if(NULL != particle_pos_Y) free(particle_pos_Y);
    if(NULL != particle_temp) free(particle_temp);
    if(NULL != scalars) free(scalars);
    if(NULL != div_constants) free(div_constants);
    if(NULL != u_acc) free(u_acc);
    if(NULL != v_acc) free(v_acc);
}

void FluidSolver::initializeConjGradMat(SparseMatrix<float> *conj_grad_A, float k, float dt, float h){

    float L = (dt*k)/(h*h);

    int diag_pos=0;

    *conj_grad_A = SparseMatrix<float>(DIM*DIM, DIM*DIM);

    for(int i=0; i<DIM; i++){
        for(int j=0; j<DIM; j++){

            if (i==0 && j==0) {
                conj_grad_A->coeffRef(diag_pos, diag_pos) = 1 + 2*L;
                conj_grad_A->coeffRef(diag_pos, diag_pos+1) = -L;
                conj_grad_A->coeffRef(diag_pos, diag_pos+DIM) = -L;
            }
            else if (i==0 && j==DIM-1) {
                conj_grad_A->coeffRef(diag_pos, diag_pos) = 1 + 2*L;
                conj_grad_A->coeffRef(diag_pos, diag_pos-1) = -L;
                conj_grad_A->coeffRef(diag_pos, diag_pos+DIM) = -L;
            }
            else if (i==DIM-1 && j==0) {
                conj_grad_A->coeffRef(diag_pos, diag_pos) = 1 + 2*L;
                conj_grad_A->coeffRef(diag_pos, diag_pos+1) = -L;
                conj_grad_A->coeffRef(diag_pos, diag_pos-DIM) = -L;
            }
            else if (i==DIM-1 && j==DIM-1) {
                conj_grad_A->coeffRef(diag_pos, diag_pos) = 1 + 2*L;
                conj_grad_A->coeffRef(diag_pos, diag_pos-1) = -L;
                conj_grad_A->coeffRef(diag_pos, diag_pos-DIM) = -L;
            }
            else if (i==0) {
                conj_grad_A->coeffRef(diag_pos, diag_pos) = 1 + 3*L;
                conj_grad_A->coeffRef(diag_pos, diag_pos+1) = -L;
                conj_grad_A->coeffRef(diag_pos, diag_pos-1) = -L;
                conj_grad_A->coeffRef(diag_pos, diag_pos+DIM) = -L;
            }
            else if (j==0){
                conj_grad_A->coeffRef(diag_pos, diag_pos) = 1 + 3*L;
                conj_grad_A->coeffRef(diag_pos, diag_pos-DIM) = -L;
                conj_grad_A->coeffRef(diag_pos, diag_pos+1) = -L;
                conj_grad_A->coeffRef(diag_pos, diag_pos+DIM) = -L;
            }
            else if (i==DIM-1){
                conj_grad_A->coeffRef(diag_pos, diag_pos) = 1 + 3*L;
                conj_grad_A->coeffRef(diag_pos, diag_pos-1) = -L;
                conj_grad_A->coeffRef(diag_pos, diag_pos+1) = -L;
                conj_grad_A->coeffRef(diag_pos, diag_pos-DIM) = -L;
            }
            else if (j==DIM-1){
                conj_grad_A->coeffRef(diag_pos, diag_pos) = 1 + 3*L;
                conj_grad_A->coeffRef(diag_pos, diag_pos-1) = -L;
                conj_grad_A->coeffRef(diag_pos, diag_pos+DIM) = -L;
                conj_grad_A->coeffRef(diag_pos, diag_pos-DIM) = -L;
            }
            else {
                conj_grad_A->coeffRef(diag_pos, diag_pos) = 1 + 4*L;
                conj_grad_A->coeffRef(diag_pos, diag_pos-1) = -L;
                conj_grad_A->coeffRef(diag_pos, diag_pos+1) = -L;
                conj_grad_A->coeffRef(diag_pos, diag_pos-DIM) = -L;
                conj_grad_A->coeffRef(diag_pos, diag_pos+DIM) = -L;
            }

            diag_pos++;
        }
    }

}

void FluidSolver::initializeConjGradProj(){

//    conj_grad_A_proj = SparseMatrix<float>(DIM*DIM, DIM*DIM);

//    int diag_pos = 0;

//    float L = 1.0;

//    for(int i=0; i<DIM; i++){
//        for(int j=0; j<DIM; j++){

//            if (i==0 && j==0) {
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 2*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+1) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+DIM) = -1*L;
//            }
//            else if (i==0 && j==DIM-1) {
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 2*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-1) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+DIM) = -1*L;
//            }
//            else if (i==DIM-1 && j==0) {
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 2*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+1) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-DIM) = -1*L;
//            }
//            else if (i==DIM-1 && j==DIM-1) {
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 2*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-1) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-DIM) = -1*L;
//            }
//            else if (i==0) {
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 3*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+1) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-1) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+DIM) = -1*L;
//            }
//            else if (j==0){
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 3*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-DIM) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+1) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+DIM) = -1*L;
//            }
//            else if (i==DIM-1){
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 3*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-1) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+1) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-DIM) = -1*L;
//            }
//            else if (j==DIM-1){
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 3*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-1) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+DIM) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-DIM) = -1*L;
//            }
//            else {
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 4*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-1) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+1) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-DIM) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+DIM) = -1*L;
//            }

//            diag_pos++;
//        }
//    }
    conj_grad_A_proj = SparseMatrix<float>(SIZE, SIZE);

    int diag_pos = 0;

    float L = 1.0;

    for(int i=0; i<DIM+2; i++){
        for(int j=0; j<DIM+2; j++){

            if (i==0 && j==0) {
 //               conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 1;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 2*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+1) = -1*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+(DIM+2)) = -1*L;
            }
//            else if (i==0 && j==1) {
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 3*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+1) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+(DIM+2)) = -1*L;
//            }
//            else if (i==1 && j==0) {
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 3*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+1) = -1*L;
//                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+DIM+2) = -1*L;
//            }
            else if (i==0 && j==DIM+1) {
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 2*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-1) = -1*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+(DIM+2)) = -1*L;
            }
            else if (i==DIM+1 && j==0) {
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 2*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+1) = -1*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-(DIM+2)) = -1*L;
            }
            else if (i==DIM+1 && j==DIM+1) {
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 2*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-1) = -1*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-(DIM+2)) = -1*L;
            }
            else if (i==0) {
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 3*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+1) = -1*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-1) = -1*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+(DIM+2)) = -1*L;
            }
            else if (j==0){
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 3*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-(DIM+2)) = -1*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+1) = -1*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+(DIM+2)) = -1*L;
            }
            else if (i==DIM+1){
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 3*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-1) = -1*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+1) = -1*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-(DIM+2)) = -1*L;
            }
            else if (j==DIM+1){
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 3*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-1) = -1*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+(DIM+2)) = -1*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-(DIM+2)) = -1*L;
            }
            else {
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos) = 4*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-1) = -1*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+1) = -1*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos-(DIM+2)) = -1*L;
                conj_grad_A_proj.coeffRef(diag_pos, diag_pos+DIM+2) = -1*L;
            }

            diag_pos++;
        }
    }

    if(DIM < 5) cout << conj_grad_A_proj;
}

void FluidSolver::arrToVec(VectorXf* vec, float* arr, int size){
    for(int i=0; i<size; i++) (*vec)(i) = arr[i];
}

void FluidSolver::vecToArr(float* arr, VectorXf* vec, int size){
    for(int i=0; i<size; i++)  arr[i] = (*vec)(i);
}

void FluidSolver::add_source(int n, float *x, float *s, float dt){

    int i, size = (n+2)*(n+2);

    for(i=0; i<size; i++) {
        x[i] += dt*s[i];
        if (x[i] > 1) x[i] = 1;
    }

}

void FluidSolver::diffuse(int n, int b, float *x, float *x0, float diff, float dt){

    int i, j, k;
    float a=dt*diff*n*n;

    for(k=0; k<20; k++){
        for(i=1; i<n+1; i++) {
            for(j=1; j<n+1; j++) {
                x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/(1+4*a);
            }
        }
    }
    set_bnd(n,b,x);
}

void FluidSolver::diffuse_cg(int n, int b, float *x, float *x0, SparseMatrix<float> * conj_grad_A){

    for(int i=0; i<DIM; i++) {
        for(int j=0; j<DIM; j++){
            conj_grad_b(i*DIM+j) = x0[(i+1)*(DIM+2)+j+1];
        }

    }

    ConjugateGradient<SparseMatrix<float> > cg;
    cg.compute(*conj_grad_A);

    //cg.setTolerance(cg_tolerance);

//    cg.setMaxIterations(1);

//    int i = 0;

//    do {
//        conj_grad_x = cg.solveWithGuess(conj_grad_b, conj_grad_x);
//        i++;
//    } while(cg.info() != Success && i<CG_ITERATIONS);


    conj_grad_x = cg.solve(conj_grad_b);

    for(int i=0; i<DIM; i++) {
        for(int j=0; j<DIM; j++){
            x[(i+1)*(DIM+2)+j+1] = conj_grad_x(i*DIM+j);
        }

    }
   set_bnd(n,b,x);
}

void FluidSolver::diffuse_bad(int n, int b, float * x, float * x0, float diff, float dt){

    int i, j;
    float a=dt*diff*n*n;

    for ( i=1 ; i<=n ; i++ ) {
        for ( j=1 ; j<=n ; j++ ) {
            x[IX(i,j)] = x0[IX(i,j)] + a*(x0[IX(i-1,j)]+x0[IX(i+1,j)]+x0[IX(i,j-1)]+x0[IX(i,j+1)]-4*x0[IX(i,j)]);
        }
    }
    set_bnd (n,b,x);
}

void FluidSolver::advect ( int N, int b, float * d, float * d0, float * u, float * v, float dt )
{
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;
    dt0 = dt*N;

    for ( i=1 ; i<=N ; i++ ) {
        for ( j=1 ; j<=N ; j++ ) {
            x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];
            if (x<0.5) x=0.5; if (x>N+0.5) x=N+ 0.5; i0=(int)x; i1=i0+1;
            if (y<0.5) y=0.5; if (y>N+0.5) y=N+ 0.5; j0=(int)y; j1=j0+1;
            s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
            d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
        }
    }
    set_bnd ( N, b, d );
}

void FluidSolver::advect_back ( int N, int b, float * d, float * d0, float * u, float * v, float dt )
{
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;
    dt0 = dt*N;

    for ( i=1 ; i<=N ; i++ ) {
        for ( j=1 ; j<=N ; j++ ) {
            x = i+dt0*u[IX(i,j)]; y = j+dt0*v[IX(i,j)];
            if (x<0.5) x=0.5; if (x>N+0.5) x=N+ 0.5; i0=(int)x; i1=i0+1;
            if (y<0.5) y=0.5; if (y>N+0.5) y=N+ 0.5; j0=(int)y; j1=j0+1;
            s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
            d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
        }
    }
    set_bnd ( N, b, d );
}

void FluidSolver::advect_bfecc(int n, int b, float *d, float *d0, float* d_bfecc_1, float* d_bfecc_2, float *u, float *v, float dt){
    advect(n, b, d_bfecc_1, d0, u, v, dt);
    advect_back(n, b, d_bfecc_2, d_bfecc_1, u, v, dt);
    for(int i=0; i<SIZE; i++) d_bfecc_1[i] = d0[i] + (d0[i]-d_bfecc_2[i])/2.0;
    advect(n, b, d, d_bfecc_1, u, v, dt);
}

void FluidSolver::project ( int N, float * u, float * v, float * p, float * div )
{
    int i, j, k;
    float h;
    h = 1.0/N;
    for ( i=1 ; i<=N ; i++ ) {
        for ( j=1 ; j<=N ; j++ ) {
            div[IX(i,j)] = -0.5*h*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)]) + div_constants[IX(i,j)];
            p[IX(i,j)] = 0;
        }
    }
    set_bnd ( N, 0, div ); set_bnd ( N, 0, p );
    for ( k=0 ; k<20 ; k++ ) {
        for ( i=1 ; i<=N ; i++ ) {
            for ( j=1 ; j<=N ; j++ ) {
                p[IX(i,j)] = (div[IX(i,j)]+p[IX(i-1,j)]+p[IX(i+1,j)]+p[IX(i,j-1)]+p[IX(i,j+1)])/4;
            }
        }
    set_bnd ( N, 0, p );
    }
    for ( i=1 ; i<=N ; i++ ) {
        for ( j=1 ; j<=N ; j++ ) {
            u[IX(i,j)] -= 0.5*(p[IX(i+1,j)]-p[IX(i-1,j)])/h;
            v[IX(i,j)] -= 0.5*(p[IX(i,j+1)]-p[IX(i,j-1)])/h;
        }
    }
    set_bnd ( N, 1, u ); set_bnd ( N, 2, v );

    for(i = 0; i<SIZE; i++) {
        div_constants[i] -= div_falloff*dt;
        if (div_constants[i] < 0) div_constants[i] = 0;
    }

    for ( i=1 ; i<=N ; i++ ) {
        for ( j=1 ; j<=N ; j++ ) {
            div[IX(i,j)] = -0.5*h*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)]) + div_constants[IX(i,j)];
            p[IX(i,j)] = 0;
        }
    }
    set_bnd ( N, 0, div ); set_bnd ( N, 0, p );
}

void FluidSolver::project_cg(int n, float *u, float *v, float *p, float *div){

//    int i, j;
//    float h;
//    h = 1.0/(n+2);

//    for ( i=1 ; i<=n ; i++ ) {
//        for ( j=1 ; j<=n ; j++ ) {
//            div[IX(i,j)] = -0.5*h*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)]) + div_constants[IX(i,j)];
//        }
//    }


//    set_bnd ( n, 0, div ); set_bnd ( n, 0, p );

//    for(int i=0; i<DIM; i++) {
//        for(int j=0; j<DIM; j++){
//            conj_grad_b(i*DIM+j) = div[(i+1)*(DIM+2)+j+1];
//        }

//    }


//    ConjugateGradient<SparseMatrix<float> > cg;
//    cg.setMaxIterations(CG_ITERATIONS);
//    //cg.setTolerance(1e-2);
//    cg.compute(conj_grad_A_proj);

//    conj_grad_x = cg.solve(conj_grad_b);

//    VectorXf res = conj_grad_A_proj*conj_grad_x - conj_grad_b;
//    printf("Norm: %f, Iterations: %d, Error: %f\n", res.norm(), cg.iterations(), cg.error()    );

//    for(int i=0; i<DIM; i++) {
//        for(int j=0; j<DIM; j++){
//            p[(i+1)*(DIM+2)+j+1] = conj_grad_x(i*DIM+j);
//        }
//    }

//    set_bnd ( n, 0, p );

//    for ( i=1 ; i<=n ; i++ ) {
//        for ( j=1 ; j<=n ; j++ ) {
//            u[IX(i,j)] -= 0.5*(p[IX(i+1,j)]-p[IX(i-1,j)])/h;
//            v[IX(i,j)] -= 0.5*(p[IX(i,j+1)]-p[IX(i,j-1)])/h;
//        }
//    }
//    set_bnd ( n, 1, u ); set_bnd ( n, 2, v );

//    for(i = 0; i<SIZE; i++) {
//        div_constants[i] -= div_falloff*dt;
//        if (div_constants[i] < 0) div_constants[i] = 0;
//    }

//    for ( i=1 ; i<=n ; i++ ) {
//        for ( j=1 ; j<=n ; j++ ) {
//            div[IX(i,j)] = -0.5*h*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)]) + div_constants[IX(i,j)];
//        }
//    }


//    set_bnd ( n, 0, div ); set_bnd ( n, 0, p );

    VectorXf conj_grad_b = VectorXf(SIZE);
    VectorXf conj_grad_x = VectorXf(SIZE);

    int i, j;
    float h;
    h = 1.0/(n+2);

    for ( i=1 ; i<=n ; i++ ) {
        for ( j=1 ; j<=n ; j++ ) {
            div[IX(i,j)] = -0.5*h*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)]) + div_constants[IX(i,j)];
        }
    }


    set_bnd ( n, 0, div );

    for(int i=0; i<DIM+2; i++) {
        for(int j=0; j<DIM+2; j++){
            conj_grad_b(IX(i,j)) = div[IX(i,j)];
        }
    }

    ConjugateGradient<SparseMatrix<float> > cg;
    cg.setMaxIterations(CG_ITERATIONS);
    cg.setTolerance(1e-3);
    cg.compute(conj_grad_A_proj);

    conj_grad_x = cg.solve(conj_grad_b);

    VectorXf res = conj_grad_A_proj*conj_grad_x - conj_grad_b;
    printf("Norm: %f, Iterations: %d, Error: %f\n", res.norm(), cg.iterations(), cg.error() );

    for(int i=0; i<DIM+2; i++) {
        for(int j=0; j<DIM+2; j++){
            p[IX(i,j)] = conj_grad_x(IX(i,j));
        }
    }

    //set_bnd ( n, 0, p );

    for ( i=1 ; i<=n ; i++ ) {
        for ( j=1 ; j<=n ; j++ ) {
            u[IX(i,j)] -= 0.5*(p[IX(i+1,j)]-p[IX(i-1,j)])/h;
            v[IX(i,j)] -= 0.5*(p[IX(i,j+1)]-p[IX(i,j-1)])/h;
        }
    }
    set_bnd ( n, 1, u ); set_bnd ( n, 2, v );

    for(i = 0; i<SIZE; i++) {
        div_constants[i] -= div_falloff*dt;
        if (div_constants[i] < 0) div_constants[i] = 0;
    }

//    for ( i=1 ; i<=n ; i++ ) {
//        for ( j=1 ; j<=n ; j++ ) {
//            div[IX(i,j)] = -0.5*h*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)]) + div_constants[IX(i,j)];
//        }
//    }


    //set_bnd ( n, 0, div ); set_bnd ( n, 0, p );
}

void FluidSolver::dens_step(int n, float * x, float * x0, float * u, float * v, float diff, float dt){

    add_source(n, x, x0, dt);
    SWAP (x, x0);
    if(GAUSS_SEIDEL==method) {
        diffuse(n, 0, x, x0, diff, dt);
        SWAP (x, x0);
        if (bfecc) advect_bfecc(n, 0, x, x0, d_bfecc_1, d_bfecc_2, u, v, dt);
        else advect(n, 0, x, x0, u, v, dt);
    }
    else {
        diffuse_cg(n, 0, x, x0, &conj_grad_A_diff);
        SWAP (x, x0);
        if(bfecc) advect_bfecc(n, 0, x, x0, d_bfecc_1, d_bfecc_2, u, v, dt);
        else advect(n, 0, x, x0, u, v, dt);
    }
    memset(x0, 0, SIZE*sizeof(float));
}

void FluidSolver::temp_step(){

    add_source(DIM, t, t_prev, dt);

    if(GAUSS_SEIDEL==method) {

        if (temp_diff_enable) {
            SWAP (t, t_prev);
            diffuse(DIM, 0, t, t_prev, heat_diff, dt);
        }

        if (temp_adv_enable){
            SWAP (t, t_prev);
            if (bfecc) advect_bfecc(DIM, 0, t, t_prev, t_bfecc_1, t_bfecc_2, u, v, dt);
            else advect(DIM, 0, t, t_prev, u, v, dt);
        }
    }

    else {

        if (temp_diff_enable) {
            SWAP (t, t_prev);
            diffuse_cg(DIM, 0, t, t_prev, &conj_grad_A_temp);
        }

        if (temp_adv_enable) {
            SWAP (t, t_prev);
            if(bfecc) advect_bfecc(DIM, 0, t, t_prev, t_bfecc_1, t_bfecc_2, u, v, dt);
            else advect(DIM, 0, t, t_prev, u, v, dt);
        }

    }

    if (temp_env_enable) {
        SWAP(t, t_prev);
        environment_loss(DIM, t, t_prev, t_amb, t_max, cr, dt);
    }
    /* Clamp temps */
    for(int i=0; i<SIZE; i++) {
        if (t[i] < 0) t[i] = 0;
        else if (t[i] > 1) t[i] = 1;
    }
    memset(t_prev, 0, SIZE*sizeof(float));
}

void FluidSolver::environment_loss(int n, float *x, float *x0, float x_amb, float x_max, float cr, float dt){

    float loss;
    for(int i=0; i<SIZE; i++) {
        loss = (x0[i] - x_amb)/(x_max - x_amb);
        x[i] = x0[i] - cr*dt*loss*loss*loss*loss;
    }
}

void FluidSolver::vel_step(int n, float *u, float *v, float *u0, float *v0, float visc, float dt){

    //method = JACOBI;

    switch(method){

        case GAUSS_SEIDEL:
            SWAP ( u0, u ); diffuse(n, 1, u, u0, visc, dt);
            SWAP ( v0, v ); diffuse(n, 2, v, v0, visc, dt);
            project ( n, u, v, u0, divergence );
            SWAP ( u0, u ); SWAP ( v0, v );
            if(bfecc) {
                advect_bfecc ( n, 1, u, u0, u_bfecc_1, u_bfecc_2, u0, v0, dt );
                advect_bfecc ( n, 2, v, v0, v_bfecc_1, v_bfecc_2, u0, v0, dt );
            }
            else {
                advect ( n, 1, u, u0, u0, v0, dt ); advect ( n, 2, v, v0, u0, v0, dt );
            }
            project ( n, u, v, u0, divergence );
            break;

        case CONJ_GRAD:
            SWAP ( u0, u ); diffuse_cg ( n, 1, u, u0, &conj_grad_A_visc);
            SWAP ( v0, v ); diffuse_cg ( n, 2, v, v0, &conj_grad_A_visc);
            project_cg ( n, u, v, u0, divergence );
            SWAP ( u0, u ); SWAP ( v0, v );
            if(bfecc) {
                advect_bfecc ( n, 1, u, u0, u_bfecc_1, u_bfecc_2, u0, v0, dt );
                advect_bfecc ( n, 2, v, v0, v_bfecc_1, v_bfecc_2, u0, v0, dt );
            }
            else {
                advect ( n, 1, u, u0, u0, v0, dt ); advect ( n, 2, v, v0, u0, v0, dt );
            }
            project_cg ( n, u, v, u0, divergence );
            break;
    }

    if(sum_divergence && iterations < DIV_ITERATIONS) {

        float total_div = 0;
        for(int i=0; i<SIZE; i++){
            total_div += divergence[i]*divergence[i];
        }
        energy[iterations] = total_div;
        iterations++;
    }
    else if(sum_divergence){
        printf("OUTPUTING DIVERGENCE\n");
        ofstream energy_file;
        char val[50];
        char curr_file[512];
        sprintf(curr_file, DIVERGENCE_FILE, file_num);
        energy_file.open(curr_file);
        //energy_file << visc << "\n";
        for(int i=0; i<1000; i++){
            sprintf(val, "%f,", energy[i]);
            energy_file << val << "\n";
        }
        energy_file.close();
        file_num++;
        sum_divergence = false;
    }

//    if(sum_velocities && iterations < 1000) {
//        if(iterations == 2) {
//            iterations = 0;
//            simEnabled = false;
//            sum_velocities = false;
//        }
//        float total_energy = 0;
//        for(int i=0; i<SIZE; i++){
//            total_energy += u[i]*u[i];
//            total_energy += v[i]*v[i];
//        }
//        energy[iterations] = total_energy;
//        iterations++;
//    }
//    else if(sum_velocities){
//        printf("OUTPUTING ENERGIES");
//        ofstream energy_file;
//        char val[50];
//        energy_file.open(ENERGY_FILE);
//        energy_file << visc << "\n";
//        for(int i=0; i<1000; i++){
//            sprintf(val, "%f,", energy[i]);
//            energy_file << val << "\n";
//        }
//        energy_file.close();
//        sum_velocities = false; }
}

//void FluidSolver::set_bnd(int n, int b, float * x){
//    int i;
//    for ( i=1 ; i<=n ; i++ ) {
//        if(b==1) x[IX(0 ,i)] = -x[IX(1,i)]; else x[IX(0 ,i)] = x[IX(1,i)];
//        if(b==1) x[IX(n+1,i)] = -x[IX(n,i)]; else x[IX(n+1,i)] = x[IX(n,i)];
//        if(b==2) x[IX(i,0)] = -x[IX(i,1)]; else x[IX(i,0)] = x[IX(i,1)];
//        if(b==2) x[IX(i,n+1)] = -x[IX(i,n+1)]; else x[IX(i,n+1)] = x[IX(i,n)];
//    }
//    x[IX(0 ,0 )] = 0.5*(x[IX(1,0 )]+x[IX(0 ,1)]);
//    x[IX(0 ,n+1)] = 0.5*(x[IX(1,n+1)]+x[IX(0 ,n )]);
//    x[IX(n+1,0 )] = 0.5*(x[IX(n,0 )]+x[IX(n+1,1)]);
//    x[IX(n+1,n+1)] = 0.5*(x[IX(n,n+1)]+x[IX(n+1,n)]);
//}

void FluidSolver::set_bnd ( int N, int b, float * x )
{
int i;

for ( i=1 ; i<=N ; i++ ) {
    if(b==1){
        x[IX(0,i)] = -x[IX(1,i)];
        x[IX(N+1,i)] = -x[IX(N,i)];
    } else {
        x[IX(0 ,i)] = x[IX(1,i)];
        x[IX(N+1,i)] = x[IX(N,i)];

    }
    if(b==2){
        x[IX(i,0 )] =-x[IX(i,1)];
        x[IX(i,N+1)] =-x[IX(i,N)];
    }
    else{
        x[IX(i,0 )] =x[IX(i,1)];
        x[IX(i,N+1)] = x[IX(i,N)];
    }

}
x[IX(0 ,0 )] = 0.5*(x[IX(1,0 )]+x[IX(0 ,1)]);
x[IX(0 ,N+1)] = 0.5*(x[IX(1,N+1)]+x[IX(0 ,N )]);
x[IX(N+1,0 )] = 0.5*(x[IX(N,0 )]+x[IX(N+1,1)]);
x[IX(N+1,N+1)] = 0.5*(x[IX(N,N+1)]+x[IX(N+1,N )]);
}


void FluidSolver::step(){

    isDone = false;
    if(simEnabled || oneStep){
        external_forces();
        vel_step ( DIM, u, v, u_prev, v_prev, visc, dt );
        update_particles();
        cfl_check();
        dens_step(DIM, d, d_prev, u, v, diff, dt);
        temp_step();
    }
    oneStep = false;
    isDone = true;
}

void FluidSolver::diffusion_test() {

    reset();

    int i,j;

    // A square that is one third of the window width/height
    for(i=DIM/3; i<=2*DIM/3; i++){
        for(j=DIM/3; j<=2*DIM/3.0; j++){
            d[IX(i,j)] = 1;
        }
    }

    iterations = 0;
    //simEnabled =true;
}

void FluidSolver::advection_test() {

    reset();

    u[IX(DIM/2, DIM/2)]  = 100;
    v[IX(DIM/2, DIM/2)]  = 0;

    int i,j;

    // A square that is one third of the window width/height
    for(i=DIM/3; i<=2*DIM/3; i++){
        for(j=DIM/3; j<=2*DIM/3.0; j++){
            //d_prev[IX(i,j)] = 1;
            d[IX(i,j)] = 1;
        }
    }

    sum_velocities = true;
    iterations = 0;
}

void FluidSolver::cig_test() {
    reset();
    u[IX(1,DIM/2)]  = 100;
}

void FluidSolver::reset(){
    memset(d_prev, 0, SIZE*sizeof(float));
    memset(u_prev, 0, SIZE*sizeof(float));
    memset(v_prev, 0, SIZE*sizeof(float));
    memset(t_prev, 0, 2*SIZE*sizeof(float));
    memset(d, 0, SIZE*sizeof(float));
    memset(u, 0, SIZE*sizeof(float));
    memset(v, 0, SIZE*sizeof(float));
    memset(t, 0, 2*SIZE*sizeof(float));
    memset(div_constants, 0, SIZE*sizeof(float));
    memset(u_acc, 0, SIZE*sizeof(float));
    memset(v_acc, 0, SIZE*sizeof(float));
    memset(particle_temp, 0, num_particles*sizeof(float));
    num_particles = 0;
    //d0 = 0;
}

float FluidSolver::myAbs(float x) {
    if (x < 0) return -x;
    return x;
}

float FluidSolver::max_vel(){
    int i;
    float vert, horz, max = 0;
    for(i=0; i<SIZE; i++){
        vert = myAbs(u[i]);
        horz = myAbs(v[i]);
        if (vert > max) max = vert;
        if (horz > max) max = horz;
    }
    return max;
}

void FluidSolver::set_dt(float new_dt) {
    dt = new_dt;
    initializeConjGradMat(&conj_grad_A_diff, diff, dt, h_space);
    initializeConjGradMat(&conj_grad_A_visc, visc, dt, h_space);
    initializeConjGradMat(&conj_grad_A_temp, diff, dt, h_space);
    initializeConjGradProj();
}

void FluidSolver::set_diff(float new_diff) {
    diff = new_diff;
    initializeConjGradMat(&conj_grad_A_diff, diff, dt, h_space);
}

void FluidSolver::set_heat_diff(float new_diff){
    heat_diff = new_diff;
    initializeConjGradMat(&conj_grad_A_temp, heat_diff, dt, h_space);
}

void FluidSolver::set_visc(float new_visc) {
    visc = new_visc;
    initializeConjGradMat(&conj_grad_A_visc, visc, dt, h_space);
}

void FluidSolver::scalar_field(int field){
    float maxVel = max_vel();

    if(field == 1){
        memcpy(scalars, v, SIZE*sizeof(float));
        for(int i=0; i<SIZE; i++) scalars[i] /= maxVel;
    }
    // Velocities (Vert)
    else if (field ==2){
        memcpy(scalars, u, SIZE*sizeof(float));
        for(int i=0; i<SIZE; i++) scalars[i] /= maxVel;
    }
    // Divergence
    else if (field==3){
        float max = 0;
        for(int i=0; i<SIZE; i++) if (myAbs(divergence[i]) > max) max = myAbs(divergence[i]);
        memcpy(scalars,divergence, SIZE*sizeof(float));
        for(int i=0; i<SIZE; i++) scalars[i] /= max;

        printf("%f\n", max);
    }

}

void FluidSolver::vel_field(int field){
    memset(u, 0, sizeof(float)*SIZE);
    memset(v, 0, sizeof(float)*SIZE);
    switch (field){
    case 1:
        for(int i=0; i<SIZE; i++) u[i] = .5;
        break;
    case 2:
        for(int i=0; i<SIZE; i++) v[i] = .005;
        break;
    case 3:
        for(int i=0; i<SIZE; i++) u[i] = -.5;
        break;
    case 4:
        for(int i=0; i<SIZE; i++) v[i] = -.5;
        break;
    case 5:
        for(int i=0; i<SIZE; i++) {
            v[i] = ((float)rand())/INT_MAX-.5;
            u[i] = ((float)rand())/INT_MAX-.5;
        }
        break;
    }

}

void FluidSolver::update_particles(){

    float pos_horz, pos_vert, weight_horz, weight_vert;
    float u_vel = 0, v_vel = 0, temp = 0;
    int grid_i = 0, grid_j = 0;

    int i0, j0;

    switch(vel_method) {

    case FORWARD_EULER:

        for(int i=0; i<num_particles; i++){

            grid_j = (int)(particle_pos_X[i]/h_space);
            grid_i = (int)(particle_pos_Y[i]/h_space);

            pos_horz = particle_pos_X[i]/h_space - grid_j;
            pos_vert = particle_pos_Y[i]/h_space - grid_i;

            if(num_particles <10){
            printf("horz: %f\n", pos_horz);
            printf("vert: %f\n", pos_vert);
            }

            if(pos_horz > 0.5) {
                j0 = grid_j+1;
                weight_horz = 1 - (pos_horz-0.5);
            }
            else {
                j0 = grid_j-1;
                weight_horz = 1 - (0.5-pos_horz);
            }
            if(pos_vert > 0.5) {
                i0 = grid_i+1;
                weight_vert = 1 - (pos_vert-0.5);
            }
            else {
                i0 = grid_i - 1;
                weight_vert = 1 - (0.5-pos_vert);
            }
            if (i0 < 0 || i0 > DIM+1) i0 = grid_i;
            if (j0 < 0 || j0 > DIM+1) j0 = grid_j;

            if(num_particles <10){
            printf("weight horz: %f\n", weight_horz);
            printf("weight vert: %f\n", weight_vert);
            }

            u_vel =    weight_vert*((weight_horz)*u[IX(grid_j, grid_i)]+(1-weight_horz)*u[IX(j0, grid_i)]) +
                   (1-weight_vert)*((weight_horz)*u[IX(grid_j, i0)]+(1-weight_horz)*u[IX(j0, i0)]);

            v_vel =  weight_vert*((weight_horz)*v[IX(grid_j, grid_i)]+(1-weight_horz)*v[IX(j0, grid_i)]) +
                    (1-weight_vert)*((weight_horz)*v[IX(grid_j, i0)]+(1-weight_horz)*v[IX(j0, i0)]);

            particle_pos_X[i] += dt*u_vel;
            particle_pos_Y[i] += dt*v_vel;


            if(particle_pos_X[i]>=1) particle_pos_X[i] = .99;
            else if (particle_pos_X[i]<=0) particle_pos_X[i] = .01;
            if(particle_pos_Y[i]>=1) particle_pos_Y[i] = .99;
            else if (particle_pos_Y[i]<=0) particle_pos_Y[i] = .01;

            if(!exploding_particles || (exploding_particles && particle_temp[i] < temp_threshold)){

                temp =  weight_vert*((weight_horz)*t[IX(grid_j, grid_i)]+(1-weight_horz)*t[IX(j0, grid_i)]) +
                        (1-weight_vert)*((weight_horz)*t[IX(grid_j, i0)]+(1-weight_horz)*t[IX(j0, i0)]);

                particle_temp[i] += dt*ah*r*r*(temp-particle_temp[i])/cm;
                if(particle_temp[i]>1) particle_temp[i] = 1;
                else if (particle_temp[i]<0) particle_temp[i] = 0;
            }
            else if (exploding_particles && particle_life[i] > 0){
                float particle_heat = Bh*z*dt;
                float particle_div = .05*dt;
                float h0, h1, v0, v1;
                int i0=grid_i, i1, j0=grid_j, j1;

                if(pos_horz > 0.5){
                    j1 = j0 + 1;
                    h0 = 1 - (pos_horz-.5);
                } else {
                    j1 = j0 - 1;
                    h0 = 1 - (.5-pos_horz);
                }
                h1 = 1 - h0;

                if(pos_vert > 0.5){
                    i1 = i0 + 1;
                    v0 = 1 - (pos_vert-.5);
                } else {
                    i1 = i0 - 1;
                    v0 = 1 - (.5-pos_vert);
                }
                v1 = 1 - v0;

                if (i1 < 0 || i1 > DIM+1) i1 = i0;
                if (j1 < 0 || j1 > DIM+1) j1 = j0;

                t[IX(j0,i0)] = particle_heat*h0*v0;
                t[IX(j0,i1)] = particle_heat*h0*v1;
                t[IX(j1,i1)] = particle_heat*h1*v1;
                t[IX(j1,i0)] = particle_heat*h1*v0;
//                t[IX(j0,i0)] = particle_heat/4;
//                t[IX(j0,i1)] = particle_heat/4;
//                t[IX(j1,i1)] = particle_heat/4;
//                t[IX(j1,i0)] = particle_heat/4;

                if(particle_life[i] == 1 && explosiveParticles) {
                    div_constants[IX(j0,i0)] = particle_div*h0*v0;
                    div_constants[IX(j0,i1)] = particle_div*h0*v1;
                    div_constants[IX(j1,i1)] = particle_div*h1*v1;
                    div_constants[IX(j1,i0)] = particle_div*h1*v0;
//                    div_constants[IX(j0,i0)] = particle_div/4;
//                    div_constants[IX(j0,i1)] = particle_div/4;
//                    div_constants[IX(j1,i1)] = particle_div/4;
//                    div_constants[IX(j1,i0)] = particle_div/4;
                }

                particle_life[i] -= z*dt;
            }

        }
        break;

    case MIDPOINT:

        float x_pos, y_pos;

        for(int i=0; i<num_particles; i++){

            grid_i = (int)(particle_pos_X[i]/h_space);
            grid_j = (int)(particle_pos_Y[i]/h_space);

            weight_vert = particle_pos_X[i]/h_space - grid_i;
            weight_horz = particle_pos_Y[i]/h_space - grid_j;

            u_vel =    weight_vert*((weight_horz)*u[IX(grid_i, grid_j)]+(1-weight_horz)*u[IX(grid_i, grid_j+1)]) +
                    (1-weight_vert)*((weight_horz)*u[IX(grid_i+1, grid_j)]+(1-weight_horz)*u[IX(grid_i+1, grid_j+1)]);

            v_vel = weight_vert*((weight_horz)*v[IX(grid_i, grid_j)]+(1-weight_horz)*v[IX(grid_i, grid_j+1)]) +
                 (1-weight_vert)*((weight_horz)*v[IX(grid_i+1, grid_j)]+(1-weight_horz)*v[IX(grid_i+1, grid_j+1)]);

            if(!exploding_particles || (exploding_particles && particle_temp[i] < temp_threshold)){
                temp = weight_vert*((weight_horz)*t[IX(grid_i, grid_j)]+(1-weight_horz)*t[IX(grid_i, grid_j+1)]) +
                        (1-weight_vert)*((weight_horz)*t[IX(grid_i+1, grid_j)]+(1-weight_horz)*t[IX(grid_i+1, grid_j+1)]);
                particle_temp[i] += dt*ah*r*r*(temp-particle_temp[i])/cm;
                if(particle_temp[i]>1) particle_temp[i] = 1;
                else if (particle_temp[i]<0) particle_temp[i] = 0;
            }
            else if (exploding_particles && particle_life[i] > 0){
                particle_life[i] -= z*dt;
                float particle_heat = Bh*z*dt;
                int i0, j0;
                if(weight_horz >= h_space/2) i0 = grid_i+1; else i0 = grid_i-1;
                if(weight_vert >= h_space/2) j0 = grid_j+1; else j0 = grid_j-1;
                if(i0 < 0) i0 = 0; if(i0>DIM+2) i0=DIM+2;
                if(j0 < 0) j0 = 0; if(j0>DIM+2) j0=DIM+2;
                t[IX(grid_i,grid_j)] += particle_heat/4.0;
                t[IX(grid_i,j0)] += particle_heat/4.0;
                t[IX(i0,j0)] += particle_heat/4.0;
                t[IX(i0,grid_j)] += particle_heat/4.0;
            }

            x_pos = particle_pos_X[i]+dt*u_vel*0.5;
            y_pos = particle_pos_Y[i]+dt*v_vel*0.5;


            grid_i = (int)(x_pos/h_space);
            grid_j = (int)(y_pos/h_space);

            weight_horz = x_pos/h_space - grid_i;
            weight_vert = y_pos/h_space - grid_j;

            u_vel =    weight_vert*((weight_horz)*u[IX(grid_i, grid_j)]+(1-weight_horz)*u[IX(grid_i, grid_j+1)]) +
                    (1-weight_vert)*((weight_horz)*u[IX(grid_i+1, grid_j)]+(1-weight_horz)*u[IX(grid_i+1, grid_j+1)]);

            v_vel = weight_vert*((weight_horz)*v[IX(grid_i, grid_j)]+(1-weight_horz)*v[IX(grid_i, grid_j+1)]) +
                 (1-weight_vert)*((weight_horz)*v[IX(grid_i+1, grid_j)]+(1-weight_horz)*v[IX(grid_i+1, grid_j+1)]);

            particle_pos_X[i] += dt*u_vel;
            particle_pos_Y[i] += dt*v_vel;

            if(particle_pos_X[i]>1) particle_pos_X[i] = 1;
            else if (particle_pos_X[i]<0) particle_pos_X[i] = 0;
            if(particle_pos_Y[i]>1) particle_pos_Y[i] = 1;
            else if (particle_pos_Y[i]<0) particle_pos_Y[i] = 0;
        }
        break;
    }
}

void *FluidSolver::update_particles_ST(void * arg){

    struct thread_data * data;
    data = (struct thread_data *)arg;

    float weight_horz=0, weight_vert=0;
    float u_vel = 0, v_vel = 0;
    int grid_i = 0, grid_j = 0;

    int begin_index = data->tid*data->solver->num_particles/NUM_THREADS;
    int end_index = (data->tid+1)*data->solver->num_particles/NUM_THREADS;

    switch(data->solver->vel_method) {

    case FORWARD_EULER:

        for(int i=begin_index; i<end_index; i++){

            grid_i = (int)(data->solver->particle_pos_X[i]/data->solver->h_space);
            grid_j = (int)(data->solver->particle_pos_Y[i]/data->solver->h_space);

            weight_horz = data->solver->particle_pos_X[i]/data->solver->h_space - grid_i;
            weight_vert = data->solver->particle_pos_Y[i]/data->solver->h_space - grid_j;


            u_vel =    weight_vert*((weight_horz)*data->solver->u[IX(grid_i, grid_j)]+(1-weight_horz)*data->solver->u[IX(grid_i, grid_j+1)]) +
                    (1-weight_vert)*((weight_horz)*data->solver->u[IX(grid_i+1, grid_j)]+(1-weight_horz)*data->solver->u[IX(grid_i+1, grid_j+1)]);

            v_vel = weight_vert*((weight_horz)*data->solver->v[IX(grid_i, grid_j)]+(1-weight_horz)*data->solver->v[IX(grid_i, grid_j+1)]) +
                 (1-weight_vert)*((weight_horz)*data->solver->v[IX(grid_i+1, grid_j)]+(1-weight_horz)*data->solver->v[IX(grid_i+1, grid_j+1)]);

            data->solver->particle_pos_X[i] += data->solver->dt*u_vel;
            data->solver->particle_pos_Y[i] += data->solver->dt*v_vel;

            if (data->solver->particle_pos_X[i] > 1) data->solver->particle_pos_X[i] = 1;
            else if (data->solver->particle_pos_X[i] < 0 ) data->solver->particle_pos_X[i] = 0;
            if (data->solver->particle_pos_Y[i] > 1) data->solver->particle_pos_Y[i] = 1;
            else if (data->solver->particle_pos_Y[i] < 0 ) data->solver->particle_pos_Y[i] = 0;
        }
        break;

    case MIDPOINT:

        float x_pos, y_pos;

        for(int i=begin_index; i<end_index; i++){

            grid_i = (int)(data->solver->particle_pos_X[i]/data->solver->h_space);
            grid_j = (int)(data->solver->particle_pos_Y[i]/data->solver->h_space);

            weight_horz = data->solver->particle_pos_X[i]/data->solver->h_space - grid_i;
            weight_vert = data->solver->particle_pos_Y[i]/data->solver->h_space - grid_j;

            u_vel =    weight_vert*((1-weight_horz)*data->solver->u[IX(grid_i, grid_j)]+(weight_horz)*data->solver->u[IX(grid_i, grid_j+1)]) +
                    (1-weight_vert)*((1-weight_horz)*data->solver->u[IX(grid_i+1, grid_j)]+(weight_horz)*data->solver->u[IX(grid_i+1, grid_j+1)]);

            v_vel = weight_vert*((weight_horz)*data->solver->v[IX(grid_i, grid_j)]+(1-weight_horz)*data->solver->v[IX(grid_i, grid_j+1)]) +
                 (1-weight_vert)*((weight_horz)*data->solver->v[IX(grid_i+1, grid_j)]+(1-weight_horz)*data->solver->v[IX(grid_i+1, grid_j+1)]);

            x_pos = data->solver->particle_pos_X[i]+data->solver->dt*u_vel*0.5;
            y_pos = data->solver->particle_pos_Y[i]+data->solver->dt*v_vel*0.5;


            grid_i = (int)(x_pos/data->solver->h_space);
            grid_j = (int)(y_pos/data->solver->h_space);

            weight_horz = x_pos/data->solver->h_space - grid_i;
            weight_vert = y_pos/data->solver->h_space - grid_j;

            u_vel =    weight_vert*((1-weight_horz)*data->solver->u[IX(grid_i, grid_j)]+(weight_horz)*data->solver->u[IX(grid_i, grid_j+1)]) +
                    (1-weight_vert)*(1-(weight_horz)*data->solver->u[IX(grid_i+1, grid_j)]+(weight_horz)*data->solver->u[IX(grid_i+1, grid_j+1)]);

            v_vel = weight_vert*((1-weight_horz)*data->solver->v[IX(grid_i, grid_j)]+(weight_horz)*data->solver->v[IX(grid_i, grid_j+1)]) +
                 (1-weight_vert)*((1-weight_horz)*data->solver->v[IX(grid_i+1, grid_j)]+(weight_horz)*data->solver->v[IX(grid_i+1, grid_j+1)]);

            data->solver->particle_pos_X[i] += data->solver->dt*u_vel;
            data->solver->particle_pos_Y[i] += data->solver->dt*v_vel;

            if (data->solver->particle_pos_X[i] > 1) data->solver->particle_pos_X[i] = 1;
            else if (data->solver->particle_pos_X[i] < 0 ) data->solver->particle_pos_X[i] = 0;
            if (data->solver->particle_pos_Y[i] > 1) data->solver->particle_pos_Y[i] = 1;
            else if (data->solver->particle_pos_Y[i] < 0 ) data->solver->particle_pos_Y[i] = 0;

        }
        break;
    }
    pthread_exit(NULL);
}

void FluidSolver::update_particles_MT(){

    pthread_t threads[NUM_THREADS];
    struct thread_data data[NUM_THREADS];
    int i, rc;

    printf("HERE");

    for(i=0; i<NUM_THREADS; i++){

        data[i].tid = i;
        data[i].solver = this;

        rc = pthread_create(&threads[i], NULL, &FluidSolver::update_particles_ST, &data[i]);

        if (rc){
             cout << "Error:unable to create thread," << rc << endl;
             exit(-1);
        }
    }
    pthread_exit(NULL);
}

void FluidSolver::buoyancy_test(){
    reset();
    for(int j=0; j<DIM+2; j++) {
        for(int i=0; i<(DIM+2)/2; i++) {
            t[IX(i,j)] = 1;
            d[IX(i,j)] = .9;
        }
    }
    t_amb = 0.0;
}

void FluidSolver::external_forces(){
    if(buoyancy_enabled){
        for(int i=0; i<SIZE; i++) {
            u_acc[i] = buoy*(t[i]-t_amb);
            u[i] += u_acc[i]*dt;
            v[i] += v_acc[i]*dt;
        }
    }
}

void FluidSolver::set_blast_rad(float radius){
    blast_radius = radius;
}

void FluidSolver::set_blast_intensity(float intensity){
    blast_intensity = intensity;
}

void FluidSolver::buoyancyEnable(bool enabled){
    buoyancy_enabled = enabled;
}

void FluidSolver::cfl_check(){
    cfl_violations = 0;
    for(int i=0; i<SIZE; i++) {
        if(u[i]*dt/h_space > CFL_LIMIT) cfl_violations++;
        if(v[i]*dt/h_space > CFL_LIMIT) cfl_violations++;
    }
}

void FluidSolver::set_buoy(float new_buoy) {
    buoy = new_buoy;
}

void FluidSolver::setExplodingParticles(bool enabled) {
    explosiveParticles = enabled;
}
