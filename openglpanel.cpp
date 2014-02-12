#include "openglpanel.h"
#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <fstream>

#include <math.h>
#define PI 3.14

#include <boost/log/trivial.hpp>

#define DENSITY_FILE "/Users/jason/densities.csv"
#define VELOCITY_FILE "/Users/jason/velocities.csv"
#define FRAME_FILE "/Users/jason/Desktop/Sim_Frames/sim_%06d.png"

#define NUM_STEPS 1

float red[300000];
float green[300000];
float blue[300000];

using namespace std;

int currentScalar = 1;
int currentFrame = 0;

OpenGLPanel::OpenGLPanel(QWidget *parent) :
    QGLWidget(parent)
{
    setFormat(QGL::DoubleBuffer | QGL::DepthBuffer);

    allocateMemory();


    solver = new FluidSolver();

    float fps = 30;
    timer = new QTimer(this);
    timer->start(1/fps*1000);

    du = 0;
    dv = 0;

    simEnabled = true;

    showGrid =  true;
    showDensities = false;
    showVelocities = false;
    showParticles = true;
    showTemperatures = false;

    buoyancyColoring = false;

    insertDensities = false;
    insertVelocities = false;
    insertHeat = false;
    ignite = false;

    addDensities = false;
    subtractDensities = false;

    boom_trigger = false;
    boom = false;

    connect(timer, SIGNAL(timeout()), this, SLOT(updateGL()));

    setMouseTracking(true);

    record = false;
}

OpenGLPanel::~OpenGLPanel(){
    deallocateMemory();
}

void OpenGLPanel::allocateMemory(){
    sources = (float *)calloc(SIZE, sizeof(float));
    if(NULL == sources) {
        //handle error
        exit(-1);
    }
}

void OpenGLPanel::deallocateMemory(){
    if(NULL != sources) free(sources);
    delete solver;
}

void OpenGLPanel::initializeGL()
{
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_COLOR_MATERIAL);
    glEnable(GL_BLEND);
    glEnable(GL_POLYGON_SMOOTH);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(0, 0, 0, 0);
}


void OpenGLPanel::resizeGL(int w, int h) {

    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, w, 0, h); // set origin to bottom left corner
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void OpenGLPanel::paintGL()
{
    int i, j;
    float c = 0, horz, vert, max;
    float hSpace = width()/(float)(DIM+2), vSpace = height()/(float)(DIM+2);

    glClear(GL_COLOR_BUFFER_BIT);


    for(i=0; i<NUM_STEPS; i++){
        if(mouseDown){
            int x = (int)(mouse_x/(hSpace));
            int y = (int)(DIM+2-mouse_y/(vSpace));
            if (insertDensities) {
                solver->d_prev[IX(y,x)] = 1000;
            }
            if (insertVelocities) {
                float u = du*solver->dt;
                float v = dv*solver->dt;
//                if( u > 1) u = 1; else if (u < -1) u = -1;
//                if( v > 1) v = 1; else if (v < -1) v = -1;
                solver->u[IX(y,x)] = u;
                solver->v[IX(y,x)] = v;
            }
            if (insertHeat){
                solver->t[IX(y,x)] += 1000*solver->dt;
                if(solver->t[IX(y,x)] > 1) solver->t[IX(y,x)] = 1;
            }
            if(ignite){
                solver->t[IX(y,x)] = 100;
            }
            if(addDensities) {
                solver->d[IX(y,x)] += 1*solver->dt;
                if(solver->d[IX(y,x)] > 1) solver->d[IX(y,x)] = 1;
            }
            if(subtractDensities){
                solver->d[IX(y,x)] -= 1*solver->dt;
                if (solver->d[IX(y,x)] < 0) solver->d[IX(y,x)] = 0;
            }

        }

        if(solver->simEnabled || solver->oneStep) explosions.clear();

        solver->step();

        char violations[20];
        sprintf(violations, "%d", solver->cfl_violations);
        emit(cflViolations(violations));

    }

    glBegin(GL_QUADS);

    glColor3f(c,c,c);
    glVertex2d(0,0);
    glVertex2d(0,height());
    glVertex2d(width(),height());
    glVertex2d(width(),0);

    if(showDensities && !showScalars){
        for(i=0; i<DIM+1; i++) {
            for(j=0; j<DIM+1; j++){
                glColor4f(1,1,1,solver->d[IX(i,j)]); glVertex2d(j*hSpace, i*vSpace);
                glColor4f(1,1,1,solver->d[IX(i+1,j)]); glVertex2d(j*hSpace, (i+1)*vSpace);
                glColor4f(1,1,1,solver->d[IX(i+1,j+1)]); glVertex2d((j+1)*hSpace, (i+1)*vSpace);
                glColor4f(1,1,1,solver->d[IX(i,j+1)]); glVertex2d((j+1)*hSpace, i*vSpace);
            }
        }
    }

    if(showTemperatures && !showScalars){

        float d, d0, r, g=0, b=0;
        float alpha = 1;
        d0 = solver->t_amb;
        for(i=0; i<DIM+2; i++) {
            for(j=0; j<DIM+2; j++){

                d = solver->t[IX(i,j)];
                r = (d-d0)/(1-d0);
                glColor4f(r,g,b,alpha);
                glVertex2d(j*hSpace, i*vSpace);

                d = solver->t[IX(i+1,j)];
                r = (d-d0)/(1-d0);
                glColor4f(r,g,b,alpha);
                glVertex2d(j*hSpace, (i+1)*vSpace);

                d = solver->t[IX(i+1,j+1)];
                r = (d-d0)/(1-d0);
                glColor4f(r,g,b,alpha);
                glVertex2d((j+1)*hSpace, (i+1)*vSpace);

                d = solver->t[IX(i,j+1)];
                r = (d-d0)/(1-d0);
                glColor4f(r,g,b,alpha);
                glVertex2d((j+1)*hSpace, i*vSpace);
            }
        }
    }


    glEnd();

    float value;
    if(showScalars){
        glBegin(GL_QUADS);
        solver->scalar_field(currentScalar);
        for(i=0; i<DIM+1; i++) {
            for(j=0; j<DIM+1; j++){
                value = solver->scalars[IX(i,j)];
                if(value < 0){
                    glColor4f(-value,0,0,1);
                    glVertex2d(j*hSpace, i*vSpace);
                    glVertex2d(j*hSpace, (i+1)*vSpace);
                    glVertex2d((j+1)*hSpace, (i+1)*vSpace);
                    glVertex2d((j+1)*hSpace, i*vSpace);
                }
                else {
                    glColor4f(0,0,value,1);
                    glVertex2d(j*hSpace, i*vSpace);
                    glVertex2d(j*hSpace, (i+1)*vSpace);
                    glVertex2d((j+1)*hSpace, (i+1)*vSpace);
                    glVertex2d((j+1)*hSpace, i*vSpace);
                }
            }
        }
        glEnd();
    }

    glBegin(GL_LINES);
    glColor3f(1,0,0);
    if(showVelocities && !showScalars){
        max = solver->max_vel();
        for(i=0; i<DIM+2; i++) {
            for(j=0; j<DIM+2; j++){
                horz = solver->v[IX(i,j)]/max;
                vert = solver->u[IX(i,j)]/max;
                glVertex2d(j*hSpace+hSpace/2, i*vSpace+vSpace/2);
                glVertex2d(j*hSpace+hSpace/2+horz*hSpace/2, i*vSpace+vSpace/2+vert*vSpace/2);
            }
        }
    }
    glEnd();
    //Reset sources
    memset(sources,0,SIZE*sizeof(float));

    //Particles
    if(showParticles && !showScalars){

//        glColor3f(1.0f,0.8f,0.0f);
            glBegin(GL_POINTS);
            float x, y;
            for(int i=0; i<solver->num_particles; i++){
                float r,g,b,life = solver->particle_life[i];
                if (life > 0){
                    if(life == 1) {r = 1; g = solver->particle_temp[i]; b=0;}
                    else{r=life, g=life, b=life;}
                    glColor3f(r, g, b);
                    //glColor3f(1.0f,rand()/(1.0*INT_MAX), 0);
                    x = this->width()*solver->particle_pos_Y[i];
                    y = this->height()*solver->particle_pos_X[i];
                    glVertex2d(x, y);
                }
            }
            glEnd();
//        } else {
//            GLUquadric *quadric = gluNewQuadric();
//            gluQuadricDrawStyle(quadric, GLU_FILL);

//            for(int i=0; i<solver->num_particles; i++){
//                glColor3f(1, green[i], 0);
//                //glColor3f(rand()/(1.0f*INT_MAX),rand()/(1.0f*INT_MAX),rand()/(1.0f*INT_MAX));
//                glPushMatrix();
//                glTranslatef(this->width()*solver->particle_pos_Y[i], this->height()*solver->particle_pos_X[i], 0);
//                gluDisk(quadric, 0, 3, 4 , 4);
//                glPopMatrix();
//            }
//        }
    }

    if(boom){
        int x, center_x = (int)(mouse_x/(width()/DIM)+0.5);
        int y, center_y = (int)(DIM-mouse_y/(height()/DIM)+0.5);
        int radius = (int)solver->blast_radius;

        if(solver->blast_shape == CIRCLE){
            for(x = center_x-radius; x<=center_x+radius; x++){
                for(y = center_y-radius; y<=center_y+radius; y++){
                    if (x >=0 && x <= DIM+1 && y >= 0 && y <= DIM+1 &&
                            ((x-center_x)*(x-center_x) + (y-center_y)*(y-center_y))<radius*radius) {
//                        solver->div_constants[IX(y,x)] = solver->blast_intensity;
                        solver->t[IX(y,x)] = 10000*solver->blast_intensity;
                        //solver->d[IX(y,x)] = 1;
                    }
                }
            }
        }
        else if(solver->blast_shape == DIAMOND){
            int y_radius = 0;
            for(x = center_x-radius; x<=center_x+radius; x++){

                for(y = center_y-y_radius; y<=center_y+y_radius; y++){
                    if (x >=0 && x <= DIM+1 && y >= 0 && y <= DIM+1) {
                        //solver->div_constants[IX(y,x)] = solver->blast_intensity;
                        solver->t[IX(y,x)] = 10000*solver->blast_intensity;
                    }
                }
                if(x < center_x) y_radius++;
                else y_radius--;
            }
        }
        else if(solver->blast_shape == SQUARE){
            for(x = center_x-radius; x<=center_x+radius; x++){
                for(y = center_y-radius; y<=center_y+radius; y++){
                    if (x >=0 && x <= DIM+1 && y >= 0 && y <= DIM+1) {
                        //solver->div_constants[IX(y,x)] = solver->blast_intensity;
                        solver->t[IX(y,x)] = 10000*solver->blast_intensity;
                    }
                }
            }
        }
        boom = false;
    }


    // Boom positions
    GLUquadric *quadric = gluNewQuadric();
    gluQuadricDrawStyle(quadric, GLU_FILL);
    int rad;
    glColor3f(1, 1, 0);


    list<BoomPosition>::iterator it;
    for(it = explosions.begin(); it != explosions.end(); it++){
        if(it->shape == CIRCLE){
            glPushMatrix();
            glTranslatef(it->x, this->height()-it->y, 0);
            gluDisk(quadric, 0, it->radius*2, 30, 30);
            glPopMatrix();
        }
        else if(it->shape == DIAMOND) {
            glBegin(GL_QUADS);
            rad = 2*it->radius;
            glVertex2f(it->x-rad, this->height()-it->y);
            glVertex2f(it->x, this->height()-it->y-rad);
            glVertex2f(it->x+rad, this->height()-it->y);
            glVertex2f(it->x, this->height()-it->y+rad);
            glEnd();

        }
        else if(it->shape == SQUARE) {
            glBegin(GL_QUADS);
            rad = 2*it->radius;
            glVertex2f(it->x-rad, this->height()-it->y-rad);
            glVertex2f(it->x-rad, this->height()-it->y+rad);
            glVertex2f(it->x+rad, this->height()-it->y+rad);
            glVertex2f(it->x+rad, this->height()-it->y-rad);
            glEnd();
        }
    }

    if(showGrid) drawGrid(DIM+2,DIM+2);

    emit mouseX(QString::number((int)mouse_x));
    emit mouseY(QString::number((int)mouse_y));

    if(solver->sum_divergence){
        emit recordingDivergence("Recording Divergence...");
        emit recordingIterations(QString("Iterations: ")+QString::number(solver->iterations));
    } else if (!solver->sum_divergence && !key_pressed) {
        emit recordingDivergence("");
        emit recordingIterations("");
    }

    if(record){
        char curr_frame[512];
        sprintf(curr_frame, FRAME_FILE, currentFrame);
        frame = grabFrameBuffer();
        frame.save(curr_frame);
        currentFrame++;
    }
}

void OpenGLPanel::drawGrid(int rows, int cols) {

    float width = this->width();
    float height = this->height();
    float vSpace = height/rows;
    float hSpace = width/cols;

    glBegin(GL_LINES);
    glLineWidth(0.1);
    glColor3f(0.3f,0.3f,0.3f);

    for(float i=vSpace; i<height; i += vSpace){
        glVertex2f(0, i);
        glVertex2f(width, i);
    }
    for(float j=hSpace; j<width; j += hSpace){
        glVertex2f(j, 0);
        glVertex2f(j, height);
    }

    glEnd();
}

void OpenGLPanel::setGrid(bool enabled) {
    showGrid = enabled;
}

void OpenGLPanel::mouseMoveEvent(QMouseEvent * e){
    float oldX = mouse_x;
    float oldY = mouse_y;
    mouse_x = e->x();
    mouse_y = e->y();
    if(mouse_x < 0) {mouse_x = 0; mouse_y=oldY;}
    if(mouse_y < 0) {mouse_y = 0; mouse_x=oldX;}
    if(mouse_x > width()) {mouse_x = width(); mouse_y = oldY;}
    if(mouse_y > height()) {mouse_y = height(); mouse_x = oldX;}
    du = -(mouse_y - oldY);
    dv = mouse_x - oldX;
}

void OpenGLPanel::mousePressEvent(QMouseEvent * e){
    float oldX = mouse_x;
    float oldY = mouse_y;
    mouse_x = e->x();
    mouse_y = e->y ();
    if(mouse_x < 0) {mouse_x = 0; mouse_y=oldY;}
    if(mouse_y < 0) {mouse_y = 0; mouse_x=oldX;}
    if(mouse_x > width()) {mouse_x = width(); mouse_y = oldY;}
    if(mouse_y > height()) {mouse_y = height(); mouse_x = oldX;}
    mouseDown = true;
}

void OpenGLPanel::mouseReleaseEvent(QMouseEvent *e){
    mouseDown = false;
    insertDensities = false;
    insertVelocities = false;
    if(boom_trigger) {
        boom = true;
        BoomPosition p = {mouse_x, mouse_y, solver->blast_radius, solver->blast_shape};
        explosions.push_front(p);
    }
}

void OpenGLPanel::keyPressEvent(QKeyEvent * e){
    key_pressed = true;
    if(e->key() == (int)'d' || e->key() == (int)'D') {
        insertDensities = true;
        emit recordingDivergence("Inserting densities...");
    }
    if(e->key() == (int)'v' || e->key() == (int)'V') {
        insertVelocities = true;
        emit recordingDivergence("Inserting velocities...");
    }
    if(e->key() == (int)'h' || e->key() == (int)'H') {
        insertHeat = true;
        emit recordingDivergence("Inserting heat...");
    }
    if(e->key() == (int)'b' || e->key() == (int)'B') {
        boom_trigger = true;
        emit recordingDivergence("Inserting charges...");
    }
    if(e->key() == (int)'a' || e->key() == (int)'A') {
        addDensities = true;
        emit recordingDivergence("Adding densities...");
    }
    if(e->key() == (int)'s' || e->key() == (int)'S') {
        subtractDensities = true;
        emit recordingDivergence("Subtracting densities...");
    }
    if(e->key() == (int)'i' || e->key() == (int)'I') {
        ignite = true;
        emit recordingDivergence("Igniting...");
    }
    if(e->key() == (int)'6') {
        solver->sum_divergence = true;
        solver->iterations = 0;
    }
    if(e->key() == (int)'7') {
        solver->sum_divergence = false;
    }

    if(e->key() == Qt::Key_Space && solver->simEnabled) {
        solver->simEnabled = false;
        emit simState("Simulation: PAUSED");
    }
    else if(e->key() == Qt::Key_Space && !solver->simEnabled) {
        solver->simEnabled = true;
        emit simState("Simulation: ON");
    }
    if(e->key() == Qt::Key_Right) solver->oneStep = true;

    if(e->key() == (int)'1') solver->vel_field(1);
    if(e->key() == (int)'2') solver->vel_field(2);
    if(e->key() == (int)'3') solver->vel_field(3);
    if(e->key() == (int)'4') solver->vel_field(4);
    if(e->key() == (int)'5') solver->vel_field(5);
}

void OpenGLPanel::keyReleaseEvent(QKeyEvent * e) {
    key_pressed = false;
    insertDensities = false;
    addDensities = false;
    subtractDensities = false;
    insertVelocities = false;
    insertHeat = false;
    ignite = false;
    boom_trigger = false;
    emit recordingDivergence("");
}

void OpenGLPanel::set_dt(int new_dt){
    solver->set_dt(new_dt/100.0);
}

void OpenGLPanel::set_diff(int new_diff){
    solver->set_diff(new_diff/1000000.0);
}

void OpenGLPanel::set_heat_diff(int new_diff){
    solver->set_heat_diff(new_diff/1000000.0);
}

void OpenGLPanel::set_visc(int new_visc){
    solver->set_visc(new_visc/500000.0);
}


void OpenGLPanel::diffusionTest(){
    reset();
    solver->diffusion_test();
}

void OpenGLPanel::advectionTest(){
    reset();
    solver->advection_test();
}

void OpenGLPanel::divergenceTest() {

//    int rad = DIM/8;
//    for(int i = DIM/2-rad; i<DIM/2+rad; i++) {
//        for(int j = DIM/2-rad; j<DIM/2+rad; j++) {
//            solver->div_constants[IX(i,j)] = .004;
//        }
//    }
    int x, center_x = DIM/2 - 5;
    int y, center_y = DIM/2 - 5;

    int r;
    float theta;

    for(r=0; r<5; r++){
        for(theta=0; theta<2*PI; theta += PI/50){
            x = center_x + r*cos(theta);
            y = center_y + r*sin(theta);
            if (x >=0 && x <= DIM+1 && y >= 0 && y <= DIM+1) {
                solver->div_constants[IX(y,x)] = .004;
                //solver->d[IX(y,x)] = 1;
            }
        }
    }
}

void OpenGLPanel::cigTest(){
    reset();
    cigTestOn = true;
    solver->cig_test();
}

void OpenGLPanel::reset(){
    cigTestOn = false;
    buoyancyColoring = false;
    explosions.clear();
    solver->reset();
}

void OpenGLPanel::setVelocities(bool enabled){
    showVelocities = enabled;
}

void OpenGLPanel::setDensities(bool enabled){
    showDensities = enabled;
}

void OpenGLPanel::setBadDiffusion(bool enabled){
    solver->bad_diffusion = enabled;
}

void OpenGLPanel::setSimEnable(bool enabled){
    solver->simEnabled = enabled;
}

void OpenGLPanel::setScalars(bool enabled){
    showScalars = enabled;
}

void OpenGLPanel::setParticles(bool enabled){
    showParticles = enabled;
}

void OpenGLPanel::scalarField(QString newField){
    if (newField == "Velocities (Horz)"){
        currentScalar = 1;
        solver->scalar_field(1);
    }
    else if (newField == "Velocities (Vert)"){
        currentScalar =2;
        solver->scalar_field(2);
    }
    else if (newField == "Divergence"){
        currentScalar= 3;
        solver->scalar_field(3);
    }
}

void OpenGLPanel::velocityField(QString newField){
    if (newField == "Field 1"){
        solver->vel_field(1);
    }
    else if (newField == "Field 2"){
        solver->vel_field(2);
    }
    else if (newField == "Field 3"){
        solver->vel_field(3);
    }
    else if (newField == "Field 4"){
        solver->vel_field(4);
    }
    else if (newField == "Field 5"){
        solver->vel_field(5);
    }
}

void OpenGLPanel::printDensities(){
    ofstream dens_file;
    char val[50];
    dens_file.open(DENSITY_FILE);
    for(int i=0; i<DIM+2; i++){
        for(int j=0; j<DIM+2; j++){
            sprintf(val, "%f,", solver->d[IX(i,j)]);
            dens_file << val;
        }
        dens_file << "\n";
    }
    dens_file.close();
}

void OpenGLPanel::printVelocities(){
    ofstream vel_file;
    char val[50];
    vel_file.open(VELOCITY_FILE);
    vel_file << "VERTICAL:\n";
    for(int i=0; i<DIM+2; i++){
        for(int j=0; j<DIM+2; j++){
            sprintf(val, "%f,", solver->u[IX(i,j)]);
            vel_file << val;
        }
        vel_file << "\n";
    }
    vel_file << "\n\nHORZONTAL:\n";
    for(int i=0; i<DIM+2; i++){
        for(int j=0; j<DIM+2; j++){
            sprintf(val, "%f,", solver->v[IX(i,j)]);
            vel_file << val;
        }
        vel_file << "\n";
    }
    vel_file.close();
}

void OpenGLPanel::setMethod(QString new_method) {

    if("Gauss-Seidel" == new_method) {
        solver->method = GAUSS_SEIDEL;
    }
    else if("Conjugate Gradient" == new_method) {
        solver->method = CONJ_GRAD;
    }
}

void OpenGLPanel::setBFECC(bool enabled){
    solver->bfecc = enabled;
}

void OpenGLPanel::startRecord(bool enabled){
    currentFrame = 0;
    record = !enabled;
}

void OpenGLPanel::stopRecord(bool disabled){
    currentFrame = 0;
    record = disabled;
}

void OpenGLPanel::insertParticles(){

    reset();

    int row_num = 200, col_num = 200;
    solver->num_particles = row_num*col_num;

    float start_x = this->width()/2.0 - this->width()/10.0;
    float start_y = this->height()/2.0 - this->height()/10.0;
    float x_increment = (this->width()/(5.0))/row_num;
    float y_increment = (this->height()/(5.0))/col_num;

    for(int i=0; i<row_num; i++){
        for(int j=0; j<col_num; j++){
            solver->particle_pos_X[i*col_num+j] = (start_x + i*x_increment)/this->width();
            solver->particle_pos_Y[i*col_num+j] = (start_y + j*y_increment)/this->height();
            solver->particle_life[i*col_num+j] = 1;
            solver->particle_temp[i*col_num+j] = 0;
            red[i*col_num+j] = rand()/(1.0*INT_MAX);
            green[i*col_num+j] = rand()/(1.0*INT_MAX);
            blue[i*col_num+j] = rand()/(1.0*INT_MAX);
        }
    }
}


void OpenGLPanel::insertRandomParticles(){

    reset();

    solver->num_particles = 100000;

    for(int i=0; i<solver->num_particles; i++){
            solver->particle_pos_X[i] = (rand()/(1.0*INT_MAX))*(1-2*solver->h_space) + solver->h_space;
            solver->particle_pos_Y[i] = (rand()/(1.0*INT_MAX))*(1-2*solver->h_space) + solver->h_space;
            solver->particle_life[i] = 1;
            solver->particle_temp[i] = 0;
            red[i] = rand()/(1.0*INT_MAX);
            green[i] = rand()/(1.0*INT_MAX);
            blue[i] = rand()/(1.0*INT_MAX);
    }
    for(int i=0; i<(DIM+2)/2; i++){
        for(int j=0; j<DIM+2; j++){
            solver->t[IX(j,i)] = .7;
        }
    }
}

void OpenGLPanel::buoyancyTest(){
    reset();
    buoyancyColoring = true;
    solver->buoyancy_test();
}

void OpenGLPanel::set_blast_rad(int radius){
    solver->set_blast_rad(radius);
}

void OpenGLPanel::set_blast_intensity(int intensity){
    solver->set_blast_intensity(intensity/1000.0);
}

void OpenGLPanel::buoyancyEnable(bool enabled){
    solver->buoyancyEnable(enabled);
}

void OpenGLPanel::set_buoy(int new_buoy){
    solver->set_buoy(new_buoy/100.0);
}

void OpenGLPanel::set_temp_diff_enable(bool enable){
    solver->temp_diff_enable = enable;
}

void OpenGLPanel::set_temp_adv_enable(bool enable){
    solver->temp_adv_enable = enable;
}

void OpenGLPanel::set_temp_env_enable(bool enable){
    solver->temp_env_enable = enable;
}

void OpenGLPanel::set_blast_shape(QString new_shape) {
    if(new_shape == "Circle") solver->blast_shape = CIRCLE;
    if(new_shape == "Diamond") solver->blast_shape = DIAMOND;
    if(new_shape == "Square") solver->blast_shape = SQUARE;
}

void OpenGLPanel::temperatureTest() {

    reset();
    solver->reset();

//    // A square that is one third of the window width/height
//    for(i=(int)((DIM+2)/3.0); i<=(int)(2*(DIM+2)/3.0); i++){
//        for(j=(int)((DIM+2)/3.0); j<=(int)(2*(DIM+2)/3.0); j++){
//            solver->t[IX(j,i)] = 1;
//        }
//    }

    solver->num_particles = 100000;
    float dx, dy, deg, rad, rad_max = .2;
    for(int i=0; i<solver->num_particles; i++){
//            rad = (-2*rad_max*(rand()/(1.0*INT_MAX))+rad_max);
//            deg = 2*PI*(rand()/(1.0*INT_MAX));
//            dx = (float)(rad*cos(deg));
//            dy = (float)(rad*sin(deg));
//            solver->particle_pos_X[i] = .5+dx;
//            solver->particle_pos_Y[i] = .5+dy;
//            solver->particle_life[i] = 1;
            solver->particle_pos_X[i] = (float)(rand()/(1.0*INT_MAX));
            solver->particle_pos_Y[i] = (float)(rand()/(1.0*INT_MAX));
            solver->particle_life[i] = 1;
            //solver->particle_temp[i] = .8;
    }

    solver->particle_pos_X[solver->num_particles-1] = 0.5+solver->h_space/10;
    solver->particle_pos_Y[solver->num_particles-1] = 0.5+solver->h_space/10;
    //solver->particle_temp[solver->num_particles-1] = 1;

}

void OpenGLPanel::setTemperatures(bool enabled){
    showTemperatures = enabled;
}

void OpenGLPanel::setExplodingParticles(bool enabled){
    solver->setExplodingParticles(enabled);
}
