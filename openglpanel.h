#ifndef OPENGLPANEL_H
#define OPENGLPANEL_H

#include <QGLWidget>
#include <QtOpenGL>
#include <QImage>
#include <glu.h>
#include <cmath>
#include "fluidsolver.h"
#include <list>

using namespace std;

typedef struct {
    float x;
    float y;
    float radius;
    Blast_Shape shape;
} BoomPosition;

class OpenGLPanel : public QGLWidget
{
    Q_OBJECT

private:
    QTimer *timer;
    QImage frame;
    bool showGrid;
    bool showVelocities;
    bool showDensities;
    bool showParticles;
    bool showTemperatures;
    bool showScalars;
    bool simEnabled;
    bool record;
    bool cigTestOn;
    bool insertDensities;
    bool addDensities;
    bool subtractDensities;
    bool insertVelocities;
    bool insertHeat;
    bool ignite;
    bool key_pressed;
    bool boom_trigger;
    bool buoyancyColoring;
    bool boom;
    bool mouseDown;
    float du;
    float dv;

    float * sources;

    float mouse_x, mouse_y;

    list<BoomPosition> explosions;

    void drawGrid(int rows, int cols);

public:
    explicit OpenGLPanel(QWidget *parent = 0);
    ~OpenGLPanel();

    FluidSolver * solver;

protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
    void mouseMoveEvent(QMouseEvent *);
    void mousePressEvent(QMouseEvent *);
    void mouseReleaseEvent(QMouseEvent *);
    void keyPressEvent(QKeyEvent *);
    void keyReleaseEvent(QKeyEvent *);

public slots:
    void setGrid(bool);
    void setVelocities(bool);
    void setDensities(bool);
    void setScalars(bool);
    void setTemperatures(bool);
    void setParticles(bool);
    void scalarField(QString);
    void velocityField(QString);
    void setMethod(QString);
    void set_dt(int);
    void set_diff(int);
    void set_visc(int);
    void set_blast_rad(int);
    void set_blast_shape(QString);
    void set_blast_intensity(int);
    void diffusionTest();
    void advectionTest();
    void divergenceTest();
    void temperatureTest();
    void buoyancyTest();
    void insertParticles();
    void insertRandomParticles();
    void cigTest();
    void reset();
    void setBadDiffusion(bool);
    void setSimEnable(bool);
    void setBFECC(bool);
    void printDensities();
    void printVelocities();
    void startRecord(bool);
    void stopRecord(bool);
    void buoyancyEnable(bool);
    void set_buoy(int);
    void set_heat_diff(int);
    void setExplodingParticles(bool);


    void set_temp_diff_enable(bool);
    void set_temp_adv_enable(bool);
    void set_temp_env_enable(bool);


signals:
    void mouseX(QString);
    void mouseY(QString);
    void simState(QString);
    void cflViolations(QString);
    void recordingDivergence(QString);
    void recordingIterations(QString);

private:
    void allocateMemory();
    void deallocateMemory();

};

#endif // OPENGLPANEL_H
