#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->velocityBox->addItem("Field 1");
    ui->velocityBox->addItem("Field 2");
    ui->velocityBox->addItem("Field 3");
    ui->velocityBox->addItem("Field 4");
    ui->velocityBox->addItem("Field 5");

    ui->scalarBox->addItem("Velocities (Horz)");
    ui->scalarBox->addItem("Velocities (Vert)");
    ui->scalarBox->addItem("Divergence");

    ui->methodBox->addItem("Conjugate Gradient");
    ui->methodBox->addItem("Gauss-Seidel");

    ui->blastBox->addItem("Circle");
    ui->blastBox->addItem("Diamond");
    ui->blastBox->addItem("Square");
}

MainWindow::~MainWindow()
{
    delete ui;
}
