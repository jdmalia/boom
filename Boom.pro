#-------------------------------------------------
#
# Project created by QtCreator 2013-08-21T13:06:09
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += x86_64

TARGET = Boom
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    openglpanel.cpp \
    fluidsolver.cpp

HEADERS  += mainwindow.h \
    openglpanel.h \
    fluidsolver.h

INCLUDEPATH += "/Users/jason/Downloads/eigen-eigen-ffa86ffb5570/Eigen/"
INCLUDEPATH +=  "/Users/jason/Downloads/glog-0.3.3/src/glog/"
INCLUDEPATH += "/Users/jason/Downloads/boost_1_54_0/"

LIBS += "-L/Users/jason/Downloads/boost_1_54_0/stage/lib/"

FORMS    += mainwindow.ui
