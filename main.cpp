#include "mainwindow.h"
#include <QApplication>

#include "logging.h"

int main(int argc, char *argv[])
{
    //google::InitGoogleLogging(argv[0]);

    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
