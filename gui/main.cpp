#include "mainwindow.h"

#include <QApplication>
#include <QSurfaceFormat>
#include <QCommandLineParser>
#include <gmsh.h>
#include <omp.h>
#include <iostream>

//#include "equationofmotionsolver.h"

int main(int argc, char *argv[])
{
    std::cout << "num_threads " << omp_get_max_threads() << std::endl;
    std::cout << "testing threads" << std::endl;
    int nthreads, tid;
#pragma omp parallel
    { std::cout << omp_get_thread_num(); }
    std::cout << std::endl;

    //EquationOfMotionSolver eq; eq.TestSolve(); return 0;

    gmsh::initialize();

    QApplication a(argc, argv);
    QApplication::setApplicationName("icFlow");
    QApplication::setApplicationVersion("3.0");

    // parse command line options
    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addVersionOption();
    // An option with a value
    QCommandLineOption initialConfigOption(QStringList() << "f" << "file",
            QCoreApplication::translate("main", "Load initial configuration from a JSON file <fileName>"),
            QCoreApplication::translate("main", "fileName"));

    parser.addOption(initialConfigOption);
    parser.process(a);

    QSurfaceFormat fmt = QVTKOpenGLNativeWidget::defaultFormat();
    QSurfaceFormat::setDefaultFormat(fmt);

    MainWindow w;
//    w.show();
    w.showMaximized();
    return a.exec();
}
