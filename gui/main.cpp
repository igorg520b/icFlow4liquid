#include "mainwindow.h"

#include <QApplication>
#include <QSurfaceFormat>
#include <QCommandLineParser>
#include <gmsh.h>
#include "spdlog/spdlog.h"


int main(int argc, char *argv[])
{
#ifndef QT_DEBUG
    spdlog::info("testing threads {}", omp_get_max_threads());
    int nthreads, tid;
#pragma omp parallel
    {     spdlog::info("{}", omp_get_thread_num()); }
#endif
    spdlog::set_pattern("%v");
    gmsh::initialize();

    QApplication a(argc, argv);
    QApplication::setApplicationName("icFlow");
    QApplication::setApplicationVersion("3.0");

    QSurfaceFormat fmt = QVTKOpenGLNativeWidget::defaultFormat();
    QSurfaceFormat::setDefaultFormat(fmt);

    MainWindow w;
    w.showMaximized();
    return a.exec();
}
