#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <QFileDialog>
#include <QSizePolicy>
#include <QPushButton>
#include <QSplitter>
#include <QLabel>
#include <QVBoxLayout>
#include <QTreeWidget>
#include <QProgressBar>
#include <QMenu>
#include <QList>
#include <QDebug>
#include <QComboBox>
#include <QMetaEnum>
#include <QDir>

#include <QtCharts>
#include <QVTKOpenGLNativeWidget.h>

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkDataSetMapper.h>
#include <vtkCamera.h>

#include <vtkOBJReader.h>
#include <vtkNamedColors.h>
#include <vtkProperty.h>
#include <vtkVersion.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCompositePolyDataMapper2.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkAxesActor.h>

#include <vtkDataSetSurfaceFilter.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>

#include <vtkAreaPicker.h>
#include <vtkPointPicker.h>
#include <vtkProp3DCollection.h>
#include <vtkRenderWindowInteractor.h>
//#include <vtkCallbackCommand.h>
//#include <vtkInteractorStyleRubberBandPick.h>

#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkJPEGWriter.h>

//#include <vtkInteractorStyleRubberBand2D.h>
#include "SpecialSelector2D.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <cstdint>

#include "objectpropertybrowser.h"
//#include "preferences_gui.h"

#include "model.h"
#include "backgroundworker.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT
private:
    Ui::MainWindow *ui;

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void showEvent( QShowEvent* event ) override;
    void closeEvent( QCloseEvent* event ) override;




private slots:
    void sliderValueChanged(int val);
    void comboboxIndexChanged_visualizations(int index);

    void on_action_quit_triggered();
    void on_action_simulation_single_step_triggered();
    void on_action_simulation_start_triggered(bool checked);
    void on_action_camera_reset_triggered();

    void updateGUI();   // when simulation is started/stopped or when a step is advanced
    void render_results();

    void background_worker_paused();

    void on_actionSetup_2_triggered();

private:
//    PreferencesGUI prefsGUI;
    icy::Model model;
    BackgroundWorker *worker;

    QString m_sSettingsFile = "ic4_config";
    QLabel *statusLabel;                    // statusbar
    QLabel *statusLabelStepFactor;             // attempt # at advancing a time step
    QLabel *labelStepCount;
    QLabel *labelAvgSeparationDistance;
    QComboBox *comboBox_visualizations;
    QSlider *slider1;

    // splitter and the right side window
    QSplitter *splitter_main;

    // properties
    ObjectPropertyBrowser *pbrowser;

    // VTK
    vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;

    QVTKOpenGLNativeWidget *qt_vtk_widget;

    vtkNew<vtkRenderer> renderer;
    vtkNew<SpecialSelector2D> specialSelector2D;

    vtkNew<vtkScalarBarActor> scalarBar;

    const QString wtitle = "icFlow4: Finite Element Simulation";

    friend class SpecialSelector2D;


};
#endif // MAINWINDOW_H
