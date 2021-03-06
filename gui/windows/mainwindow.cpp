#include <QFileDialog>
#include <algorithm>
#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "mesh.h"

MainWindow::~MainWindow() {delete ui;}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    worker = new BackgroundWorker(&model);

    connect(&model, SIGNAL(fractureProgress()),SLOT(updateGUI()));
    connect(&model, SIGNAL(stepAborted()),SLOT(updateGUI()));
    connect(&model, SIGNAL(stepCompleted()), SLOT(updateGUI()));
    connect(worker, SIGNAL(workerPaused()), SLOT(background_worker_paused()));

    // property browser
    pbrowser = new ObjectPropertyBrowser(this);

    // VTK
    qt_vtk_widget = new QVTKOpenGLNativeWidget();
    qt_vtk_widget->setRenderWindow(renderWindow);

    renderer->SetBackground(1.0,1.0,1.0);
    renderer->AddActor(model.mesh.actor_collisions);
    renderer->AddActor(model.mesh.actor_mesh_deformable);
    renderer->AddActor(model.mesh.actor_boundary_all);
    renderer->AddActor(model.mesh.actor_boundary_intended_indenter);
    renderer->AddActor(model.mesh.actor_czs);

    renderWindow->AddRenderer(renderer);
    renderWindow->GetInteractor()->SetInteractorStyle(specialSelector2D);
    specialSelector2D->mw = this;
    renderer->AddActor(specialSelector2D->actor);


    // splitter
    splitter_main = new QSplitter(Qt::Orientation::Horizontal);
    splitter_main->addWidget(pbrowser);
    splitter_main->addWidget(qt_vtk_widget);
    //splitter_main->setStretchFactor(0,100);
    //splitter_main->setStretchFactor(1,500);
    splitter_main->setSizes(QList<int>({100, 500}));
    setCentralWidget(splitter_main);

    // toolbar - comboboxes
    comboBox_visualizations = new QComboBox();
    ui->toolBar->addWidget(comboBox_visualizations);

    // populate combobox
    QMetaEnum qme = QMetaEnum::fromType<icy::Model::VisOpt>();
    for(int i=0;i<qme.keyCount();i++) comboBox_visualizations->addItem(qme.key(i));

    connect(comboBox_visualizations, QOverload<int>::of(&QComboBox::currentIndexChanged),
            [=](int index){ comboboxIndexChanged_visualizations(index); });


    // slider
    slider1 = new QSlider(Qt::Horizontal);
    ui->toolBar->addWidget(slider1);
    slider1->setTracking(true);
    slider1->setMinimum(0);
    slider1->setMaximum(1000);
    connect(slider1, SIGNAL(valueChanged(int)), this, SLOT(sliderValueChanged(int)));

    labelStepCount = new QLabel();

    // statusbar
    statusLabel = new QLabel("-");
    statusLabelStepFactor = new QLabel(" --- ");
    labelVolumeChange = new QLabel(" --- ");
    labelElapsedTime = new QLabel(" --- ");

    QSizePolicy sp;
    sp.setHorizontalPolicy(QSizePolicy::Fixed);
    statusLabelStepFactor->setSizePolicy(sp);
    statusLabelStepFactor->setFixedWidth(100);
    labelStepCount->setSizePolicy(sp);
    labelStepCount->setFixedWidth(100);
    labelVolumeChange->setSizePolicy(sp);
    labelVolumeChange->setFixedWidth(100);
    labelElapsedTime->setSizePolicy(sp);
    labelElapsedTime->setFixedWidth(100);

    ui->statusbar->addWidget(statusLabel);
    ui->statusbar->addPermanentWidget(labelElapsedTime);
    ui->statusbar->addPermanentWidget(labelStepCount);
    ui->statusbar->addPermanentWidget(statusLabelStepFactor);
    ui->statusbar->addPermanentWidget(labelVolumeChange);

    // read/restore saved settings
    QSettings settings(m_sSettingsFile);
    splitter_main->restoreState(settings.value("splitter_main").toByteArray());

    vtkCamera* camera = renderer->GetActiveCamera();
    camera->ParallelProjectionOn();

    QVariant var = settings.value("camData");
    if(!var.isNull()) {
        QByteArray arr = var.toByteArray();
        double *vec = (double*)arr.constData();
        camera->SetPosition(vec[0],vec[1],vec[2]);
        camera->SetFocalPoint(vec[3],vec[4],vec[5]);
        camera->SetViewUp(vec[6],vec[7],vec[8]);
        camera->SetViewAngle(vec[9]);
    }
    camera->Modified();

    renderer->AddActor(scalarBar);
    scalarBar->SetLookupTable(model.mesh.hueLut);

    scalarBar->SetMaximumWidthInPixels(130);
    scalarBar->SetBarRatio(0.07);
    scalarBar->SetMaximumHeightInPixels(400);
    scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
    scalarBar->GetPositionCoordinate()->SetValue(0.90,0.015, 0.0);
    scalarBar->SetLabelFormat("%.1e");
    scalarBar->GetLabelTextProperty()->BoldOff();
    scalarBar->GetLabelTextProperty()->ItalicOff();
    scalarBar->GetLabelTextProperty()->ShadowOff();
    scalarBar->GetLabelTextProperty()->SetColor(0.1,0.1,0.1);

    renderWindow->Render();

    pbrowser->setActiveObject(&model.prms);

    model.Reset(settings.value("setup_option",0).toInt());
    comboBox_visualizations->setCurrentIndex(settings.value("vis_option").toInt());

    slider1->setValue(310);
    updateGUI();
}

void MainWindow::showEvent( QShowEvent*)
{
}


void MainWindow::on_action_quit_triggered() { this->close(); }

void MainWindow::closeEvent( QCloseEvent* event )
{
    // save settings and stop simulation
    QSettings settings(m_sSettingsFile);
    qDebug() << "MainWindow: closing main window; " << settings.fileName();

    settings.setValue("splitter_main", splitter_main->saveState());

    double data[10];
    renderer->GetActiveCamera()->GetPosition(&data[0]);
    renderer->GetActiveCamera()->GetFocalPoint(&data[3]);
    renderer->GetActiveCamera()->GetViewUp(&data[6]);
    data[9]=renderer->GetActiveCamera()->GetViewAngle();

    QByteArray arr((char*)&data[0], sizeof(double)*10);
    settings.setValue("camData", arr);
    settings.setValue("vis_option", comboBox_visualizations->currentIndex());
    settings.setValue("setup_option", model.mesh.typeOfSetup);

    // kill backgroundworker
    worker->Finalize();
    event->accept();
}

void MainWindow::on_action_simulation_start_triggered(bool checked)
{
    if(!worker->running && checked){
        qDebug() << "start button - starting";
        statusLabel->setText("starting simulation");
        model.Prepare();

        worker->Resume();
    }
    else if(worker->running && !checked)
    {
        statusLabel->setText("pausing simulation");
        qDebug() << "start button - pausing";
        worker->Pause();
        ui->action_simulation_start->setEnabled(false);
    }

}

void MainWindow::background_worker_paused()
{
    // enable the "Start" button
    qDebug() << "MainWindow::background_worker_paused()";
    ui->action_simulation_start->blockSignals(true);
    ui->action_simulation_start->setEnabled(true);
    ui->action_simulation_start->setChecked(false);
    ui->action_simulation_start->blockSignals(false);
    statusLabel->setText("stopped");
}




void MainWindow::on_action_simulation_single_step_triggered()
{
    qDebug() << "take one step";
    model.Prepare();
    model.Step();
}

void MainWindow::on_action_camera_reset_triggered()
{
    qDebug() << "MainWindow::on_action_camera_reset_triggered()";
    vtkCamera* camera = renderer->GetActiveCamera();
    camera->ParallelProjectionOn();
    camera->SetClippingRange(1e1,1e3);
    camera->SetFocalPoint(0.0, 0.0, 0.0);
    camera->SetPosition(0.0, 0.0, 50.0);
    camera->SetViewUp(0.0, 1.0, 0.0);
    camera->Modified();
    renderWindow->Render();
}



void MainWindow::updateGUI()
{
    bool r = worker->running;
    ui->action_simulation_single_step->setEnabled(!r);

    labelStepCount->setText(QString{"step: %1"}.arg(model.currentStep));
    statusLabelStepFactor->setText(QString{"%1"}.arg(model.timeStepFactor, 6, 'f', 3, '0'));
    labelElapsedTime->setText(QString{"%1"}.arg(model.simulationTime, 5, 'f', 3, '0'));
    labelVolumeChange->setText(QString{"%1%"}.arg((model.mesh.area_current/model.mesh.area_initial)*100, 6, 'f', 3, '0'));

    render_results();
}


void MainWindow::render_results()
{
    model.UnsafeSynchronizeVTK();
    renderWindow->Render();
}

void MainWindow::comboboxIndexChanged_visualizations(int index)
{
    qDebug() << "comboboxIndexChanged_visualizations " << index;
    model.ChangeVisualizationOption((icy::Model::VisOpt)index);
    scalarBar->SetVisibility(index != 0);
    renderWindow->Render();
}



void MainWindow::sliderValueChanged(int val)
{
    model.SetIndenterPosition(0.001 * val);
    render_results();
}


void MainWindow::on_actionSave_Mesh_triggered()
{
    /*
    qDebug() << "Save mesh";

    QDir outputDirectory;
    QString outputPathName = "output";
    outputDirectory.setPath(outputPathName);
    if(!outputDirectory.exists()) {
        bool result = QDir::current().mkdir(outputPathName);
        if(!result) throw std::runtime_error("could not create output directory");
        outputDirectory.setPath(QDir::current().path() + "/output");
        if(!outputDirectory.exists()) throw std::runtime_error("could not open output directory");
    }
    qDebug() << "output directory: " << outputDirectory.path();

    QFileDialog qfd(this, "Save Mesh", outputDirectory.path(), "HDF5 files (*.h5)");
    qfd.setDefaultSuffix("h5");

    qfd.setAcceptMode(QFileDialog::AcceptSave);
    qfd.setFileMode(QFileDialog::AnyFile);

    bool result = qfd.exec();
    if(!result) return;

    QString fileName = qfd.selectedFiles()[0];

    if (fileName.isEmpty()) return;
    model.mesh.fragments[1]->SaveFragment(fileName.toStdString());
*/
}





void MainWindow::on_actionRemesh_triggered()
{
}


void MainWindow::on_actionSwap_Buffers_triggered()
{
}


void MainWindow::on_actionClear_Velocity_triggered()
{
    for(icy::Node *nd : model.mesh.allNodes) nd->vn.setZero();
}

void MainWindow::on_actionUse_Initial_State_triggered()
{
    model.mesh.showDeformation = icy::Mesh::ShowDeformationOption::initial;
    model.displacementsInvalid = true;
    render_results();
}

void MainWindow::on_actionCurrent_Space_triggered()
{
    model.mesh.showDeformation = icy::Mesh::ShowDeformationOption::current;
    model.displacementsInvalid = true;
    render_results();
}


void MainWindow::on_actionIndentation_triggered()
{
    model.Reset(0);
    updateGUI();
}

void MainWindow::on_actionShear_triggered()
{
    model.Reset(1);
    updateGUI();
}

void MainWindow::on_actionStretch_triggered()
{
    model.Reset(2);
    updateGUI();
}

void MainWindow::on_actionSelf_collision_triggered()
{
    model.Reset(3);
    updateGUI();
}


void MainWindow::on_actionCZs_triggered()
{
    model.Reset(4);
    updateGUI();
}

void MainWindow::on_actionCZs_fracture_triggered()
{
    model.Reset(5);
    updateGUI();
}

