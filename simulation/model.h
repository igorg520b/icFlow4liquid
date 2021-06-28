#ifndef MESHCOLLECTION_H
#define MESHCOLLECTION_H

#include <QFileInfo>
#include <QObject>
#include <QMutex>

#include <vector>
#include <algorithm>
#include <chrono>
#include <unordered_set>

#include "parameters_sim.h"
#include "equationofmotionsolver.h"
#include "modelcontrollerinterface.h"

#include <Eigen/Core>

namespace icy { class Model; class Mesh; class Node; class Element; }

class icy::Model : public QObject, public ModelControllerInterface
{
    Q_OBJECT

    // ModelController
public:
    void Prepare(void) override;        // invoked once, at simulation start
    bool Step(void) override;           // either invoked by Worker or via GUI
    void RequestAbort(void) override;   // invoked from GUI
private:
    bool abortRequested = false;
    void Aborting();       // called before exiting Step() if aborted
signals:
    void stepCompleted();
    void stepAborted();

    // Model
public:
    SimParams prms;
    icy::Mesh *mesh;
    EquationOfMotionSolver eqOfMotion;

    double avgSeparationDistance = -1;
    int currentStep = 0;
    double timeStepFactor = 1;

    Model();
    ~Model();
    void Reset(SimParams &prms);
    void InitialGuess(SimParams &prms, double timeStep, double timeStepFactor);
    bool AssembleAndSolve(SimParams &prms, double timeStep);    // true if solved; false -> parameter is time step mult
    void AcceptTentativeValues(double timeStep);
    void UnsafeUpdateGeometry();
    void PositionIndenter(double offset);

signals:
    void requestGeometryUpdate(); // UnsafeUpdateGeometry() is invoked from the main thread

    // VTK visualization
public:
    enum VisOpt { none, elem_area, energy_density, stress_xx, stress_yy, stress_hydrostatic, non_symm_measure,
                ps1, ps2, shear_stress, volume_change, velocity_div, velocity_div_nd};
    Q_ENUM(VisOpt)
    void ChangeVisualizationOption(icy::Model::VisOpt option);
private:
    QMutex vtk_update_mutex; // to prevent modifying mesh data while updating VTK representation
    bool vtk_update_requested = false;  // true when signal has been already emitted to update vtk geometry

};

#endif // MESHCOLLECTION_H
