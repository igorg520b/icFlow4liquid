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
    void Prepare() override;        // invoked once, at simulation start
    bool Step() override;           // either invoked by Worker or via GUI
    void RequestAbort() override;   // invoked from GUI
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

    int currentStep;
    double timeStepFactor, simulationTime;

    Model();
    ~Model();
    void Reset(SimParams &prms);
    void InitialGuess(SimParams &prms, double timeStep, double timeStepFactor);
    bool AssembleAndSolve(SimParams &prms, double timeStep, bool restShape = false);  // return true if solved
    bool AcceptTentativeValues(double timeStep);    // return true if plastic deformation occurred
    void GetNewMaterialPosition();

    void UnsafeUpdateGeometry();
    void PositionIndenter(double offset);

    void AttachSpring(double X, double Y, double radius);   // attach spring to nodes
    void ReleaseSpring();
    void AdjustSpring(double dX, double dY);

signals:
//    void requestGeometryUpdate(); // request the main thread to redraw

    // VTK visualization
public:
    enum VisOpt { none, elem_area, energy_density, stress_xx, stress_yy, stress_hydrostatic,
                ps1, ps2, shear_stress, volume_change, velocity_div, elem_group, node_group,
                vel_mag, Green_strain_xx, Green_strain_yy, Green_strain_xy, plasticity_norm,
                adj_elems_count};
    Q_ENUM(VisOpt)
    void ChangeVisualizationOption(icy::Model::VisOpt option);
private:
    QMutex vtk_update_mutex; // to prevent modifying mesh data while updating VTK representation
    bool vtk_update_requested = false;  // true when signal has been already emitted to update vtk geometry

    Eigen::Vector2d spring;

};

#endif // MESHCOLLECTION_H
