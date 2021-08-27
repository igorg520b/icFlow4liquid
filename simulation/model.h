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
#include "mesh.h"

#include <Eigen/Core>

namespace icy { class Model; }

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
    constexpr static int colWidth = 12;    // table column width when logging
signals:
    void stepCompleted();
    void stepAborted();

    // Model
public:
    SimParams prms;
    icy::Mesh mesh;
    EquationOfMotionSolver eqOfMotion;

    int currentStep;
    double timeStepFactor, simulationTime;

    void Reset(unsigned setup);
    void InitialGuess(double timeStep, double timeStepFactor);
    bool AssembleAndSolve(double timeStep, bool restShape = false);  // return true if solved
    bool AcceptTentativeValues(double timeStep);    // return true if plastic deformation occurred
    void GetNewMaterialPosition();      // relax the mesh to a new rest state

    void UnsafeUpdateGeometry();
    void SetIndenterPosition(double position);

    void AttachSpring(double X, double Y, double radius);   // attach spring to nodes
    void ReleaseSpring();
    void AdjustSpring(double dX, double dY);

    void Fracture(double timeStep);

signals:
//    void requestGeometryUpdate(); // request the main thread to redraw

    // VTK visualization
public:
    enum VisOpt { none, plasticity_norm, plasticity_gamma, plasticity_tau_ratio,
                  stress_xx, stress_yy, stress_hydrostatic, ps1, ps2, shear_stress,
                  Green_strain_xx, Green_strain_yy, Green_strain_xy,
                  QM1, avg_edge_len,

                  elem_area, energy_density, volume_change, velocity_div, elem_group, node_group,
                vel_mag,
                adj_elems_count_nd,  nd_max_normal_traction};
    Q_ENUM(VisOpt)
    void ChangeVisualizationOption(icy::Model::VisOpt option);
private:
    QMutex vtk_update_mutex; // to prevent modifying mesh data while updating VTK representation
    bool vtk_update_requested = false;  // true when signal has been already emitted to update vtk geometry

    Eigen::Vector2d spring;

};

#endif // MESHCOLLECTION_H
