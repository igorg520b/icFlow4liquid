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
    void Reset(unsigned setup);
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
    void fractureProgress();    // signals that the VTK view of mesh is not in sync with internal representation

    // Model
public:
    SimParams prms;
    icy::Mesh mesh;
    EquationOfMotionSolver eqOfMotion;

    int currentStep;
    double timeStepFactor, simulationTime;

    void SetIndenterPosition(double position);
    void AttachSpring(double X, double Y, double radius);   // attach spring to nodes
    void ReleaseSpring();
    void AdjustSpring(double dX, double dY);

private:
    void Fracture(double timeStep);
    void Fracture_LocalSubstep();    // part of Fracture()
    void InitialGuess(double timeStep, double timeStepFactor);
    bool AssembleAndSolve(double timeStep, bool enable_collisions, bool enable_spring,
                          std::vector<icy::Node*> &nodes, std::vector<icy::Element*> &elems,
                          std::vector<icy::CohesiveZone*> &czs);  // return true if solved
    bool AcceptTentativeValues(double timeStep);    // return true if plastic deformation occurred
    void GetNewMaterialPosition();      // relax the mesh to a new rest state

    // VTK visualization
public:
    enum VisOpt { none, plasticity_norm, plasticity_gamma, plasticity_tau_ratio, nd_max_normal_traction,
                  nd_isBoundary, elem_isBoundary,
                  stress_xx, stress_yy, stress_hydrostatic, ps1, ps2, shear_stress,
                  Green_strain_xx, Green_strain_yy, Green_strain_xy,
                  quality_measure, avg_edge_len,
                  elem_area, energy_density, volume_change, velocity_div, node_group, elem_group,
                  node_traversal,
                  vel_mag, adj_elems_count_nd};
    Q_ENUM(VisOpt)

    void UnsafeSynchronizeVTK();    // synchronize what VTK shows with internal mesh representation; invoke from the main thread
    void ChangeVisualizationOption(icy::Model::VisOpt option);  // invoke from the main thread
    bool topologyInvalid = true;
    bool displacementsInvalid = true;

private:
    QMutex vtk_update_mutex; // to prevent modifying mesh data while updating VTK representation
    bool vtk_update_requested = false;  // true when signal has been already emitted to update vtk geometry

    Eigen::Vector2d spring;

};

#endif // MESHCOLLECTION_H
