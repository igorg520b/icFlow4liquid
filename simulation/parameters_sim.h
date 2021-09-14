#ifndef P_SIM_H
#define P_SIM_H

#include <QObject>
#include <QDebug>

#include <Eigen/Core>
#include <iostream>

// variables related to the formulation of the model

namespace icy { class SimParams; }

class icy::SimParams : public QObject
{
    Q_OBJECT

    // general
    Q_PROPERTY(int s_MaxSteps MEMBER MaxSteps NOTIFY propertyChanged)
    // meshing
    Q_PROPERTY(double s_ElemSize MEMBER CharacteristicLength NOTIFY propertyChanged)

    // integration
    Q_PROPERTY(double in_InitialTimeStep MEMBER InitialTimeStep NOTIFY propertyChanged)
    Q_PROPERTY(double in_ConvergenceEpsilon MEMBER ConvergenceEpsilon NOTIFY propertyChanged)
    Q_PROPERTY(double in_ConvergenceCutoff MEMBER ConvergenceCutoff NOTIFY propertyChanged)
    Q_PROPERTY(int in_MinIter MEMBER MinIter NOTIFY propertyChanged)
    Q_PROPERTY(int in_MaxIter MEMBER MaxIter NOTIFY propertyChanged)

    // material parameters and physical constants
    Q_PROPERTY(double p_Gravity MEMBER Gravity NOTIFY propertyChanged)
    Q_PROPERTY(double p_Density MEMBER Density NOTIFY propertyChanged)
    Q_PROPERTY(double p_YoungsModulus READ getYoungsModulus WRITE setYoungsModulus)
    Q_PROPERTY(double p_PoissonsRatio READ getPoissonsRatio WRITE setPoissonsRatio)
    Q_PROPERTY(double p_Kappa READ getKappa)
    Q_PROPERTY(double p_Thickness MEMBER Thickness NOTIFY propertyChanged)

    Q_PROPERTY(double p_PlasticYieldThreshold MEMBER PlasticYieldThreshold NOTIFY propertyChanged)
    Q_PROPERTY(double p_PlasticFlowRate MEMBER PlasticFlowRate NOTIFY propertyChanged)
    Q_PROPERTY(bool f_EnablePlasticity MEMBER EnablePlasticity NOTIFY propertyChanged)

    Q_PROPERTY(double p_InteractionDistance MEMBER InteractionDistance NOTIFY propertyChanged)
    Q_PROPERTY(bool f_EnableCollisions MEMBER EnableCollisions NOTIFY propertyChanged)

    // fracture
    Q_PROPERTY(double f_WeakeningCoeff MEMBER FractureWeakeningCoeff NOTIFY propertyChanged)
    Q_PROPERTY(bool f_EnableFracture MEMBER EnableFracture NOTIFY propertyChanged)
    Q_PROPERTY(double f_TemporalAttenuation MEMBER FractureTemporalAttenuation NOTIFY propertyChanged)
    Q_PROPERTY(int f_MaxSubsteps MEMBER FractureMaxSubsteps NOTIFY propertyChanged)
    Q_PROPERTY(double f_TractionThreshold MEMBER FractureTractionThreshold NOTIFY propertyChanged)
    Q_PROPERTY(unsigned f_SubstepLevels MEMBER FractureSubstepLevels NOTIFY propertyChanged)
    Q_PROPERTY(unsigned f_TimerLevels MEMBER FractureTimerLevels NOTIFY propertyChanged)
    Q_PROPERTY(double f_AreaThreshold MEMBER FractureAreaThreshold NOTIFY propertyChanged)
    Q_PROPERTY(double f_AngleThreshold MEMBER FractureAngleThreshold NOTIFY propertyChanged)

public:
    int MaxSteps, MinIter, MaxIter;
    double InitialTimeStep;
    double Gravity, Density, PoissonsRatio, YoungsModulus, Thickness;
    double CharacteristicLength;

    double ConvergenceEpsilon, ConvergenceCutoff;

    // collisions
    double InteractionDistance;
    bool EnableCollisions;

    // plasticity
    double PlasticYieldThreshold, PlasticFlowRate;
    bool EnablePlasticity;

    // fracture
    double FractureWeakeningCoeff;
    bool EnableFracture;
    double FractureTemporalAttenuation;
    int FractureMaxSubsteps;
    double FractureTractionThreshold;
    unsigned FractureSubstepLevels, FractureTimerLevels;
    double FractureAngleThreshold, FractureAreaThreshold;

    double lambda, mu, Kappa;
    double getKappa() {return Kappa;}

    double getYoungsModulus() {return YoungsModulus;}
    void setYoungsModulus(double ym)
    {
        YoungsModulus=ym;
        RecomputeLamdaMuAndKappa();
    }

    double getPoissonsRatio() {return PoissonsRatio;}
    void setPoissonsRatio(double nu)
    {
        PoissonsRatio=nu;
        RecomputeLamdaMuAndKappa();
    }

    void RecomputeLamdaMuAndKappa()
    {
        lambda = (YoungsModulus*PoissonsRatio)/((1.0+PoissonsRatio)*(1.0-2.0*PoissonsRatio)); // Lamé's first parameter
        mu = YoungsModulus/(2*(1+PoissonsRatio));                 // Lamé's second parameter - shear modulus
        Kappa = 0.1*lambda/CharacteristicLength;
        emit propertyChanged();
    }


    SimParams() { Reset(); }

    void Reset()
    {
        MaxSteps = 20000;
        InitialTimeStep = 0.05;

        // material parameters and physical constants
        Gravity = 9.81;
        Density = 1;
        Thickness = 0.1;

        PoissonsRatio = 0.3;
        YoungsModulus = 5000;

        InteractionDistance = 0.01;

        CharacteristicLength = 0.07;
        ConvergenceEpsilon = 1e-2;
        ConvergenceCutoff = 1e-7;

        MinIter = 3;
        MaxIter = 6;
        RecomputeLamdaMuAndKappa();

        PlasticYieldThreshold = 100;
        PlasticFlowRate = 1;

        // fracture
        FractureWeakeningCoeff = 0.75;
        FractureTemporalAttenuation = 0.025;
        FractureMaxSubsteps = 1000;
        FractureTractionThreshold = 500;
        FractureSubstepLevels = 4;
        FractureTimerLevels = 10;
        FractureAngleThreshold = 10;    // in degrees
        FractureAreaThreshold = 1e-4;

        // features
        EnableCollisions = true;
        EnablePlasticity = false;
        EnableFracture = true;
    }


signals:
    void propertyChanged();
};


#endif
