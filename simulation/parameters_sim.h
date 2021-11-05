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

    // features
    Q_PROPERTY(bool f_EnableFracture MEMBER EnableFracture NOTIFY propertyChanged)
    Q_PROPERTY(bool f_InsertCZs MEMBER EnableInsertCZs NOTIFY propertyChanged)
    Q_PROPERTY(bool f_EnableCollisions MEMBER EnableCollisions NOTIFY propertyChanged)
    Q_PROPERTY(bool f_EnablePlasticity MEMBER EnablePlasticity NOTIFY propertyChanged)
    Q_PROPERTY(bool f_EnableCZs MEMBER EnableCZs NOTIFY propertyChanged)

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

    Q_PROPERTY(double p_InteractionDistance MEMBER InteractionDistance NOTIFY propertyChanged)

    // fracture
    Q_PROPERTY(double f_WeakeningCoeff MEMBER FractureWeakeningCoeff NOTIFY propertyChanged)
    Q_PROPERTY(double f_TemporalAttenuation MEMBER FractureTemporalAttenuation NOTIFY propertyChanged)
    Q_PROPERTY(int f_MaxSubsteps MEMBER FractureMaxSubsteps NOTIFY propertyChanged)
    Q_PROPERTY(double f_TractionThreshold MEMBER FractureTractionThreshold NOTIFY propertyChanged)
    Q_PROPERTY(unsigned f_SubstepLevels MEMBER FractureSubstepLevels NOTIFY propertyChanged)
    Q_PROPERTY(unsigned f_TimerLevels MEMBER FractureTimerLevels NOTIFY propertyChanged)
    Q_PROPERTY(double f_AreaThreshold MEMBER FractureAreaThreshold NOTIFY propertyChanged)
    Q_PROPERTY(double f_AngleThreshold MEMBER FractureAngleThreshold NOTIFY propertyChanged)

    // cohesive zones
    Q_PROPERTY(double cz_alpha WRITE set_cz_alpha READ get_cz_alpha)
    Q_PROPERTY(double cz_beta WRITE set_cz_beta READ get_cz_beta)
    Q_PROPERTY(double cz_lambda_n WRITE set_cz_lambda_n READ get_cz_lambda_n)
    Q_PROPERTY(double cz_lambda_t WRITE set_cz_lambda_t READ get_cz_lambda_t)
    Q_PROPERTY(double cz_phi_n WRITE set_cz_phi_n READ get_cz_phi_n)
    Q_PROPERTY(double cz_phi_t WRITE set_cz_phi_t READ get_cz_phi_t)
    Q_PROPERTY(double cz_sigma_max WRITE set_cz_sigma_max READ get_cz_sigma_max)
    Q_PROPERTY(double cz_tau_max WRITE set_cz_tau_max READ get_cz_tau_max)
    Q_PROPERTY(double cz_deln READ get_cz_del_n)
    Q_PROPERTY(double cz_delt READ get_cz_del_t)


public:
    SimParams() { Reset(); }

    // enable/disable features
    bool EnableCZs;
    bool EnableFracture;
    bool EnablePlasticity;
    bool EnableCollisions;
    bool EnableInsertCZs;

    // integrator / simulation
    int MaxSteps, MinIter, MaxIter;
    double ConvergenceEpsilon, ConvergenceCutoff;
    double InitialTimeStep;
    double CharacteristicLength;


    // general parameters of the material and the integrator
    double Gravity, Density, PoissonsRatio, YoungsModulus, Thickness;
    double lambda, mu;  // Lamé parameters

    // collisions
    double InteractionDistance, Kappa;

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


    // plasticity
    double PlasticYieldThreshold, PlasticFlowRate;


    // fracture
    double FractureWeakeningCoeff;
    double FractureTemporalAttenuation;
    int FractureMaxSubsteps;
    double FractureTractionThreshold;
    unsigned FractureSubstepLevels, FractureTimerLevels;
    double FractureAngleThreshold, FractureAreaThreshold;


    // cz parameters
    double cz_alpha, cz_beta;           // brittle / ductile
    double cz_lambda_n, cz_lambda_t;    // max traction location
    double cz_phi_n, cz_phi_t;          // fracture energy
    double cz_sigma_max, cz_tau_max;    // traction thresholds

    // computed cz variables
    double cz_del_n, cz_del_t;
    double cz_p_m, cz_p_n;
    double cz_pMtn, cz_pMnt; // < phi_t - phi_n >, < phi_n - phi_t >
    double cz_gam_n, cz_gam_t;
    double cz_nThreshold, cz_tThreshold; // CZ peak separation point

    void set_cz_alpha(double value) { cz_alpha=value; RecomputeCZParams(); emit propertyChanged(); }
    double get_cz_alpha() {return cz_alpha;}

    void set_cz_beta(double value) { cz_beta=value; RecomputeCZParams(); emit propertyChanged(); }
    double get_cz_beta() {return cz_beta;}

    void set_cz_lambda_n(double value) {cz_lambda_n = value; RecomputeCZParams(); emit propertyChanged(); }
    double get_cz_lambda_n() {return cz_lambda_n;}

    void set_cz_lambda_t(double value) {cz_lambda_t = value; RecomputeCZParams(); emit propertyChanged(); }
    double get_cz_lambda_t() {return cz_lambda_t;}

    void set_cz_phi_n(double value) {cz_phi_n = value; RecomputeCZParams(); emit propertyChanged(); }
    double get_cz_phi_n() {return cz_phi_n;}
    void set_cz_phi_t(double value) {cz_phi_t = value; RecomputeCZParams(); emit propertyChanged(); }
    double get_cz_phi_t() {return cz_phi_t;}

    void set_cz_sigma_max(double value) {cz_sigma_max = value; RecomputeCZParams(); emit propertyChanged(); }
    double get_cz_sigma_max() {return cz_sigma_max;}
    void set_cz_tau_max(double value) {cz_tau_max = value; RecomputeCZParams(); emit propertyChanged(); }
    double get_cz_tau_max() {return cz_tau_max;}
    double get_cz_del_n() {return cz_del_n;}
    double get_cz_del_t() {return cz_del_t;}


    void Reset()
    {
        // features
        EnableCollisions = true;
        EnablePlasticity = false; //true;
        EnableFracture = true;
        EnableInsertCZs = true;
        EnableCZs = true;

        MaxSteps = 20000;
        InitialTimeStep = 0.0005;//0.05;

        // material parameters and physical constants
        Gravity = 9.81;
        Density = 1;
        Thickness = 0.1;

        PoissonsRatio = 0.3;
        YoungsModulus = 3e9; // 5000;

        InteractionDistance = 0.01;

        CharacteristicLength = 0.035;//0.07;
        ConvergenceEpsilon = 1e-2;
        ConvergenceCutoff = 1e-7;

        MinIter = 3;
        MaxIter = 6;
        RecomputeLamdaMuAndKappa();

        PlasticYieldThreshold = 100;
        PlasticFlowRate = 1;

        // fracture
        FractureWeakeningCoeff = 0.75;
        FractureTemporalAttenuation = InitialTimeStep*2;//0.025;
        FractureMaxSubsteps = 1000;
        FractureTractionThreshold = 1e6;//1500;
        FractureSubstepLevels = 4;
        FractureTimerLevels = 10;
        FractureAngleThreshold = 10;    // in degrees
        FractureAreaThreshold = 1e-4;

        // cohesive zones
        cz_alpha = 4;
        cz_beta = 4;
        cz_lambda_n = 0.02;
        cz_lambda_t = 0.02;
        cz_phi_n = 30;
        cz_phi_t = 30; // fracture energy
        cz_sigma_max = 5e5;//500;
        cz_tau_max = 5e5;//500;
        RecomputeCZParams();
    }



    double Macaulay(double a, double b) { if (a > b) return a - b; else return 0; }

    void RecomputeCZParams()
    {
        cz_pMnt = Macaulay(cz_phi_n, cz_phi_t);
        cz_pMtn = Macaulay(cz_phi_t, cz_phi_n);

        double rn_sq = cz_lambda_n * cz_lambda_n;
        double rt_sq = cz_lambda_t * cz_lambda_t;
        cz_p_m = (cz_alpha * (cz_alpha - 1.0) * rn_sq) / (1.0 - cz_alpha * rn_sq);
        cz_p_n = (cz_beta * (cz_beta - 1.0) * rt_sq) / (1.0 - cz_beta * rt_sq);

        if (cz_phi_n < cz_phi_t)
        {
            cz_gam_n = pow(cz_alpha / cz_p_m, cz_p_m);
            cz_gam_t = -cz_phi_t * pow(cz_beta / cz_p_n, cz_p_n);
        }
        else
        {
            cz_gam_n = -cz_phi_n * pow(cz_alpha / cz_p_m, cz_p_m);
            cz_gam_t = pow(cz_beta / cz_p_n, cz_p_n);
        }

        cz_del_n = (cz_phi_n / cz_sigma_max) * cz_alpha * cz_lambda_n *
                   pow((1.0 - cz_lambda_n), (cz_alpha - 1.0)) * ((cz_alpha / cz_p_m) + 1.0) *
                   pow(((cz_alpha / cz_p_m) * cz_lambda_n + 1.0), (cz_p_m - 1.0));
        cz_del_t = (cz_phi_t / cz_tau_max) * cz_beta * cz_lambda_t *
                   pow((1.0 - cz_lambda_t), (cz_beta - 1.0)) * ((cz_beta / cz_p_n) + 1.0) *
                   pow(((cz_beta / cz_p_n) * cz_lambda_t + 1.0), (cz_p_n - 1.0));

        cz_nThreshold = cz_del_n * cz_lambda_n;
        cz_tThreshold = cz_del_t * cz_lambda_t;
    }

signals:
    void propertyChanged();
};


#endif
