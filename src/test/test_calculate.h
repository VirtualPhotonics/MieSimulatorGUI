#ifndef TEST_CALCULATE_H
#define TEST_CALCULATE_H

#include <QObject>
#include "../app/parameters.h"
#include "../app/calc/calculate.h"

class TestCalculate : public QObject
{
    Q_OBJECT

public:
    TestCalculate();
    ~TestCalculate();

private slots:    
    void init();
    void cleanup();

    void test_DoSimulationdata();
    void test_DoSimulation_RealVsComplexRefIndex();
    void test_DoSimulation_MultiWavelength();
    void test_DoSimulation_ZeroDensitySafety();
    void test_DiameterRangeSetting_logNormal();
    void test_DiameterRangeSetting_gaussian();
    void test_DiameterRangeSetting_monoDisperse();
    void test_SetSphereRadiusAndRefIndex_monoDisperse();
    void test_SetSphereRadiusAndRefIndex_poly_logNormal_volfrac();
    void test_SetSphereRadiusAndRefIndex_poly_gaussian_conc();
    void test_CalculateG();
    void test_CalculateForwardBackward();
    void test_CalculatePowerLawAutoFitSimple();
    void test_CalculatePowerLawAutoFitComplex();
    void test_PhaseFunctionNormalization();
    void test_ComputeMuspAtRefWavel_sanityCheck();
    void test_CheckIndependentScattering_HighConcentration();
    void test_CheckIndependentScattering_LowConcSmallSize();
    void test_CheckIndependentScattering_TransitionalDependent();
    void test_CheckIndependentScattering_LowConcLargeSizeIndependent();
    void test_CheckIndependentScattering_Polydisperse();

private:
    Parameters *mPara;
    Calculate *mCalc;
};

#endif // TEST_CALCULATE_H
