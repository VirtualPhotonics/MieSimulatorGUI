#ifndef TEST_CALCULATE_H
#define TEST_CALCULATE_H

#include <QObject>
#include <QTest>
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
    void test_ComputeMuspAtRefWavel_sanityCheck();

private:
    parameters *mPara;
    calculate *mCalc;
};

#endif // TEST_CALCULATE_H
