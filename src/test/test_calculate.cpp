#include <QDebug>
#include <QtMath>
#include <complex>

#include "test_calculate.h"

QTEST_MAIN(TestCalculate)

TestCalculate::TestCalculate()
{
}

TestCalculate::~TestCalculate()
{
}

//Initialize values for tests
void TestCalculate::init()
{
    mPara = new parameters();
    mCalc = new calculate();

    mPara->radArray = new double[mPara->nRadius];
    mPara->numDensityArray = new double[mPara->nRadius];
    mPara->scatRefRealArray = new double[mPara->nRadius];
    mPara->scatRefImagArray = new double[mPara->nRadius];
    mPara->wavelArray = new double[mPara->nWavel];
    mPara->mus = new double[mPara->nWavel];
    mPara->g = new double[mPara->nWavel];
    mPara->cSca = new double[mPara->nWavel];
    mPara->cExt = new double[mPara->nWavel];
    mPara->cBack = new double[mPara->nWavel];
    mPara->SizePara = new double[mPara->nWavel];
    mPara->forward = new double[mPara->nWavel];
    mPara->backward = new double[mPara->nWavel];

    mPara->S1 = new std::complex<double>*[mPara->nWavel];
    mPara->S2 = new std::complex<double>*[mPara->nWavel];
    mPara->phaseFunctionAve = new double*[mPara->nWavel];
    mPara->phaseFunctionPara = new double*[mPara->nWavel];
    mPara->phaseFunctionPerp = new double*[mPara->nWavel];

    for (unsigned int w = 0; w < mPara->nWavel; ++w)
    {
        mPara->S1[w] = new std::complex<double>[mPara->nTheta];
        mPara->S2[w] = new std::complex<double>[mPara->nTheta];
        mPara->phaseFunctionAve[w] = new double[mPara->nTheta];
        mPara->phaseFunctionPara[w] = new double[mPara->nTheta];
        mPara->phaseFunctionPerp[w] = new double[mPara->nTheta];
    }
}

// clean up resources
void TestCalculate::cleanup()
{
    if (mPara->S1) {
        for (unsigned int w = 0; w < mPara->nWavel; ++w)
            delete[] mPara->S1[w];
        delete[] mPara->S1;
        mPara->S1 = nullptr;
    }
    if (mPara->S2) {
        for (unsigned int w = 0; w < mPara->nWavel; ++w)
            delete[] mPara->S2[w];
        delete[] mPara->S2;
        mPara->S2 = nullptr;
    }
    if (mPara->phaseFunctionAve) {
        for (unsigned int w = 0; w < mPara->nWavel; ++w)
            delete[] mPara->phaseFunctionAve[w];
        delete[] mPara->phaseFunctionAve;
        mPara->phaseFunctionAve = nullptr;
    }
    if (mPara->phaseFunctionPara) {
        for (unsigned int w = 0; w < mPara->nWavel; ++w)
            delete[] mPara->phaseFunctionPara[w];
        delete[] mPara->phaseFunctionPara;
        mPara->phaseFunctionPara = nullptr;
    }
    if (mPara->phaseFunctionPerp) {
        for (unsigned int w = 0; w < mPara->nWavel; ++w)
            delete[] mPara->phaseFunctionPerp[w];
        delete[] mPara->phaseFunctionPerp;
        mPara->phaseFunctionPerp = nullptr;
    }

    if (mPara->radArray) delete[] mPara->radArray;
    if (mPara->numDensityArray) delete[] mPara->numDensityArray;
    if (mPara->scatRefRealArray) delete[] mPara->scatRefRealArray;
    if (mPara->scatRefImagArray) delete[] mPara->scatRefImagArray;
    if (mPara->wavelArray) delete[] mPara->wavelArray;
    if (mPara->mus) delete[] mPara->mus;
    if (mPara->g) delete[] mPara->g;
    if (mPara->cSca) delete[] mPara->cSca;
    if (mPara->cExt) delete[] mPara->cExt;
    if (mPara->cBack) delete[] mPara->cBack;
    if (mPara->SizePara) delete[] mPara->SizePara;
    if (mPara->forward) delete[] mPara->forward;
    if (mPara->backward) delete[] mPara->backward;

    mPara->radArray = nullptr;
    mPara->numDensityArray = nullptr;
    mPara->scatRefRealArray = nullptr;
    mPara->scatRefImagArray = nullptr;
    mPara->wavelArray = nullptr;
    mPara->mus = nullptr;
    mPara->g = nullptr;
    mPara->cSca = nullptr;
    mPara->cExt = nullptr;
    mPara->cBack = nullptr;
    mPara->SizePara = nullptr;
    mPara->forward = nullptr;
    mPara->backward = nullptr;

    if (mCalc) delete mCalc;
    if (mPara) delete mPara;
    mCalc = nullptr;
    mPara = nullptr;
}

// Test case for DoSimulationdata
void TestCalculate::test_DoSimulationdata()
{
    // Set up a simple mono-disperse case at 632.8nm and 7
    mPara->nRadius = 1;
    mPara->nTheta = 361;
    mPara->nWavel = 1;
    mPara->meanRadius = 0.5;
    mPara->scatRefReal = 1.5;
    mPara->scatRefImag = -0.5;
    mPara->medRef = 1.0;
    mPara->numDensityArray[0] = 0.1;
    mPara->startWavel = 632.8;
    mPara->endWavel = 632.8;
    mPara->stepWavel = 10;

    mPara->wavelArray[0] = 632.8;

    mCalc->SetSphereRadiusAndRefIndex(mPara, 3, false);
    QLabel mockLabel;
    mCalc->DoSimulation(&mockLabel, mPara);

    double expectedMus1 = 93.3690054109;

    QCOMPARE(mPara->mus[0], expectedMus1);

    // Verify that the g factor is not NaN
    QVERIFY(!std::isnan(mPara->g[0]));

    // Verify that forward/backward scattering is non-zero
    QVERIFY(mPara->forward[0] > 0.0);
    QVERIFY(mPara->backward[0] > 0.0);
}

// Test case for DiameterRangeSetting with log-normal distribution
void TestCalculate::test_DiameterRangeSetting_logNormal()
{
    mPara->meanRadius = 1.0;
    mPara->stdDev = 0.2;
    mPara->nRadius = 31;
    mCalc->DiameterRangeSetting(mPara, 0); // 0 = Log Normal

    // Check range
    QVERIFY(mPara->minRadius > 0);
    QVERIFY(mPara->minRadius < mPara->meanRadius);
    QVERIFY(mPara->maxRadius > mPara->minRadius);   
    QVERIFY(mPara->maxRadius > mPara->meanRadius);
}

// Test case for DiameterRangeSetting with Gaussian distribution
void TestCalculate::test_DiameterRangeSetting_gaussian()
{
    mPara->meanRadius = 1.0;
    mPara->stdDev = 0.2;
    mPara->nRadius = 31;
    mCalc->DiameterRangeSetting(mPara, 1); // 1 = Gaussian

    // Check range
    QVERIFY(mPara->minRadius > 0);
    QVERIFY(mPara->minRadius < mPara->meanRadius);
    QVERIFY(mPara->maxRadius > mPara->minRadius);
    QVERIFY(mPara->maxRadius > mPara->meanRadius);
    // For large radii with small stDev, the difference between max and mean should be
    // roughly equal to mean and min.
    QCOMPARE(mPara->maxRadius - mPara->meanRadius, mPara->meanRadius - mPara->minRadius);
}

void TestCalculate::test_DiameterRangeSetting_monoDisperse()
{
    mPara->meanRadius = 5.0;
    mCalc->DiameterRangeSetting(mPara, 3);
    QCOMPARE(mPara->minRadius, 5.0);
    QCOMPARE(mPara->maxRadius, 5.0);
}

// Test case for SetSphereRadiusAndRefIndex for mono-disperse spheres
void TestCalculate::test_SetSphereRadiusAndRefIndex_monoDisperse()
{
    mPara->nRadius = 1;
    mPara->meanRadius = 1.0;
    mPara->sphNumDensity = 5e6;
    mPara->scatRefReal = 1.5;
    mPara->scatRefImag = 0.1;

    mCalc->SetSphereRadiusAndRefIndex(mPara, 0, true); // The distribution index doesn't matter for nRadius=1.

    //The elements in Array[0] are the same as the initial values
    QCOMPARE(mPara->radArray[0], mPara->meanRadius);
    QCOMPARE(mPara->numDensityArray[0], mPara->sphNumDensity);
    QCOMPARE(mPara->scatRefRealArray[0], mPara->scatRefReal);
    QCOMPARE(mPara->scatRefImagArray[0], mPara->scatRefImag);
}

// Test case for SetSphereRadiusAndRefIndex for poly-disperse with volume fraction
void TestCalculate::test_SetSphereRadiusAndRefIndex_poly_logNormal_volfrac()
{
    mPara->nRadius = 21;
    mPara->meanRadius = 0.5;
    mPara->stdDev = 0.2;
    mPara->volFraction = 0.1;

    // Set min/max radius
    mCalc->DiameterRangeSetting(mPara, 0);
    mCalc->SetSphereRadiusAndRefIndex(mPara, 0, true); // 0 = log-normal, true = use volume fraction.

    // Check if the arrays have been populated.
    QVERIFY(mPara->radArray[0] > 0);
    QVERIFY(mPara->numDensityArray[0] > 0);
    // Check that the sum of the number densities is greater than zero.
    double sumNumDensity = 0.0;
    for (unsigned int i = 0; i < mPara->nRadius; ++i) {
        sumNumDensity += mPara->numDensityArray[i];
    }
    QVERIFY(sumNumDensity > 0.0);
}

// Test case for SetSphereRadiusAndRefIndex for poly-disperse with concentration.
void TestCalculate::test_SetSphereRadiusAndRefIndex_poly_gaussian_conc()
{
    mPara->nRadius = 21;
    mPara->meanRadius = 0.5;
    mPara->stdDev = 0.2;
    mPara->sphNumDensity = 1e8;

    // Set min/max radius
    mCalc->DiameterRangeSetting(mPara, 1);
    mCalc->SetSphereRadiusAndRefIndex(mPara, 1, false); // 1 = Gaussian, false = use concentration.

    // Check if the arrays have been populated.
    QVERIFY(mPara->radArray[0] > 0);
    QVERIFY(mPara->numDensityArray[0] > 0);
    // Check that the sum of the number densities is not zero.
    double sumNumDensity = 0.0;
    for (unsigned int i = 0; i < mPara->nRadius; ++i) {
        sumNumDensity += mPara->numDensityArray[i];
    }
    QVERIFY(qFuzzyCompare(sumNumDensity, mPara->sphNumDensity));
}

// Test case for CalculateG.
void TestCalculate::test_CalculateG()
{
    std::complex<double> S1[361], S2[361];
    for (int i = 0; i < 361; ++i) {
        // Mock a simple phase function where S1 and S2 are constant.
        S1[i] = {1.0, 0.0};
        S2[i] = {1.0, 0.0};
    }
    // With constant S1/S2, g should be near 0
    double g_value = mCalc->CalculateG(S1, S2, mPara);

    //g should be close to 0.0
    QVERIFY(g_value< 1e-12);
}

// Test case for CalculateForwardBackward.
void TestCalculate::test_CalculateForwardBackward()
{
    std::complex<double> S1[361], S2[361];
    mPara->nTheta =361;
    for (unsigned int i = 0; i < mPara->nTheta; ++i) {
        // Mock a simple phase function where S1 and S2 are constant.
        S1[i] = {1.0, 0.0};
        S2[i] = {1.0, 0.0};
    }
    // The forward and backward scattering should be equal for a symmetric phase function.
    double forward = mCalc->CalculateForwardBackward(S1, S2, mPara, 0, ((mPara->nTheta-1)/2)+1);
    double backward = mCalc->CalculateForwardBackward(S1, S2, mPara, ((mPara->nTheta-1)/2), mPara->nTheta);

    // With constant S1/S2, (forward - backward) should be near 0
    QVERIFY(fabs(forward - backward) <1e-12);
}

// Test case for CalculatePowerLawAutoFitSimple.
void TestCalculate::test_CalculatePowerLawAutoFitSimple()
{
    //Test for complex power law with b=2.
    mPara->refWavel = 800;
    mPara->refWavelIdx = 3;
    mPara->muspAtRefWavel[3] = 10.0;
    mPara->nWavel = 9;

    mPara->wavelArray = new double[mPara->nWavel];
    mPara->mus = new double[mPara->nWavel];
    mPara->g = new double[mPara->nWavel];

    for (unsigned int i = 0; i < mPara->nWavel; ++i) {
        mPara->wavelArray[i] = 600.0 + 50.0 * i; // Wavelengths from 600 to 1000 nm.
        mPara->mus[i] = mPara->muspAtRefWavel[3] * pow(mPara->wavelArray[i] / mPara->refWavel, -2.0);
        mPara->g[i] = 0; //Assume g= 0;
    }
    mCalc->CalculatePowerLawAutoFitSimple(mPara);

    // The calculated bMie should be very close to 2.0.
    QCOMPARE(mPara->bMie, 2.0);
}

// Test case for CalculatePowerLawAutoFitComplex.
void TestCalculate::test_CalculatePowerLawAutoFitComplex()
{
    //Test for complex power law with fRay=0.5 and bMie=1.5.
    mPara->refWavel = 800;
    mPara->refWavelIdx = 3;
    mPara->muspAtRefWavel[3] = 10.0;
    mPara->nWavel = 9;

    mPara->wavelArray = new double[mPara->nWavel];
    mPara->mus = new double[mPara->nWavel];
    mPara->g = new double[mPara->nWavel];

    for (unsigned int i = 0; i < mPara->nWavel; ++i) {
        mPara->wavelArray[i] = 600.0 + 50.0 * i; // Wavelengths from 600 to 1000 nm.
        double x = mPara->wavelArray[i] / mPara->refWavel;
        mPara->mus[i] = mPara->muspAtRefWavel[3] * (0.5 * pow(x, -4.0) + (1.0 - 0.5) * pow(x, -1.5));
        mPara->g[i] = 0; //Assume g= 0;
    }
    mCalc->CalculatePowerLawAutoFitComplex(mPara);

    // The calculated values should be close to the expected values.
    QCOMPARE(mPara->fRay, 0.5);
    QCOMPARE(mPara->bMie, 1.5);
}

// Sanity check for the ComputeMuspAtRefWavel method.
void TestCalculate::test_ComputeMuspAtRefWavel_sanityCheck()
{
    mCalc->SetSphereRadiusAndRefIndex(mPara, 0, true);
    mCalc->ComputeMuspAtRefWavel(mPara);

    // The muspAtRefWavel array should contain values
    bool all_zero = true;
    for (int i = 0; i < 6; ++i) {
        if (mPara->muspAtRefWavel[i] != 0.0) {
            all_zero = false;
            break;
        }
    }
    QVERIFY(!all_zero);
}
