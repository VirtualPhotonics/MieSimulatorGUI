/**********************************************************************
** Tests for functions in calculate.cpp
**********************************************************************/

#include <QDebug>
#include <QtMath>
#include <QTest>
#include <complex>
#include "test_calculate.h"
#include "utilities.h"

TestCalculate::TestCalculate()
{
}

TestCalculate::~TestCalculate()
{
}

//Initialize values for tests
void TestCalculate::init()
{
    mPara = new Parameters();
    mCalc = new Calculate();

    mPara->radArray = new double[mPara->nRadius];
    mPara->numDensityArray = new double[mPara->nRadius];
    mPara->scatRefRealArray = new double[mPara->nRadius];
    mPara->scatRefImagArray = new double[mPara->nRadius];
    mPara->medRefArray = new double[mPara->nRadius];
    mPara->wavelArray = new double[mPara->nWavel];
    mPara->mus = new double[mPara->nWavel];
    mPara->g = new double[mPara->nWavel];
    mPara->cSca = new double[mPara->nWavel];
    mPara->cExt = new double[mPara->nWavel];
    mPara->cBack = new double[mPara->nWavel];
    mPara->sizePara = new double[mPara->nWavel];
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
    if (mPara->S1)
    {
        for (unsigned int w = 0; w < mPara->nWavel; ++w)
            delete[] mPara->S1[w];
        delete[] mPara->S1;
        mPara->S1 = nullptr;
    }
    if (mPara->S2)
    {
        for (unsigned int w = 0; w < mPara->nWavel; ++w)
            delete[] mPara->S2[w];
        delete[] mPara->S2;
        mPara->S2 = nullptr;
    }
    if (mPara->phaseFunctionAve)
    {
        for (unsigned int w = 0; w < mPara->nWavel; ++w)
            delete[] mPara->phaseFunctionAve[w];
        delete[] mPara->phaseFunctionAve;
        mPara->phaseFunctionAve = nullptr;
    }
    if (mPara->phaseFunctionPara)
    {
        for (unsigned int w = 0; w < mPara->nWavel; ++w)
            delete[] mPara->phaseFunctionPara[w];
        delete[] mPara->phaseFunctionPara;
        mPara->phaseFunctionPara = nullptr;
    }
    if (mPara->phaseFunctionPerp)
    {
        for (unsigned int w = 0; w < mPara->nWavel; ++w)
            delete[] mPara->phaseFunctionPerp[w];
        delete[] mPara->phaseFunctionPerp;
        mPara->phaseFunctionPerp = nullptr;
    }

    if (mPara->radArray) delete[] mPara->radArray;
    if (mPara->numDensityArray) delete[] mPara->numDensityArray;
    if (mPara->scatRefRealArray) delete[] mPara->scatRefRealArray;
    if (mPara->scatRefImagArray) delete[] mPara->scatRefImagArray;
    if (mPara->medRefArray) delete[] mPara->medRefArray;
    if (mPara->wavelArray) delete[] mPara->wavelArray;
    if (mPara->mus) delete[] mPara->mus;
    if (mPara->g) delete[] mPara->g;
    if (mPara->cSca) delete[] mPara->cSca;
    if (mPara->cExt) delete[] mPara->cExt;
    if (mPara->cBack) delete[] mPara->cBack;
    if (mPara->sizePara) delete[] mPara->sizePara;
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
    mPara->sizePara = nullptr;
    mPara->forward = nullptr;
    mPara->backward = nullptr;

    if (mCalc) delete mCalc;
    if (mPara) delete mPara;
    mCalc = nullptr;
    mPara = nullptr;
}

// Test case: Check DoSimulationdata
void TestCalculate::test_DoSimulationdata()
{
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

    double expectedMus = 93.36900541;
    double expectedG = 0.861445339;
    double expectedForward = 96.201215357;
    double expectedBackward = 3.7987846427;

    QVERIFY(fabs(mPara->mus[0] - expectedMus) <1e-8);
    QVERIFY(fabs(mPara->g[0] - expectedG) <1e-8);
    QVERIFY(fabs(mPara->forward[0] - expectedForward) <1e-8);
    QVERIFY(fabs(mPara->backward[0] - expectedBackward) <1e-8);
}

//Test case: Check real path vs complex path with tiny imaginary value
void TestCalculate::test_DoSimulation_RealVsComplexRefIndex()
{
    mPara->nRadius = 1;
    mPara->nWavel = 1;
    mPara->wavelArray[0] = 632.8;
    mPara->radArray[0] = 0.5;
    mPara->medRefArray[0] = 1.33;
    mPara->scatRefRealArray[0] = 1.5;
    mPara->numDensityArray[0] = 1e6;

    QLabel mockLabel;

    // Purely real
    mPara->scatRefImagArray[0] = 0.0;
    mCalc->DoSimulation(&mockLabel, mPara);
    double musReal = mPara->mus[0];
    QVERIFY(musReal > 0);

    // Complex (very small imaginary part)
    mPara->scatRefImagArray[0] = 1e-10;
    mCalc->DoSimulation(&mockLabel, mPara);

    QVERIFY(qAbs(mPara->mus[0] - musReal) < 1e-5);
}

//Test case: Check DoSimulation with multiple wavelengths
void TestCalculate::test_DoSimulation_MultiWavelength()
{
    mPara->nRadius = 1;
    mPara->nWavel = 2;
    mPara->nTheta = 181;
    mPara->wavelArray[0] = 500.0;
    mPara->wavelArray[1] = 800.0;
    mPara->radArray[0] = 0.5;
    mPara->medRefArray[0] = 1.33;
    mPara->scatRefRealArray[0] = 1.5;
    mPara->scatRefImagArray[0] = 0.0;
    mPara->numDensityArray[0] = 1e6;

    QLabel mockLabel;
    mCalc->DoSimulation(&mockLabel, mPara);

    // Mus should generally decrease as wavelength increases for this size
    QVERIFY(mPara->mus[0] > 0);
    QVERIFY(mPara->mus[1] > 0);
    QVERIFY(mPara->mus[0] != mPara->mus[1]);
}

// Test case: Test for zero density
void TestCalculate::test_DoSimulation_ZeroDensitySafety()
{
    mPara->nRadius = 1;
    mPara->nWavel = 2;
    mPara->nTheta = 181;
    mPara->wavelArray[0] = 500.0;
    mPara->wavelArray[1] = 800.0;
    mPara->radArray[0] = 0.5;
    mPara->medRefArray[0] = 1.33;
    mPara->scatRefRealArray[0] = 1.5;
    mPara->scatRefImagArray[0] = 0.0;
    mPara->numDensityArray[0] = 0;

    QLabel mockLabel;
    mCalc->DoSimulation(&mockLabel, mPara);

    QVERIFY(qIsNaN(mPara->g[0]) || mPara->g[0] == 0.0);
}

// Test case: Check DiameterRangeSetting with log-normal distribution
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

// Test case: Check DiameterRangeSetting with Gaussian distribution
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
    QCOMPARE(mPara->maxRadius - mPara->meanRadius, mPara->meanRadius - mPara->minRadius);
}

// Test case: Check DiameterRangeSetting for mono-disperse spheres
void TestCalculate::test_DiameterRangeSetting_monoDisperse()
{
    mPara->meanRadius = 5.0;
    mCalc->DiameterRangeSetting(mPara, 3);
    QCOMPARE(mPara->minRadius, 5.0);
    QCOMPARE(mPara->maxRadius, 5.0);
}

// Test case: SetSphereRadiusAndRefIndex for mono-disperse spheres
void TestCalculate::test_SetSphereRadiusAndRefIndex_monoDisperse()
{
    mPara->nRadius = 1;
    mPara->meanRadius = 1.0;
    mPara->sphNumDensity = 5e6;
    mPara->scatRefReal = 1.5;
    mPara->scatRefImag = 0.1;

    // The distribution index doesn't matter for nRadius=1.
    mCalc->SetSphereRadiusAndRefIndex(mPara, 0, true);

    //The elements in Array[0] are the same as the initial values
    QCOMPARE(mPara->radArray[0], mPara->meanRadius);
    QCOMPARE(mPara->numDensityArray[0], mPara->sphNumDensity);
    QCOMPARE(mPara->scatRefRealArray[0], mPara->scatRefReal);
    QCOMPARE(mPara->scatRefImagArray[0], mPara->scatRefImag);
}

// Test case: SetSphereRadiusAndRefIndex for poly-disperse with volume fraction
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
    for (unsigned int i = 0; i < mPara->nRadius; ++i)
    {
        sumNumDensity += mPara->numDensityArray[i];
    }
    QVERIFY(sumNumDensity > 0.0);
}

// Test case: SetSphereRadiusAndRefIndex for poly-disperse with concentration.
void TestCalculate::test_SetSphereRadiusAndRefIndex_poly_gaussian_conc()
{
    mPara->nRadius = 21;
    mPara->meanRadius = 0.5;
    mPara->stdDev = 0.2;
    mPara->sphNumDensity = 1e8;

    mCalc->DiameterRangeSetting(mPara, 1); // 1 = Gaussian
    mCalc->SetSphereRadiusAndRefIndex(mPara, 1, false); // false = use concentration.

    QCOMPARE(mPara->radArray[0], 1e-10);
    QVERIFY(mPara->numDensityArray[0] > 0);

    double sumNumDensity = 0.0;
    for (unsigned int i = 0; i < mPara->nRadius; ++i)
    {
        sumNumDensity += mPara->numDensityArray[i];
    }

    double relativeError = std::abs(sumNumDensity - mPara->sphNumDensity) / mPara->sphNumDensity;
    QVERIFY(relativeError < 1e-6);
}

// Test case: Check CalculateG.
void TestCalculate::test_CalculateG()
{
    std::complex<double> S1[361], S2[361];
    for (int i = 0; i < 361; ++i)
    {
        // Mock a simple phase function where S1 and S2 are constant.
        S1[i] = {1.0, 0.0};
        S2[i] = {1.0, 0.0};
    }
    // With constant S1/S2, g should be near 0
    double g_value = mCalc->CalculateG(S1, S2, mPara);

    //g should be close to 0.0
    QVERIFY(g_value< 1e-12);
}

// Test case: Check CalculateForwardBackward.
void TestCalculate::test_CalculateForwardBackward()
{
    std::complex<double> S1[361], S2[361];
    mPara->nTheta =361;
    for (unsigned int i = 0; i < mPara->nTheta; ++i)
    {
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

// Test case: Check CalculatePowerLawAutoFitSimple.
void TestCalculate::test_CalculatePowerLawAutoFitSimple()
{
    //Test for complex power law with b=2.
    double b = 2.0;
    mPara->refWavel = 800;
    mPara->refWavelIdx = mPara->wavel800;
    mPara->muspAtRefWavel[3] = 10.0;
    mPara->nWavel = 9;

    mPara->wavelArray = new double[mPara->nWavel];
    mPara->mus = new double[mPara->nWavel];
    mPara->g = new double[mPara->nWavel];

    for (unsigned int i = 0; i < mPara->nWavel; ++i)
    {
        mPara->wavelArray[i] = 600.0 + 50.0 * i; // Wavelengths from 600 to 1000 nm.
        mPara->mus[i] = mPara->muspAtRefWavel[3] * pow(mPara->wavelArray[i] / mPara->refWavel, -b);
        mPara->g[i] = 0; //Assume g= 0;
    }
    mCalc->CalculatePowerLawAutoFitSimple(mPara);

    // The calculated bMie should be very close to 2.0.
    QCOMPARE(mPara->bMie, b);
}

// Test case: Check CalculatePowerLawAutoFitComplex.
void TestCalculate::test_CalculatePowerLawAutoFitComplex()
{
    //Test for complex power law with fRay=0.5 and bMie=1.5.
    double fRay = 0.5;
    double bMie = 1.5;
    mPara->refWavel = 800;
    mPara->refWavelIdx = mPara->wavel800;
    mPara->muspAtRefWavel[3] = 10.0;
    mPara->nWavel = 9;

    mPara->wavelArray = new double[mPara->nWavel];
    mPara->mus = new double[mPara->nWavel];
    mPara->g = new double[mPara->nWavel];

    for (unsigned int i = 0; i < mPara->nWavel; ++i)
    {
        mPara->wavelArray[i] = 600.0 + 50.0 * i; // Wavelengths from 600 to 1000 nm.
        double x = mPara->wavelArray[i] / mPara->refWavel;
        mPara->mus[i] = mPara->muspAtRefWavel[3] * (fRay * pow(x, -4.0) + (1.0 - fRay) * pow(x, -bMie));
        mPara->g[i] = 0; //Assume g= 0;
    }
    mCalc->CalculatePowerLawAutoFitComplex(mPara);


    // The calculated values should be close to the expected values.
    QCOMPARE(mPara->fRay, fRay);
    QCOMPARE(mPara->bMie, bMie);
}

// Test case: Check Normalized Phase function
void TestCalculate::test_PhaseFunctionNormalization()
{
    mPara->nRadius = 1;
    mPara->nWavel = 1;
    mPara->nTheta = 361;
    mPara->wavelArray[0] = 500.0;
    mPara->radArray[0] = 0.5;
    mPara->medRefArray[0] = 1.33;
    mPara->scatRefRealArray[0] = 1.5;
    mPara->scatRefImagArray[0] = 0.0;
    mPara->numDensityArray[0] = 1e6;

    QLabel mockLabel;
    mCalc->DoSimulation(&mockLabel, mPara);

    double integral = 0.0;
    double stepTheta = M_PI / (mPara->nTheta - 1);
    Utilities util;

    for (unsigned int t = 0; t < mPara->nTheta; ++t)
    {
        double theta = t * stepTheta;
        double weight = util.SimpsonsWeight(t, mPara->nTheta);
        integral += mPara->phaseFunctionAve[0][t] * sin(theta) * 2.0 * M_PI * weight;
    }
    double unity = integral * stepTheta;

    // Ideally, the normalized phase function integral should be 1.0.
    QVERIFY(qAbs(unity - 1.0) < 1e-6);
}

// Sanity check for the ComputeMuspAtRefWavel method.
void TestCalculate::test_ComputeMuspAtRefWavel_sanityCheck()
{
    mCalc->SetSphereRadiusAndRefIndex(mPara, 0, true);
    mCalc->ComputeMuspAtRefWavel(mPara);

    // The muspAtRefWavel array should contain values
    bool all_zero = true;
    for (int i = 0; i < 6; ++i)
    {
        if (mPara->muspAtRefWavel[i] != 0.0)
        {
            all_zero = false;
            break;
        }
    }
    QVERIFY(!all_zero);
}

// Test case: High concentration regime (f_v > 0.1)
void TestCalculate::test_CheckIndependentScattering_HighConcentration()
{
    double clearanceToWavelength, sizeParameter, volFraction, criticalWavelength;
    QString strRegime;
    bool flagVolFrac = false;

    mPara->nRadius = 1;
    mPara->meanRadius = 1.0;
    mPara->startWavel = 500.0;
    mPara->numDensityArray[0] = 3e7;
    mPara->medRef = 1.0;

    bool isDependent = mCalc->CheckIndependentScattering(mPara, clearanceToWavelength, sizeParameter,
                                                         volFraction, criticalWavelength, strRegime,
                                                         flagVolFrac);
    QVERIFY(isDependent == true);
}

// Test case: Low concentration, size parameter below threshold (f_v <= 0.006, sizePara <= 0.388)
void TestCalculate::test_CheckIndependentScattering_LowConcSmallSize()
{
    double clearanceToWavelength, sizeParameter, volFraction, criticalWavelength;
    QString strRegime;
    bool flagVolFrac = false;

    mPara->nRadius = 1;
    mPara->meanRadius = 0.01;
    mPara->startWavel = 1000.0;
    mPara->numDensityArray[0] = 1e5;
    mPara->medRef = 1.0;

    bool isDependent = mCalc->CheckIndependentScattering(mPara, clearanceToWavelength, sizeParameter,
                                                         volFraction, criticalWavelength, strRegime,
                                                         flagVolFrac);
    QVERIFY(isDependent == false);
}

// Test case: Transitional regime, Dependent (0.006 < f_v <= 0.1 and low clearance)
void TestCalculate::test_CheckIndependentScattering_TransitionalDependent()
{
    double clearanceToWavelength, sizeParameter, volFraction, criticalWavelength;
    QString strRegime;
    bool flagVolFrac = false;

    mPara->nRadius = 1;
    mPara->meanRadius = 1.0;
    mPara->startWavel = 1000.0;
    mPara->numDensityArray[0] = 5e6;
    mPara->medRef = 1.0;

    // interParticleDistance = 1000 / (5e6)^(1/3) = 5.84
    // clearanceToWavelength = (5.84 - 2*1.0) / 1.0 = 3.84
    // 3.84 <= 5.0 (required clearance)
    bool isDependent = mCalc->CheckIndependentScattering(mPara, clearanceToWavelength, sizeParameter,
                                                         volFraction, criticalWavelength, strRegime,
                                                         flagVolFrac);
    QVERIFY(isDependent == true);
}

// Test case: Low concentration, size parameter above threshold, Independent (High clearance)
void TestCalculate::test_CheckIndependentScattering_LowConcLargeSizeIndependent()
{
    double clearanceToWavelength, sizeParameter, volFraction, criticalWavelength;
    QString strRegime;
    bool flagVolFrac = false;

    mPara->nRadius = 1;
    mPara->meanRadius = 1.0;
    mPara->startWavel = 500.0;
    mPara->numDensityArray[0] = 1e4;
    mPara->medRef = 1.0;

    // interParticleDistance = 1000 / (1e4)^(1/3) = 46.4
    // clearanceToWavelength = (46.4 - 2) / 0.5 = 88.8
    // 88.8 > 5.0, should be independent (false)
    bool isDependent = mCalc->CheckIndependentScattering(mPara, clearanceToWavelength, sizeParameter,
                                                         volFraction, criticalWavelength, strRegime,
                                                         flagVolFrac);
    QVERIFY(isDependent == false);
}

// Test case: Polydisperse branch logic check
void TestCalculate::test_CheckIndependentScattering_Polydisperse()
{
    double clearanceToWavelength, sizeParameter, volFraction, criticalWavelength;
    QString strRegime;
    bool flagVolFrac = false;

    mPara->nRadius = 2;
    mPara->radArray[0] = 0.5;
    mPara->radArray[1] = 1.5;
    mPara->numDensityArray[0] = 1e6;
    mPara->numDensityArray[1] = 1e6;
    mPara->startWavel = 1000.0;
    mPara->medRef = 1.0;

    bool isDependent = mCalc->CheckIndependentScattering(mPara, clearanceToWavelength, sizeParameter,
                                                         volFraction, criticalWavelength, strRegime,
                                                         flagVolFrac);
    QVERIFY(isDependent == false);
}
