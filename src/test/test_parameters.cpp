/**********************************************************************
** Tests for functions in parameters.cpp
**********************************************************************/

#include "test_parameters.h"
#include <QTest>

TestParameters::TestParameters()
{
}

TestParameters::~TestParameters()
{
}

//Initialize values for tests
void TestParameters::init()
{
    mPara = new Parameters();
    mPara->scatRefReal = 1.377;
    mPara->medRef = 1.333;
    mPara->scatRefImag = 0.0;
    mPara->startWavel = 600;
    mPara->endWavel = 1000;
    mPara->stepWavel = 10;
    mPara->sphNumDensity = 1e8;
    mPara->volFraction = 0.1;
    mPara->meanRadius = 0.1;
    mPara->stdDev = 0.25;
    mPara->nRadius = 31;
}

// clean up resources
void TestParameters::cleanup()
{
    delete mPara;
    mPara = nullptr;
}

// Test case: Check CheckValidityCommonParameters
void TestParameters::test_CheckValidityCommonParameters_valid()
{
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, true, false);
    QVERIFY(result.isValid);
}

// Test case: Refractive index is zero
void TestParameters::test_CheckValidityCommonParameters_invalidRefractiveIndex()
{

    mPara->scatRefReal = 0.0;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, false, false);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Refractive index cannot be zero."));
}

// Test case: Relative refractive index is 1.0
void TestParameters::test_CheckValidityCommonParameters_invalidRelativeRefractiveIndexOne()
{    
    mPara->scatRefReal = 1.5;
    mPara->medRef = 1.5;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, false, false);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Relative refractive index cannot be 1.0."));
}

// Test case: Relative refractive index is out of range (m > 5.0)
void TestParameters::test_CheckValidityCommonParameters_invalidRelativeRefractiveIndexRange()
{    
    mPara->scatRefReal = 6.0;
    mPara->medRef = 1.0;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, false, false);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Unrealistic relative refractive index! Check sphere and medium refractive index values."));
}

// Test case: Imaginary refractive index is too large
void TestParameters::test_CheckValidityCommonParameters_invalidImaginaryRefractiveIndex()
{    
    mPara->scatRefImag = 5.0;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, false, false);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Imaginary refractive index must be negative or less than 5.0."));
}

// Test case: Starting wavelength is zero
void TestParameters::test_CheckValidityCommonParameters_invalidWavelengthZero()
{    
    mPara->startWavel = 0.0;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, false, false);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("The starting or ending wavelength cannot be zero."));
}

 // Test case: Wavelength step is zero
void TestParameters::test_CheckValidityCommonParameters_invalidWavelengthStepZero()
{   
    mPara->stepWavel = 0.0;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, false, false);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Wavelength step cannot be zero."));
}

// Test case: Starting wavelength is too small
void TestParameters::test_CheckValidityCommonParameters_invalidStartWavelengthRange()
{    
    mPara->startWavel = 49.0;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, false, false);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Current minimum wavlength is 50nm."));
}

// Test case: Ending wavelength is too large
void TestParameters::test_CheckValidityCommonParameters_invalidEndWavelengthRange()
{    
    mPara->endWavel = 3001.0;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, false, false);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Current maximum wavlength is 3000nm."));
}

// Test case: Starting wavelength is greater than ending wavelength
void TestParameters::test_CheckValidityCommonParameters_invalidWavelengthOrder()
{    
    mPara->startWavel = 1100.0;
    mPara->endWavel = 1000.0;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, false, false);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("The starting wavelength is greater than the ending wavelength."));
}

// Test case: Number density is zero when selected
void TestParameters::test_CheckValidityCommonParameters_invalidNumDensityZero()
{    
    mPara->sphNumDensity = 0.0;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, true, false);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Sphere concentration cannot be zero."));
}

// Test case: Volume fraction is zero when selected
void TestParameters::test_CheckValidityCommonParameters_invalidVolFractionZero()
{    
    mPara->volFraction = 0.0;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, false, true);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Volume Fraction cannot be zero."));
}

// Test case: Volume fraction exceeds 1.0 when selected
void TestParameters::test_CheckValidityCommonParameters_invalidVolFractionUpperLimit()
{    
    mPara->volFraction = 1.1;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, false, true);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Volume Fraction must not exceed the maximum packing factor (~0.74)."));
}

// Test case: Mean radius is too large
void TestParameters::test_CheckValidityCommonParameters_invalidMeanRadiusRange()
{    
    mPara->meanRadius = 151.0;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, false, false);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Diameter is out of range! Enter a value between 0.0001μm and 300μm."));
}

// Test case: Concentration * volume exceeds the limit for mono-disperse
void TestParameters::test_CheckValidityCommonParameters_invalidVolumeExceedsLimit()
{
    mPara->meanRadius = 100.0; // Radius in micrometers
    mPara->sphNumDensity = 239.0;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(true, true, false);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("'Concentration x Sphere Volume' exceeds the maximum packing factor! Reduce Concentration (Conc)."));
}

// Test case: Relative refractive index is too small (m < 0.05)
void TestParameters::test_CheckValidityCommonParameters_invalidRelativeRefractiveIndexSmall()
{
    mPara->scatRefReal = 0.04;
    mPara->medRef = 1.0;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, false, false);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Unrealistic relative refractive index! Check sphere and medium refractive index values."));
}

// Test case: Volume fraction exactly at the maximum packing factor limit
void TestParameters::test_CheckValidityCommonParameters_volFractionAtLimit()
{
    mPara->volFraction = 0.74048;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, false, true);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Volume Fraction must not exceed the maximum packing factor (~0.74)."));
}

// Test case: Mean radius is too small (Distribution Parameters)
void TestParameters::test_CheckValidityCommonParameters_invalidRelativeRefractiveIndexTooLow()
{
    mPara->meanRadius = 0.0004; // Below 0.0005
    ParameterValidationResult result = mPara->CheckValidityDistributionParameters(0);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Diameter is out of range! Enter a value between 0.001μm and 50μm."));
}

// Test case: Ensures a valid set of parameters passes the distribution check (for Log Normal)
void TestParameters::test_CheckValidityDistributionParameters_valid()
{    
    ParameterValidationResult result = mPara->CheckValidityDistributionParameters(0);
    QVERIFY(result.isValid);
}

// Test case:Standard deviation is zero
void TestParameters::test_CheckValidityDistributionParameters_invalidStdDevZero()
{    
    mPara->stdDev = 0.0;
    ParameterValidationResult result = mPara->CheckValidityDistributionParameters(0);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Standard Deviation is zero! Use 'Mono Disperse'."));
}

// Test case: Log Normal distribution with a large standard deviation
void TestParameters::test_CheckValidityDistributionParameters_invalidStdDevLogNormalUpper()
{    
    mPara->stdDev = 3.1;
    ParameterValidationResult result = mPara->CheckValidityDistributionParameters(0);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Large standard deviation provides an abnormal Log Normal distribution! Current limit for Log Normal is 3.0μm."));
}

// Test case: Log Normal distribution with a small standard deviation
void TestParameters::test_CheckValidityDistributionParameters_invalidStdDevLogNormalLower()
{    
    mPara->stdDev = 1e-6;
    ParameterValidationResult result = mPara->CheckValidityDistributionParameters(0);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("The standard deviation is too small! Use 'Mono Disperse'."));
}

// Test case: Gaussian distribution with a large standard deviation
void TestParameters::test_CheckValidityDistributionParameters_invalidStdDevGaussianUpper()
{    
    mPara->stdDev = 51.0;
    ParameterValidationResult result = mPara->CheckValidityDistributionParameters(1);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("The standard deviation is too large! Current limit for Gaussian is 50.0μm."));
}

// Test case: Gaussian distribution with a small standard deviation
void TestParameters::test_CheckValidityDistributionParameters_invalidStdDevGaussianLower()
{    
    mPara->stdDev = 1e-9;
    ParameterValidationResult result = mPara->CheckValidityDistributionParameters(1);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("The standard deviation is too small! Use 'Mono Disperse'."));
}

// Test case: nRadius is 1
void TestParameters::test_CheckValidityDistributionParameters_invalidNRadiusOne()
{    
    mPara->nRadius = 1;
    ParameterValidationResult result = mPara->CheckValidityDistributionParameters(0);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Discrete sphere size is 1! Use 'Mono Disperse'."));
}

// Test case: nRadius is out of range
void TestParameters::test_CheckValidityDistributionParameters_invalidNRadiusRange()
{
    mPara->nRadius = 202;
    ParameterValidationResult result = mPara->CheckValidityDistributionParameters(0);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Number of sphere sizes is out of range! Enter a value between 2 and 201."));
}

// Test case: mean radius is too large for distribution check
void TestParameters::test_CheckValidityDistributionParameters_invalidMeanRadiusDistributionRange()
{    
    mPara->meanRadius = 26.0;
    ParameterValidationResult result = mPara->CheckValidityDistributionParameters(0);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Diameter is out of range! Enter a value between 0.001μm and 50μm."));
}

// Test case: Ratio of stdDev/meanRadius is too small
void TestParameters::test_CheckValidityDistributionParameters_invalidStdDevMeanRadiusRatio()
{    
    mPara->stdDev = 1e-4;
    mPara->meanRadius = 10.0;
    ParameterValidationResult result = mPara->CheckValidityDistributionParameters(0);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Standard deviation to mean diameter ratio is smaller than 1e-5! Use 'Mono Disperse'."));
}

// Test case: Relative refractive index is too low (m < 0.05)
void TestParameters::test_CheckValidityDistributionParameters_invalidMeanRadiusLower()
{
    mPara->scatRefReal = 0.04;
    mPara->medRef = 1.0;
    ParameterValidationResult result = mPara->CheckValidityCommonParameters(false, false, false);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Unrealistic relative refractive index! Check sphere and medium refractive index values."));
}

// Test case: Polydisperse packing volume is valid
void TestParameters::test_CheckPackingVolume_invalid()
{
    double totalVolume = 0.5;
    ParameterValidationResult result = mPara->CheckValidityPackingVolume(totalVolume);
    QVERIFY(result.isValid);
}

// Test case: Polydisperse packing volume is valid
void TestParameters::test_CheckPackingVolume_valid()
{
    double totalVolume = 0.8;
    ParameterValidationResult result = mPara->CheckValidityPackingVolume(totalVolume);
    QVERIFY(!result.isValid);
    QCOMPARE(result.errorMessage, QString("Total sphere volume in 1mm³ exceeds the maximum packing factor. Reduce Concentration."));
}
