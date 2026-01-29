#ifndef TEST_PARAMETERS_H
#define TEST_PARAMETERS_H

#include "../app/parameters.h"

// Define a test class
class TestParameters : public QObject
{
    Q_OBJECT

public:
    TestParameters();
    ~TestParameters();

private slots:
    void init();
    void cleanup();

    void test_CheckValidityCommonParameters_valid();
    void test_CheckValidityCommonParameters_invalidRefractiveIndex();
    void test_CheckValidityCommonParameters_invalidRelativeRefractiveIndexOne();
    void test_CheckValidityCommonParameters_invalidRelativeRefractiveIndexRange();
    void test_CheckValidityCommonParameters_invalidImaginaryRefractiveIndex();
    void test_CheckValidityCommonParameters_invalidWavelengthZero();
    void test_CheckValidityCommonParameters_invalidWavelengthStepZero();
    void test_CheckValidityCommonParameters_invalidStartWavelengthRange();
    void test_CheckValidityCommonParameters_invalidEndWavelengthRange();
    void test_CheckValidityCommonParameters_invalidWavelengthOrder();
    void test_CheckValidityCommonParameters_invalidNumDensityZero();
    void test_CheckValidityCommonParameters_invalidVolFractionZero();
    void test_CheckValidityCommonParameters_invalidVolFractionUpperLimit();
    void test_CheckValidityCommonParameters_invalidMeanRadiusRange();
    void test_CheckValidityCommonParameters_invalidVolumeExceedsLimit();
    void test_CheckValidityCommonParameters_invalidRelativeRefractiveIndexSmall();
    void test_CheckValidityCommonParameters_volFractionAtLimit();
    void test_CheckValidityCommonParameters_invalidRelativeRefractiveIndexTooLow();

    void test_CheckValidityDistributionParameters_valid();
    void test_CheckValidityDistributionParameters_invalidStdDevZero();
    void test_CheckValidityDistributionParameters_invalidStdDevLogNormalUpper();
    void test_CheckValidityDistributionParameters_invalidStdDevLogNormalLower();
    void test_CheckValidityDistributionParameters_invalidStdDevGaussianUpper();
    void test_CheckValidityDistributionParameters_invalidStdDevGaussianLower();
    void test_CheckValidityDistributionParameters_invalidNRadiusOne();
    void test_CheckValidityDistributionParameters_invalidNRadiusRange();
    void test_CheckValidityDistributionParameters_invalidMeanRadiusDistributionRange();
    void test_CheckValidityDistributionParameters_invalidStdDevMeanRadiusRatio();
    void test_CheckValidityDistributionParameters_invalidMeanRadiusLower();

    void test_CheckPackingVolume_invalid();
    void test_CheckPackingVolume_valid();

private:
    Parameters *mPara;
};

#endif // TEST_PARAMETERS_H
