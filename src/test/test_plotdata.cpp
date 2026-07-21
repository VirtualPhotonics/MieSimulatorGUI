/**********************************************************************
** Tests for functions in dialog/plotdata.cpp
**********************************************************************/

#include <QtMath>
#include <QTest>
#include "test_plotdata.h"

//Initialize values for tests
void TestPlotData::init()
{
    mPara = new Parameters();
    mPlotData = new PlotData();
}

// clean up resources
void TestPlotData::cleanup()
{
    if (mPara->phaseFunctionPara)
    {
        for (unsigned int i = 0; i < mPara->nWavel; i++)
        {
            delete [] mPara->phaseFunctionPara[i];
            delete [] mPara->phaseFunctionPerp[i];
            delete [] mPara->phaseFunctionAve[i];
        }
        delete [] mPara->phaseFunctionPara;
        delete [] mPara->phaseFunctionPerp;
        delete [] mPara->phaseFunctionAve;
        mPara->phaseFunctionPara = nullptr;
        mPara->phaseFunctionPerp = nullptr;
        mPara->phaseFunctionAve = nullptr;
    }

    delete mPara;
    mPara = nullptr;
    delete mPlotData;
    mPlotData = nullptr;
}

// Test case: FindMinLogPolarPlot with a single wavelength.
void TestPlotData::test_FindMinLogPolarPlot_singleWavelength()
{
    mPara->nWavel = 1;
    mPara->nTheta = 3;

    mPara->phaseFunctionPara = new double*[1];
    mPara->phaseFunctionPerp = new double*[1];
    mPara->phaseFunctionAve = new double*[1];
    mPara->phaseFunctionPara[0] = new double[3]{1.0, 0.1, 10.0};
    mPara->phaseFunctionPerp[0] = new double[3]{2.0, 0.01, 5.0};
    mPara->phaseFunctionAve[0] = new double[3]{1.5, 0.001, 7.5};

    // Smallest value across all three series is 0.001 -> log10(0.001) = -3
    double result = mPlotData->FindMinLogPolarPlot(mPara);
    QCOMPARE(result, floor(log10(0.001)));
}

// Test case: FindMinLogPolarPlot where the minimum occurs in a wavelength
// other than the first one, to exercise the running-minimum comparison
// across the outer wavelength loop.
void TestPlotData::test_FindMinLogPolarPlot_minimumInLaterWavelength()
{
    mPara->nWavel = 2;
    mPara->nTheta = 2;

    mPara->phaseFunctionPara = new double*[2];
    mPara->phaseFunctionPerp = new double*[2];
    mPara->phaseFunctionAve = new double*[2];

    // Wavelength 0: larger minimum
    mPara->phaseFunctionPara[0] = new double[2]{1.0, 2.0};
    mPara->phaseFunctionPerp[0] = new double[2]{1.0, 2.0};
    mPara->phaseFunctionAve[0] = new double[2]{1.0, 2.0};

    // Wavelength 1: contains the overall smallest value
    mPara->phaseFunctionPara[1] = new double[2]{1e-6, 1.0};
    mPara->phaseFunctionPerp[1] = new double[2]{1.0, 1.0};
    mPara->phaseFunctionAve[1] = new double[2]{1.0, 1.0};

    double result = mPlotData->FindMinLogPolarPlot(mPara);
    QCOMPARE(result, floor(log10(1e-6)));
}

// Test case: RearrangePhaseFunctionData in linear mode with theta expressed
// as signed degrees (-180 to +180).
void TestPlotData::test_RearrangePhaseFunctionData_negPosDegrees_linear()
{
    mPara->nWavel = 1;
    mPara->nTheta = 3;
    int totalSize = 2 * static_cast<int>(mPara->nTheta) - 1; // 5

    mPara->phaseFunctionPara = new double*[1];
    mPara->phaseFunctionPerp = new double*[1];
    mPara->phaseFunctionAve = new double*[1];
    mPara->phaseFunctionPara[0] = new double[3]{1.0, 2.0, 3.0};
    mPara->phaseFunctionPerp[0] = new double[3]{4.0, 5.0, 6.0};
    mPara->phaseFunctionAve[0] = new double[3]{7.0, 8.0, 9.0};

    QVector<double> theta(totalSize);
    QVector<double> phaseFuncPara(totalSize);
    QVector<double> phaseFuncPerp(totalSize);
    QVector<double> phaseFuncAve(totalSize);

    bool flagLinearLog = true;         // true: no log transform
    bool flagThetaNegPosOrPos = true;  // -180 to +180 degrees

    mPlotData->RearrangePhaseFunctionData(mPara, theta, phaseFuncPara, phaseFuncPerp, phaseFuncAve,
                                          0, flagLinearLog, flagThetaNegPosOrPos);

    // i=0 -> angleDeg=0, forwardIndex=backwardIndex=2
    QCOMPARE(theta[2], 0.0);
    // i=1 -> angleDeg=90, forwardIndex=3, backwardIndex=1
    QCOMPARE(theta[3], 90.0);
    QCOMPARE(theta[1], -90.0);
    // i=2 -> angleDeg=180, forwardIndex=4, backwardIndex=0
    QCOMPARE(theta[4], 180.0);
    QCOMPARE(theta[0], -180.0);

    // Linear mode: values copied as-is, mirrored to forward/backward indices
    QCOMPARE(phaseFuncPara[3], 2.0);
    QCOMPARE(phaseFuncPara[1], 2.0);
    QCOMPARE(phaseFuncPerp[4], 6.0);
    QCOMPARE(phaseFuncAve[2], 7.0);
}

// Test case: RearrangePhaseFunctionData in log mode with theta expressed
// in radians over 0 to 2*pi.
void TestPlotData::test_RearrangePhaseFunctionData_zeroToTwoPiRadians_log()
{
    mPara->nWavel = 1;
    mPara->nTheta = 2;
    int totalSize = 2 * static_cast<int>(mPara->nTheta) - 1; // 3

    mPara->phaseFunctionPara = new double*[1];
    mPara->phaseFunctionPerp = new double*[1];
    mPara->phaseFunctionAve = new double*[1];
    mPara->phaseFunctionPara[0] = new double[2]{1.0, 100.0};
    mPara->phaseFunctionPerp[0] = new double[2]{1.0, 100.0};
    mPara->phaseFunctionAve[0] = new double[2]{1.0, 100.0};

    QVector<double> theta(totalSize);
    QVector<double> phaseFuncPara(totalSize);
    QVector<double> phaseFuncPerp(totalSize);
    QVector<double> phaseFuncAve(totalSize);

    bool flagLinearLog = false;         // false: apply log10 transform
    bool flagThetaNegPosOrPos = false;  // 0 to 2*pi radians

    mPlotData->RearrangePhaseFunctionData(mPara, theta, phaseFuncPara, phaseFuncPerp, phaseFuncAve,
                                          0, flagLinearLog, flagThetaNegPosOrPos);

    // i=0 -> angleDeg=0, forwardIndex=backwardIndex=1, special-cased to 0 (not 2*pi)
    QCOMPARE(theta[1], 0.0);

    // i=1 -> angleDeg=180, forwardIndex=2, backwardIndex=0.
    // Both radPos and radNeg reduce to pi for a 180 degree angle.
    QVERIFY(qFuzzyCompare(theta[2], M_PI));
    QVERIFY(qFuzzyCompare(theta[0], M_PI));

    // Log mode: value at theta index 1 comes from phaseFunctionPara[0][0] = 1.0 -> log10(1.0+tiny) ~ 0
    QCOMPARE(phaseFuncPara[1], 0.0);
    // value at theta index 2 (and mirrored index 0) comes from phaseFunctionPara[0][1] = 100.0 -> log10(100.0) = 2
    QVERIFY(qFuzzyCompare(phaseFuncPara[2], 2.0));
    QVERIFY(qFuzzyCompare(phaseFuncPara[0], 2.0));
}