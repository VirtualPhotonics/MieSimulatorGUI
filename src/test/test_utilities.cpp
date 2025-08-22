/**********************************************************************
** Tests for functions in utilities.cpp
**********************************************************************/

#include <QtMath>
#include <QTime>
#include <QTest>
#include <QCoreApplication>
#include "test_utilities.h"

TestUtilities::TestUtilities()
{
}

TestUtilities::~TestUtilities()
{
}

//Initialize values for tests
void TestUtilities::init()
{
    mUtil = new Utilities();
}

// clean up resources
void TestUtilities::cleanup()
{
    if (mUtil) {
        delete mUtil;
        mUtil = nullptr;
    }
}

// Test case: ComplexAbsSquared function.
void TestUtilities::test_ComplexAbsSquared()
{
    // Test with a simple real number
    std::complex<double> realComplex(3.0, 0.0);
    QCOMPARE(mUtil->ComplexAbsSquared(realComplex), 9.0);

    // Test with a simple imaginary number
    std::complex<double> imagComplex(0.0, 4.0);
    QCOMPARE(mUtil->ComplexAbsSquared(imagComplex), 16.0);

    // Test with a combination of real and imaginary parts
    std::complex<double> complexComplex(3.0, 4.0);
    QCOMPARE(mUtil->ComplexAbsSquared(complexComplex), 25.0);

    // Test with a negative number
    std::complex<double> negComplex(-2.0, -1.0);
    QCOMPARE(mUtil->ComplexAbsSquared(negComplex), 5.0);

    // Test with zero
    std::complex<double> zeroComplex(0.0, 0.0);
    QCOMPARE(mUtil->ComplexAbsSquared(zeroComplex), 0.0);
}

// Test case: ComplexAbs function.
void TestUtilities::test_ComplexAbs()
{
    // Test with a simple real number
    std::complex<double> realComplex(5.0, 0.0);
    QVERIFY(qFuzzyCompare(mUtil->ComplexAbs(realComplex), 5.0));

    // Test with a simple imaginary number
    std::complex<double> imagComplex(0.0, -12.0);
    QVERIFY(qFuzzyCompare(mUtil->ComplexAbs(imagComplex), 12.0));

    // Test with a combination of real and imaginary parts
    std::complex<double> complexComplex(3.0, 4.0);
    QVERIFY(qFuzzyCompare(mUtil->ComplexAbs(complexComplex), 5.0));

    // Test with zero
    std::complex<double> zeroComplex(0.0, 0.0);
    QVERIFY(qFuzzyCompare(mUtil->ComplexAbs(zeroComplex), 0.0));
}

// Test case: SimpsonsWeight function.
void TestUtilities::test_SimpsonsWeight()
{
    unsigned int n = 5; // A simple odd number of points (n=5, i=0 to 4)
    // First and last points
    QVERIFY(qFuzzyCompare(mUtil->SimpsonsWeight(0, n), 1.0/3.0));
    QVERIFY(qFuzzyCompare(mUtil->SimpsonsWeight(4, n), 1.0/3.0));

    // Odd indexed points (1, 3)
    QVERIFY(qFuzzyCompare(mUtil->SimpsonsWeight(1, n), 4.0/3.0));
    QVERIFY(qFuzzyCompare(mUtil->SimpsonsWeight(3, n), 4.0/3.0));

    // Even indexed points (2)
    QVERIFY(qFuzzyCompare(mUtil->SimpsonsWeight(2, n), 2.0/3.0));
}

//Test case: NiceStep for ticks
void TestUtilities::test_NiceStep()
{
    //Test for a small range
    double range1 = 1.0;
    int initialCircles1 = 10;
    QCOMPARE(mUtil->NiceStep(range1, initialCircles1), 0.1);

    // Test for fraction < 3.0
    double range2 = 25.0;
    int initialCircles2 = 10;
    QCOMPARE(mUtil->NiceStep(range2, initialCircles2), 2.0);

    // Test for fraction < 7.5
    double range3 = 55.0;
    int initialCircles3 = 10;
    QCOMPARE(mUtil->NiceStep(range3, initialCircles3), 5.0);

    // Test for fraction >= 7.5
    double range4 = 80.0;
    int initialCircles4 = 10;
    QCOMPARE(mUtil->NiceStep(range4, initialCircles4), 10.0);

    // Range is a power of 10
    double range5 = 100.0;
    int initialCircles5 = 10;
    QCOMPARE(mUtil->NiceStep(range5, initialCircles5), 10.0);

    // A zero range
    double range6 = 0.0;
    int initialCircles6 = 1;
    QCOMPARE(mUtil->NiceStep(range6, initialCircles6), 0.0);

    // A zero initialCircle
    double range7 = 1.0;
    int initialCircles7 = 0;
    QCOMPARE(mUtil->NiceStep(range7, initialCircles7), 1.0);
}

//Test case: Find minimum  or maximum within three vectors
void TestUtilities::test_FindMinMax()
{
    // Finding the minimum value
    QVector<double> yPara = {10.5, 20.0, 5.2};
    QVector<double> yPerp = {15.1, 7.8, 22.3};
    QVector<double> yAve = {3.4, 18.9, 12.0};
    bool findMin = true;
    double result = mUtil->FindMinMax(yPara, yPerp, yAve, findMin);
    QCOMPARE(result, 3.4);

    // Finding the maximum value
    yPara = {-10.5, -20.0, 5.2};
    yPerp = {15.1, 7.8, -22.3};
    yAve = {3.4, 18.9, 12.0};
    findMin = true;
    result = mUtil->FindMinMax(yPara, yPerp, yAve, findMin);
    QCOMPARE(result, -22.3);

    // Empty vector
    yPara = {10.5, 20.0, 5.2};
    QVector<double> yPerp0;
    yAve = {3.4, 18.9, 12.0};
    findMin = false;
    result = mUtil->FindMinMax(yPara, yPerp0, yAve, findMin);
    QCOMPARE(result, 20.0);

    // All vectors being empty vectors
    QVector<double> yPara0;
    QVector<double> yAve0;
    findMin = true;
    result = mUtil->FindMinMax(yPara0, yPerp0, yAve0, findMin);
    QVERIFY(std::isnan(result));
}
