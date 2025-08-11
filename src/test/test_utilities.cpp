/**********************************************************************
** Tests for functions in utilities.cpp
**********************************************************************/

#include <QtMath>
#include <QTime>
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
    mUtilities = new utilities();
}

// clean up resources
void TestUtilities::cleanup()
{
    if (mUtilities) {
        delete mUtilities;
        mUtilities = nullptr;
    }
}

// Test case: ComplexAbsSquared function.
void TestUtilities::test_ComplexAbsSquared()
{
    // Test with a simple real number
    std::complex<double> realComplex(3.0, 0.0);
    QCOMPARE(mUtilities->ComplexAbsSquared(realComplex), 9.0);

    // Test with a simple imaginary number
    std::complex<double> imagComplex(0.0, 4.0);
    QCOMPARE(mUtilities->ComplexAbsSquared(imagComplex), 16.0);

    // Test with a combination of real and imaginary parts
    std::complex<double> complexComplex(3.0, 4.0);
    QCOMPARE(mUtilities->ComplexAbsSquared(complexComplex), 25.0);

    // Test with a negative number
    std::complex<double> negComplex(-2.0, -1.0);
    QCOMPARE(mUtilities->ComplexAbsSquared(negComplex), 5.0);

    // Test with zero
    std::complex<double> zeroComplex(0.0, 0.0);
    QCOMPARE(mUtilities->ComplexAbsSquared(zeroComplex), 0.0);
}

// Test case: ComplexAbs function.
void TestUtilities::test_ComplexAbs()
{
    // Test with a simple real number
    std::complex<double> realComplex(5.0, 0.0);
    QVERIFY(qFuzzyCompare(mUtilities->ComplexAbs(realComplex), 5.0));

    // Test with a simple imaginary number
    std::complex<double> imagComplex(0.0, -12.0);
    QVERIFY(qFuzzyCompare(mUtilities->ComplexAbs(imagComplex), 12.0));

    // Test with a combination of real and imaginary parts
    std::complex<double> complexComplex(3.0, 4.0);
    QVERIFY(qFuzzyCompare(mUtilities->ComplexAbs(complexComplex), 5.0));

    // Test with zero
    std::complex<double> zeroComplex(0.0, 0.0);
    QVERIFY(qFuzzyCompare(mUtilities->ComplexAbs(zeroComplex), 0.0));
}

// Test case: SimpsonsWeight function.
void TestUtilities::test_SimpsonsWeight()
{
    unsigned int n = 5; // A simple odd number of points (n=5, i=0 to 4)
    // First and last points
    QVERIFY(qFuzzyCompare(mUtilities->SimpsonsWeight(0, n), 1.0/3.0));
    QVERIFY(qFuzzyCompare(mUtilities->SimpsonsWeight(4, n), 1.0/3.0));

    // Odd indexed points (1, 3)
    QVERIFY(qFuzzyCompare(mUtilities->SimpsonsWeight(1, n), 4.0/3.0));
    QVERIFY(qFuzzyCompare(mUtilities->SimpsonsWeight(3, n), 4.0/3.0));

    // Even indexed points (2)
    QVERIFY(qFuzzyCompare(mUtilities->SimpsonsWeight(2, n), 2.0/3.0));
}
