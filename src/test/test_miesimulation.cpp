/**********************************************************************
** Tests for functions in miesimulation.cpp
**********************************************************************/

#include <QtMath>
#include "test_miesimulation.h"
#include "calc/miesimulation.h"

TestMieSimulation::TestMieSimulation()
{
}

TestMieSimulation::~TestMieSimulation()
{
}

//Initialize values for tests
void TestMieSimulation::init()
{
    mMieSimulation = new MieSimulation();
}

// clean up resources
void TestMieSimulation::cleanup()
{
    delete mMieSimulation;
    mMieSimulation = nullptr;
}

// Test case: Sanity check for FarFieldSolutionForRealRefIndex.
void TestMieSimulation::test_FarFieldSolutionForRealRefIndex_sanityCheck()
{
    std::complex<double> cS1, cS2;
    double qSca, qExt, qBack;
    double xPara = 10.0;
    double relRef = 1.33;
    double mu = 0.5;

    mMieSimulation->FarFieldSolutionForRealRefIndex(&cS1, &cS2, &qSca, &qExt, &qBack, xPara, relRef, mu);

    // Verify that the output values are not NaN or infinity and are non-zero.
    QVERIFY(!qIsNaN(cS1.real()) && !qIsInf(cS1.real()));
    QVERIFY(!qIsNaN(cS1.imag()) && !qIsInf(cS1.imag()));
    QVERIFY(!qIsNaN(cS2.real()) && !qIsInf(cS2.real()));
    QVERIFY(!qIsNaN(cS2.imag()) && !qIsInf(cS2.imag()));
    QVERIFY(qSca > 0.0);
    QVERIFY(qExt > 0.0);
    QVERIFY(qBack > 0.0);
}

// Test case: Sanity check for FarFieldSolutionForComplexRefIndex.
void TestMieSimulation::test_FarFieldSolutionForComplexRefIndex_sanityCheck()
{
    std::complex<double> cS1, cS2;
    double qSca, qExt, qBack;
    double xPara = 10.0;
    std::complex<double> cRelRef(1.5, -0.1); // Complex refractive index
    double mu = 0.5;

    mMieSimulation->FarFieldSolutionForComplexRefIndex(&cS1, &cS2, &qSca, &qExt, &qBack, xPara, cRelRef, mu);

    // Verify that the output values are not NaN or infinity and are non-zero.
    QVERIFY(!qIsNaN(cS1.real()) && !qIsInf(cS1.real()));
    QVERIFY(!qIsNaN(cS1.imag()) && !qIsInf(cS1.imag()));
    QVERIFY(!qIsNaN(cS2.real()) && !qIsInf(cS2.real()));
    QVERIFY(!qIsNaN(cS2.imag()) && !qIsInf(cS2.imag()));
    QVERIFY(qSca > 0.0);
    QVERIFY(qExt > 0.0);
    QVERIFY(qBack > 0.0);
}

// Test case: Consistency check between real and complex functions.
void TestMieSimulation::test_Consistency_RealAndComplex()
{
    std::complex<double> realcS1, realcS2;
    double realqSca, realqExt, realqBack;

    std::complex<double> complexcS1, complexcS2;
    double complexqSca, complexqExt, complexqBack;

    double xPara = 5.0;
    double relRef = 1.4;
    std::complex<double> crelRef(relRef, 0.0); // Complex refractive index with zero imaginary part
    double mu = -0.9;

    // Call both functions with equivalent parameters
    mMieSimulation->FarFieldSolutionForRealRefIndex(&realcS1, &realcS2, &realqSca, &realqExt, &realqBack, xPara, relRef, mu);
    mMieSimulation->FarFieldSolutionForComplexRefIndex(&complexcS1, &complexcS2, &complexqSca, &complexqExt, &complexqBack, xPara, crelRef, mu);

    // Compare the results. QCOMPARE is a safe way to compare floating-point numbers.
    QCOMPARE(realcS1.real(), complexcS1.real());
    QCOMPARE(realcS1.imag(), complexcS1.imag());
    QCOMPARE(realcS2.real(), complexcS2.real());
    QCOMPARE(realcS2.imag(), complexcS2.imag());
    QCOMPARE(realqSca, complexqSca);
    QCOMPARE(realqExt, complexqExt);
    QCOMPARE(realqBack, complexqBack);
}

// Test case: Edge cases like extreme xPara values.
void TestMieSimulation::test_FarFieldSolution_EdgeCases()
{
    std::complex<double> cS1, cS2;
    double qSca, qExt, qBack;
    double xPara_small = 0.01; // Rayleigh scattering regime
    double xPara_large = 100.0; // Large size parameter
    double relRef = 1.5;
    double mu = 1.0;

    // Test with a very small xPara
    mMieSimulation->FarFieldSolutionForRealRefIndex(&cS1, &cS2, &qSca, &qExt, &qBack, xPara_small, relRef, mu);
    QVERIFY(!qIsNaN(qSca) && qSca > 0);
    QVERIFY(!qIsNaN(qBack) && qBack > 0);

    // Test with a very large xPara
    mMieSimulation->FarFieldSolutionForRealRefIndex(&cS1, &cS2, &qSca, &qExt, &qBack, xPara_large, relRef, mu);
    QVERIFY(!qIsNaN(qSca) && qSca > 0);
    QVERIFY(!qIsNaN(qBack) && qBack > 0);

    // Test for mu = 1 (pure forward scattering)
    mMieSimulation->FarFieldSolutionForRealRefIndex(&cS1, &cS2, &qSca, &qExt, &qBack, 10.0, 1.33, 1.0);
    // For mu = 1, S1 and S2 should be equal
    QCOMPARE(cS1.real(), cS2.real());
    QCOMPARE(cS1.imag(), cS2.imag());

    // Test for mu = -1 (pure backward scattering)
    mMieSimulation->FarFieldSolutionForRealRefIndex(&cS1, &cS2, &qSca, &qExt, &qBack, 10.0, 1.33, -1.0);
    // For mu = -1, S1 and S2 should be related by S1 = -S2
    QCOMPARE(cS1.real(), -cS2.real());
    QCOMPARE(cS1.imag(), -cS2.imag());
}

