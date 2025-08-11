#ifndef TEST_MIESIMULATION_H
#define TEST_MIESIMULATION_H

#include <QObject>
#include <QTest>
#include <complex>

class MieSimulation;

class TestMieSimulation : public QObject
{
    Q_OBJECT

public:
    TestMieSimulation();
    ~TestMieSimulation();

private slots:
    void init();
    void cleanup();

    void test_FarFieldSolutionForRealRefIndex_sanityCheck();
    void test_FarFieldSolutionForComplexRefIndex_sanityCheck();
    void test_Consistency_RealAndComplex();
    void test_FarFieldSolution_EdgeCases();

private:
    MieSimulation *mMieSimulation;
};

#endif // TEST_MIESIMULATION_H
