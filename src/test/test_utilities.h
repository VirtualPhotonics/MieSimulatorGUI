#ifndef TEST_UTILITIES_H
#define TEST_UTILITIES_H

#include <QObject>
#include <QtTest>
#include "../app/calc/utilities.h"

class TestUtilities : public QObject
{
    Q_OBJECT

public:
    TestUtilities();
    ~TestUtilities();

private slots:
    void init();
    void cleanup();

    void test_ComplexAbsSquared();
    void test_ComplexAbs();
    void test_SimpsonsWeight();

private:
    utilities *mUtilities;
};

#endif // TEST_UTILITIES_H
