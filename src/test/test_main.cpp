#include <QTest>
#include "test_calculate.h"
#include "test_parameters.h"
#include "test_miesimulation.h"
#include "test_utilities.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    int status = 0;

    status |= QTest::qExec(new TestParameters(), argc, argv);
    status |= QTest::qExec(new TestCalculate(), argc, argv);
    status |= QTest::qExec(new TestMieSimulation(), argc, argv);
    status |= QTest::qExec(new TestUtilities(), argc, argv);

    return status;
}
