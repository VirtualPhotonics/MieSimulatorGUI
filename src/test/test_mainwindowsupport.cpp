/**********************************************************************
** Tests for functions in dialog/mainwindowsupport.cpp
**********************************************************************/

#include <QtMath>
#include <QTest>
#include <QTemporaryDir>
#include <QFile>
#include <QTextStream>
#include "test_mainwindowsupport.h"

TestMainWindowSupport::TestMainWindowSupport()
{
}

TestMainWindowSupport::~TestMainWindowSupport()
{
}

//Initialize values for tests
void TestMainWindowSupport::init()
{
    mPara = new Parameters();
    mSupport = new MainWindowSupport();
}

// clean up resources
void TestMainWindowSupport::cleanup()
{
    // ReadCustomData allocates these directly on mPara; Parameters' own
    // destructor does not free them, so the test must do so explicitly.
    delete [] mPara->radArray;
    delete [] mPara->numDensityArray;
    delete [] mPara->scatRefRealArray;
    delete [] mPara->scatRefImagArray;
    delete [] mPara->medRefArray;

    delete mPara;
    mPara = nullptr;
    delete mSupport;
    mSupport = nullptr;
}

// Test case: ReadCustomData with a valid 4-column file (diameter, numDensity,
// scatRefReal, scatRefImag). medRefArray should fall back to para->medRef.
void TestMainWindowSupport::test_ReadCustomData_validFourColumn()
{
    QTemporaryDir dir;
    QVERIFY(dir.isValid());
    QString fileName = dir.filePath("custom_four_column.csv");

    QFile file(fileName);
    QVERIFY(file.open(QIODevice::WriteOnly | QIODevice::Text));
    QTextStream out(&file);
    out << "2.0,100,1.4,0.0\n";
    out << "4.0,200,1.5,0.01\n";
    out << "6.0,150,1.45,0.0\n";
    file.close();

    mPara->medRef = 1.333;
    bool dataValidFlag = false;
    mSupport->ReadCustomData(mPara, fileName, &dataValidFlag);

    QCOMPARE(mPara->nRadius, 3u);
    // Diameters are halved into radii
    QCOMPARE(mPara->radArray[0], 1.0);
    QCOMPARE(mPara->radArray[1], 2.0);
    QCOMPARE(mPara->radArray[2], 3.0);
    QCOMPARE(mPara->numDensityArray[1], 200.0);
    QCOMPARE(mPara->scatRefRealArray[2], 1.45);
    QCOMPARE(mPara->scatRefImagArray[1], 0.01);
    // 4-column format has no per-row medium ref index -> defaults to para->medRef
    QCOMPARE(mPara->medRefArray[0], 1.333);
    QCOMPARE(mPara->medRefArray[2], 1.333);

    QCOMPARE(mPara->minRadius, 1.0);
    QCOMPARE(mPara->maxRadius, 3.0);
    QCOMPARE(mPara->meanRadius, 2.0);
    QCOMPARE(mPara->stdDev, 1.0);
}

// Test case: ReadCustomData with a valid 5-column file, where the 5th column
// supplies a per-row medium refractive index instead of the default.
void TestMainWindowSupport::test_ReadCustomData_validFiveColumn()
{
    QTemporaryDir dir;
    QVERIFY(dir.isValid());
    QString fileName = dir.filePath("custom_five_column.csv");

    QFile file(fileName);
    QVERIFY(file.open(QIODevice::WriteOnly | QIODevice::Text));
    QTextStream out(&file);
    out << "2.0,100,1.4,0.0,1.330\n";
    out << "8.0,200,1.5,0.02,1.335\n";
    file.close();

    mPara->medRef = 1.333; // should NOT be used when a 5th column is present
    bool dataValidFlag = false;
    mSupport->ReadCustomData(mPara, fileName, &dataValidFlag);

    QCOMPARE(mPara->nRadius, 2u);
    QCOMPARE(mPara->radArray[0], 1.0);
    QCOMPARE(mPara->radArray[1], 4.0);
    QCOMPARE(mPara->medRefArray[0], 1.330);
    QCOMPARE(mPara->medRefArray[1], 1.335);

    QCOMPARE(mPara->minRadius, 1.0);
    QCOMPARE(mPara->maxRadius, 4.0);
    QCOMPARE(mPara->meanRadius, 2.5);
}

// Test case: ReadCustomData accepts comma, semicolon, and tab as delimiters.
void TestMainWindowSupport::test_ReadCustomData_delimiterVariants()
{
    QTemporaryDir dir;
    QVERIFY(dir.isValid());
    QString fileName = dir.filePath("custom_delimiters.csv");

    QFile file(fileName);
    QVERIFY(file.open(QIODevice::WriteOnly | QIODevice::Text));
    QTextStream out(&file);
    out << "2.0,100,1.4,0.0\n";   // comma
    out << "4.0;200;1.5;0.0\n";   // semicolon
    out << "6.0\t300\t1.6\t0.0\n"; // tab
    file.close();

    bool dataValidFlag = false;
    mSupport->ReadCustomData(mPara, fileName, &dataValidFlag);

    QCOMPARE(mPara->nRadius, 3u);
    QCOMPARE(mPara->radArray[0], 1.0);
    QCOMPARE(mPara->radArray[1], 2.0);
    QCOMPARE(mPara->radArray[2], 3.0);
    QCOMPARE(mPara->numDensityArray[1], 200.0);
    QCOMPARE(mPara->numDensityArray[2], 300.0);
}

// Test case: DeleteArrays frees all owned arrays and resets the array flag,
// without touching pointers it doesn't own (e.g. radArray/numDensityArray).
void TestMainWindowSupport::test_DeleteArrays_resetsArrayFlag()
{
    mPara->nWavel = 2;

    mPara->wavelArray = new double[2]{600.0, 700.0};
    mPara->cSca = new double[2]{0.0, 0.0};
    mPara->cExt = new double[2]{0.0, 0.0};
    mPara->cBack = new double[2]{0.0, 0.0};
    mPara->sizePara = new double[2]{0.0, 0.0};
    mPara->mus = new double[2]{0.0, 0.0};
    mPara->g = new double[2]{0.0, 0.0};
    mPara->forward = new double[2]{0.0, 0.0};
    mPara->backward = new double[2]{0.0, 0.0};

    mPara->phaseFunctionAve = new double*[2];
    mPara->phaseFunctionPara = new double*[2];
    mPara->phaseFunctionPerp = new double*[2];
    mPara->S1 = new std::complex<double>*[2];
    mPara->S2 = new std::complex<double>*[2];
    for (unsigned int i = 0; i < mPara->nWavel; i++)
    {
        mPara->phaseFunctionAve[i] = new double[3]{0.0, 0.0, 0.0};
        mPara->phaseFunctionPara[i] = new double[3]{0.0, 0.0, 0.0};
        mPara->phaseFunctionPerp[i] = new double[3]{0.0, 0.0, 0.0};
        mPara->S1[i] = new std::complex<double>[3];
        mPara->S2[i] = new std::complex<double>[3];
    }

    bool arrayFlag = true;
    mSupport->DeleteArrays(mPara, &arrayFlag);

    QVERIFY(!arrayFlag);

    // DeleteArrays only frees each row of the 2D arrays, not the outer
    // pointer arrays themselves; free those here to keep the test clean.
    delete [] mPara->phaseFunctionAve;
    delete [] mPara->phaseFunctionPara;
    delete [] mPara->phaseFunctionPerp;
    delete [] mPara->S1;
    delete [] mPara->S2;
    mPara->phaseFunctionAve = nullptr;
    mPara->phaseFunctionPara = nullptr;
    mPara->phaseFunctionPerp = nullptr;
    mPara->S1 = nullptr;
    mPara->S2 = nullptr;
}