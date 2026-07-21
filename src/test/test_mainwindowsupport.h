#ifndef TEST_MAINWINDOWSUPPORT_H
#define TEST_MAINWINDOWSUPPORT_H

#include <QObject>
#include "parameters.h"
#include "dialog/mainwindowsupport.h"

class TestMainWindowSupport : public QObject
{
    Q_OBJECT

public:
    TestMainWindowSupport() = default;;
    ~TestMainWindowSupport() = default;;

private slots:
    void init();
    void cleanup();

    void test_ReadCustomData_validFourColumn();
    void test_ReadCustomData_validFiveColumn();
    void test_ReadCustomData_delimiterVariants();
    void test_DeleteArrays_resetsArrayFlag();

private:
    Parameters *mPara;
    MainWindowSupport *mSupport;
};

#endif // TEST_MAINWINDOWSUPPORT_H