#ifndef MAINWINDOWSUPPORT_H
#define MAINWINDOWSUPPORT_H

#include <qmath.h>
#include "parameters.h"
#include "calc/calculate.h"
#include "ui_mainwindow.h"

class MainWindowSupport
{
public:
    MainWindowSupport(void);
    void InitializeGUI(Ui_MainWindow *ui, Parameters *para);
    void SetWidgets(Ui_MainWindow *ui, Parameters *para);
    void LoadInputData(Ui_MainWindow *ui, Parameters *para);
    void InitializeArrays(Ui_MainWindow *ui, Parameters *para, bool *arrayFlag);
    void DeleteArrays(Parameters *para, bool *arrayFlag);
    void SetWavelengthSliders(Ui_MainWindow *ui);
    void ProcessMonoDisperse(Ui_MainWindow *ui, Parameters *para);
    void ProcessPolyDisperse(Ui_MainWindow *ui, Parameters *para);
    void ProcessDistribution(Ui_MainWindow *ui, Parameters *para, unsigned int distIndex);
    void DisableEnableRealImagButtons(Ui_MainWindow *ui);
    void DisableWidgetsDuringSimulation(Ui_MainWindow *ui, Parameters *para, bool flag);
    void DisableWidgetsDuringCustomPolyDisperseData(Ui_MainWindow *ui, bool flag);
    void ReadCustomData(Parameters * para, QString fileName, bool * dataValidFlag);
    void PrepareScatteringRegimeWarning(double clearanceToWavelength, double sizeParameter,
                                        double volFraction, double criticalWavelength,
                                        QString strRegime);
    void DisplayWarning(QString warningMessage);

private:
    Calculate *mCalc;
    QPen mPen;
};

#endif // MAINWINDOWUPPORT_H
