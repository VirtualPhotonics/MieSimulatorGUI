#ifndef MAINWINDOWSUPPORT_H
#define MAINWINDOWSUPPORT_H

#include <qmath.h>
#include "lib/qcustomplot.h"
#include "parameters.h"
#include "calc/calculate.h"
#include "ui_mainwindow.h"

class MainWindowSupport
{
public:
    MainWindowSupport(void);
    void SetWidgets(Ui_MainWindow *ui, parameters *para);
    void InitializeGUI(Ui_MainWindow *ui, parameters *para);
    void LoadInputData(Ui_MainWindow *ui, parameters *para);
    void InitializeArrays(Ui_MainWindow *ui, parameters *para, bool *arrayFlag);
    void DeleteArrays(parameters *para, bool *arrayFlag);
    void SetWavelengthSliders(Ui_MainWindow *ui);
    void ProcessMonoDisperse(Ui_MainWindow *ui, parameters *para);
    void ProcessPolyDisperse(Ui_MainWindow *ui, parameters *para);
    void ProcessDistribution(Ui_MainWindow *ui, parameters *para, unsigned int distIndex);
    void DisableEnableRealImagButtons(Ui_MainWindow *ui);
    void DisableWidgetsDuringSimulation(Ui_MainWindow *ui, parameters *para, bool flag);
    void DisableWidgetsDuringCustomPolyDisperseData(Ui_MainWindow *ui, bool flag);
    void ReadCustomData(parameters * para, QString fileName, bool * dataValidFlag);

private:
    QCustomPlot *mCustomPlot;    
    calculate *mCalc;
    QPen mPen;
};

#endif // MAINWINDOWUPPORT_H
