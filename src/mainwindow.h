#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMessageBox>
#include <qmath.h>
#include "qwt/qwt_polar_grid.h"
#include "qwt/qwt_polar_curve.h"
#include "qwt/qwt_scale_engine.h"
#include "qcustomplot.h"
#include "parameters.h"
#include "calculate.h"
#include "utilities.h"
#include "optionsdialog.h"
#include "ui_mainwindow.h"
#include "displaydialog.h"
#include "mainwindowsupport.h"
#include "plotdata.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:    
    void on_pushButton_ShowDistributionAndCustom_clicked();
    void on_pushButton_RunSimulation_clicked();
    void on_pushButton_SaveData_clicked();
    void on_pushButton_DisplayData_clicked();
    void on_pushButton_BestFit_clicked();
    void on_pushButton_Close_clicked();
    void on_radioButton_LinearYAxis_clicked();
    void on_radioButton_LogYAxis_clicked();
    void on_radioButton_LinearXAxis_clicked();
    void on_radioButton_LogXAxis_clicked();
    void on_radioButton_MonoDisperse_clicked();
    void on_radioButton_PolyDisperse_clicked();
    void on_radioButton_Conc_mm3_clicked();
    void on_radioButton_VolFrac_clicked();
    void on_radioButton_PhaseLinear_clicked();
    void on_radioButton_PhaseLog_clicked();
    void on_radioButton_PhaseAverage_clicked();
    void on_radioButton_PhasePara_clicked();
    void on_radioButton_PhasePerp_clicked();
    void on_radioButton_S1_clicked();
    void on_radioButton_S2_clicked();
    void on_radioButton_S1S2_Abs_clicked();
    void on_radioButton_S1S2_Real_clicked();
    void on_radioButton_S1S2_Imag_clicked();
    void on_radioButton_Phase_DTheta0_1_clicked();
    void on_doubleSpinBox_F_valueChanged(double arg1);
    void on_doubleSpinBox_B_valueChanged(double arg1);
    void on_slider_ConcPercentChange_valueChanged(int value);
    void on_slider_WL_PFPolar_valueChanged(int value);
    void on_slider_WL_PFLinear_valueChanged(int value);
    void on_slider_WL_S1S2_valueChanged(int value);
    void on_comboBox_Distribution_currentIndexChanged(int value);

    void Initialize();
    void UpdatePhaseFunctionPolarPlot();
    void UpdatePhaseFunctionLinearPlot();
    void UpdateS1S2Plot();
    void UpdateMuspFitPlot();
    void UpdateMuspFitError();
    void MouseOverPlotScatteringCrossSection(QMouseEvent *event);
    void MouseOverPlotExtinctionCrossSection(QMouseEvent *event);
    void MouseOverPlotBackscatteringCrossSection(QMouseEvent *event);
    void MouseOverPlotSizeParameter(QMouseEvent *event);
    void MouseOverPlotMus(QMouseEvent *event);
    void MouseOverPlotMusp(QMouseEvent *event);
    void MouseOverPlotG(QMouseEvent *event);
    void MouseOverPlotForwardBackward(QMouseEvent *event);
    void MouseOverPlotS1S2(QMouseEvent *event);
    void MouseOverPlotPhaseFunctionLinear(QMouseEvent *event);
    void MouseOverPlotDistribution(QMouseEvent *event);
    void MouseOverPlotMuspPowerLaw(QMouseEvent *event);
    void DisplayCurveData(QMouseEvent *event, QCustomPlot *curPlot, QString strNameX, QString strNameY);

private:
    Ui::MainWindow *ui;
    parameters *mPara;
    bool mDistPlotFlag;         //Flag to check distributionPlot
    bool mOtherPlotsFlag;       //Flag to check other customplots
    bool mArrayFlag;            //Flag to find dynamic array status to initialize or delete
    bool mLoadCustomNoGoodFlag; //Flag to check "Custom" data
};

#endif // MAINWINDOW_H
