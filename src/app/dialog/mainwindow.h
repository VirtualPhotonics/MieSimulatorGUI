#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMessageBox>
#include <qmath.h>
#include "lib/qcustomplot.h"
#include "parameters.h"
#include "ui_mainwindow.h"


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
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
    void on_radioButton_NumDen_clicked();
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
    void on_radioButton_FittingSimple_clicked();
    void on_radioButton_FittingComplex_clicked();
    void on_radioButton_RefWavel500_clicked();
    void on_radioButton_RefWavel600_clicked();
    void on_radioButton_RefWavel700_clicked();
    void on_radioButton_RefWavel800_clicked();
    void on_radioButton_RefWavel900_clicked();
    void on_radioButton_RefWavel1000_clicked();
    void on_slider_ConcPercentChange_valueChanged(int value);
    void on_slider_WL_PFPolar_valueChanged(int value);
    void on_slider_WL_PFLinear_valueChanged(int value);
    void on_slider_WL_S1S2_valueChanged(int value);
    void on_comboBox_Distribution_currentIndexChanged(int value);
    void on_pushButton_SavePlot_Mus_clicked();
    void on_pushButton_SavePlot_Csca_clicked();
    void on_pushButton_SavePlot_Cext_clicked();
    void on_pushButton_SavePlot_Cback_clicked();
    void on_pushButton_SavePlot_Musp_clicked();
    void on_pushButton_SavePlot_MuspPowerLaw_clicked();
    void on_pushButton_SavePlot_Distribution_clicked();
    void on_pushButton_SavePlot_SizePara_clicked();
    void on_pushButton_SavePlot_S1S2_clicked();
    //void on_pushButton_SavePlot_PhaseFunctionPolar_clicked();
    void on_pushButton_SavePlot_PhaseFunctionLinear_clicked();
    void on_pushButton_SavePlot_G_clicked();
    void on_pushButton_SavePlot_FB_clicked();

private:
    void Initialize();
    void UpdatePhaseFunctionPolarPlot();
    void UpdatePhaseFunctionLinearPlot();
    void UpdateS1S2Plot();
    void UpdateMuspFitPlot();
    void UpdateMuspFitErrorDisplay();
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
    void SavePlot(QCustomPlot *curPlot, QString fileName);
    //void SavePolarPlot(QwtPolarPlot *curPlot, QString fileName);
    void RememberLastDirectory(QString fileName);    
    void Slider_F_valueChanged(int value);
    void Slider_B_valueChanged(int value);
    void DoubleSpinBox_F_valueChanged (double value);
    void DoubleSpinBox_B_valueChanged (double value);
    void CalculateFittingAndDisplay();

private:
    Ui::MainWindow *ui;
    parameters *mPara;
    bool mDistPlotFlag;         //Flag to check distributionPlot
    bool mOtherPlotsFlag;       //Flag to check other customplots
    bool mArrayFlag;            //Flag to find dynamic array status to initialize or delete
    bool mLoadCustomNoGoodFlag; //Flag to check "Custom" data
};

#endif // MAINWINDOW_H
