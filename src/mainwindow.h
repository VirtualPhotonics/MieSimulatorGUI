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
#include "displaydialog.h"
#include "mainwindowsupport.h"
#include "plotdata.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushButton_showDistributionAndCustom_clicked();
    void on_pushButton_runSimulation_clicked();
    void on_pushButton_saveData_clicked();
    void on_pushButton_displayData_clicked();
    void on_pushButton_bestFit_clicked();
    void on_pushButton_close_clicked();
    void on_radioButton_linearYaxis_clicked();
    void on_radioButton_logYaxis_clicked();
    void on_radioButton_linearXaxis_clicked();
    void on_radioButton_logXaxis_clicked();
    void on_radioButton_monoDisperse_clicked();
    void on_radioButton_polyDisperse_clicked();
    void on_radioButton_conc_mm3_clicked();
    void on_radioButton_volFrac_clicked();

    void on_radioButton_polarPlotScale_linear_clicked();
    void on_radioButton_polarPlotScale_log_clicked();
    void on_radioButton_pFunction_ave_clicked();
    void on_radioButton_pFunction_para_clicked();
    void on_radioButton_pFunction_perp_clicked();

    void on_radioButton_s1_amp_clicked();
    void on_radioButton_s2_amp_clicked();
    void on_radioButton_s1s2_abs_clicked();
    void on_radioButton_s1s2_real_clicked();
    void on_radioButton_s1s2_imag_clicked();
    void on_radioButton_dThetaStep0_1_clicked();
    void on_doubleSpinBox_powerLaw_f_valueChanged(double arg1);
    void on_doubleSpinBox_powerLaw_b_valueChanged(double arg1);
    void on_slider_concPercentChange_valueChanged(int value);
    void on_slider_pFunctionPolar_wavel_valueChanged(int value);
    void on_slider_pFunctionLinear_wavel_valueChanged(int value);
    void on_slider_s1s2_wavel_valueChanged(int value);
    void on_comboBox_distribution_currentIndexChanged(int value);

    void mouseMove_customPlot_csca(QMouseEvent *event);
    void mouseMove_customPlot_cext(QMouseEvent *event);
    void mouseMove_customPlot_cback(QMouseEvent *event);
    void mouseMove_customPlot_sizePara(QMouseEvent *event);
    void mouseMove_customPlot_mus(QMouseEvent *event);
    void mouseMove_customPlot_musp(QMouseEvent *event);
    void mouseMove_customPlot_g(QMouseEvent *event);
    void mouseMove_customPlot_forwardBackward(QMouseEvent *event);
    void mouseMove_customPlot_s1s2(QMouseEvent *event);
    void mouseMove_customPlot_pFunctionLinear(QMouseEvent *event);
    void mouseMove_customPlot_distribution(QMouseEvent *event);
    void mouseMove_customPlot_muspPowerLaw(QMouseEvent *event);

private:
    Ui::MainWindow *ui;
    parameters *_para;
    bool _distPlotFlag;            //Flag to check distributionPlot
    bool _otherPlotsFlag;          //Flag to check other customplots
    bool _arrayFlag;               //Flag to find dynamic array status to initialize or delete
    bool _loadCustomDataCheckFlag; //Flag to check loading "Custom" data

    void Initialize();
    void UpdatePhaseFunctionPolarPlot();
    void UpdatePhaseFunctionLinearPlot();
    void UpdateS1S2Plot();
    void UpdateMuspFitPlot();
    void UpdateMuspFitErrorDisplay();
    void DisplayCurveData(QMouseEvent *event, QCustomPlot *curPlot, QString strNameX, QString strNameY);
};
#endif // MAINWINDOW_H
