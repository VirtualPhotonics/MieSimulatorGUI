//**********************************************************************
//** All MainWindow functions are listed in this file.
//**********************************************************************

#include "mainwindow.h"

//***********************Starting Main window and initializations ***********************
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent), ui(new Ui::MainWindow)
{
    //Initialize UI
    ui->setupUi(this);

   //Intialize Window
    Initialize();

    //Connect signals and slots
    connect(ui->customPlot_csca, SIGNAL(mouseMove(QMouseEvent*)), SLOT(mouseMove_customPlot_csca(QMouseEvent*)));
    connect(ui->customPlot_cext, SIGNAL(mouseMove(QMouseEvent*)), SLOT(mouseMove_customPlot_cext(QMouseEvent*)));
    connect(ui->customPlot_cback, SIGNAL(mouseMove(QMouseEvent*)), SLOT(mouseMove_customPlot_cback(QMouseEvent*)));
    connect(ui->customPlot_sizePara, SIGNAL(mouseMove(QMouseEvent*)), SLOT(mouseMove_customPlot_sizePara(QMouseEvent*)));
    connect(ui->customPlot_mus, SIGNAL(mouseMove(QMouseEvent*)), SLOT(mouseMove_customPlot_mus(QMouseEvent*)));
    connect(ui->customPlot_musp, SIGNAL(mouseMove(QMouseEvent*)), SLOT(mouseMove_customPlot_musp(QMouseEvent*)));
    connect(ui->customPlot_g, SIGNAL(mouseMove(QMouseEvent*)), SLOT(mouseMove_customPlot_g(QMouseEvent*)));
    connect(ui->customPlot_forwardBackward, SIGNAL(mouseMove(QMouseEvent*)), SLOT(mouseMove_customPlot_forwardBackward(QMouseEvent*)));
    connect(ui->customPlot_s1s2, SIGNAL(mouseMove(QMouseEvent*)), SLOT(mouseMove_customPlot_s1s2(QMouseEvent*)));
    connect(ui->customPlot_pFunctionLinear, SIGNAL(mouseMove(QMouseEvent*)), SLOT(mouseMove_customPlot_pFunctionLinear(QMouseEvent*)));
    connect(ui->customPlot_distribution, SIGNAL(mouseMove(QMouseEvent*)), SLOT(mouseMove_customPlot_distribution(QMouseEvent*)));
    connect(ui->customPlot_muspPowerLaw, SIGNAL(mouseMove(QMouseEvent*)), SLOT(mouseMove_customPlot_muspPowerLaw(QMouseEvent*)));
}

MainWindow::~MainWindow()
{
    delete ui;
}

//***************************** Intialize *********************************************
void MainWindow::Initialize()
{
    //Initialize paramaters
    _otherPlotsFlag = false;
    _distPlotFlag = false;
    _arrayFlag = false;
    _loadCustomDataCheckFlag = false;
    _para = new parameters();
    _para->polarCurve = new QwtPolarCurve();

    MainWindowSupport support;
    support.InitializeGUI(ui);
    support.SetWidgets(ui);

    PlotData plot;
    plot.InitializeDistributionPlot(ui);
    plot.InitializePhaseFunctionPolarPlot(ui,_para);
    plot.InitializePhaseFunctionLinearPlot(ui);
    plot.InitializeAllOtherPlots(ui);
}

//********************************** Run Simulation ************************************
void MainWindow::on_pushButton_runSimulation_clicked()
{
    PlotData plot;
    MainWindowSupport support;

    //Initialize
    plot.ClearPlots(ui,_para);
    support.SetWidgets(ui);
    support.LoadInputData(ui,_para);
    support.SetWavelengthSliders(ui);              //set slider position

    if  (_loadCustomDataCheckFlag)
    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Error");
        msgBox.setText("No Data available.");
        msgBox.setInformativeText("Click 'Load Custom Data'");
        msgBox.exec();
    }
    else
    {

        if (!support.CheckInputParameters(ui,_para))    //sanity check
        {
            //Disable widgets
            support.DisableWidgetsDuringSimulation(ui, _para, true);
            //If _arrayFlag is TRUE, delete dynamic arrays before next simulation
            if (_arrayFlag)
                support.DeleteArrays(_para, &_arrayFlag);
            //Initialize dynamic arrays
            support.InitializeArrays(ui,_para, &_arrayFlag);

            //Mono disperse
            if (ui->radioButton_monoDisperse->isChecked())
            {
                support.ProcessDistribution(ui, _para, 3);
                support.ProcessMonoDisperse(ui,_para);
                ui->label_progress->setText("<font color=\"green\"> Completed!</font>");
                ui->slider_concPercentChange->setValue(0);
                _otherPlotsFlag = true;
                _distPlotFlag = true;
            }

            //Poly disperse
            if (ui->radioButton_polyDisperse->isChecked())
            {
                if (!support.CheckDistribution(ui,_para))   //sanity check
                {
                    support.ProcessDistribution(ui, _para, static_cast<unsigned int>(ui->comboBox_distribution->currentIndex()));
                    support.ProcessPolyDisperse(ui,_para);
                    ui->label_progress->setText("<font color=\"green\"> Completed!</font>");
                    ui->slider_concPercentChange->setValue(0);
                    _otherPlotsFlag = true;
                    _distPlotFlag = true;
                }
                if (ui->comboBox_distribution->currentIndex() ==2)
                    support.DisableWidgetsDuringCustomPolyDisperseData(ui, true);
                else
                    support.DisableWidgetsDuringCustomPolyDisperseData(ui, false);
            }
            //Enable widgets
            support.DisableWidgetsDuringSimulation(ui, _para, false);
        }
    }
}

//********************************** Show distribution ****************************************
void MainWindow::on_pushButton_showDistributionAndCustom_clicked()
{
    MainWindowSupport support;
    bool dataValidFlag = false;

    if (ui->comboBox_distribution->currentIndex() == 2)
    {
        MainWindowSupport support;
        QString fileName = QFileDialog::getOpenFileName
          (this,tr("Read Data"),"",tr("Text File (*.txt);;All Files (*)"));

        if (fileName.isEmpty())
            return;
        else
            support.ReadCustomData(_para, fileName, &dataValidFlag);
    }

    if (dataValidFlag)
    {
        PlotData plot;
        plot.ClearPlots(ui,_para);
        _otherPlotsFlag = false;
        _distPlotFlag = false;
        _loadCustomDataCheckFlag = true;
    }
    else
    {
        support.LoadInputData(ui, _para);
        if (!support.CheckInputParameters(ui,_para))    //sanity check
        {
            if (!support.CheckDistribution(ui,_para))   //sanity check
            {
                support.SetWidgets(ui);
                support.ProcessDistribution(ui, _para, static_cast<unsigned int>(ui->comboBox_distribution->currentIndex()));
                _distPlotFlag = true;
            }
        }
        _loadCustomDataCheckFlag = false;
    }
}

//**********************************  Save Data ********************************************
void MainWindow::on_pushButton_saveData_clicked()
{
    if (_otherPlotsFlag)
    {
        OptionsDialog saveOptions;
        saveOptions.SaveData(ui, _para);
    }
    else
    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Error");
        msgBox.setText("No Data available.");
        msgBox.exec();
    }
}

//********************************** Display Data ********************************************
void MainWindow::on_pushButton_displayData_clicked()
{
    if (_otherPlotsFlag)
    {
        DisplayDialog *display;
        display = new DisplayDialog();
        display->DisplayData(ui,  _para);
        display->show();
    }
    else
    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Error");
        msgBox.setText("No Data available.");
        msgBox.exec();
    }
}

//********************************** Calculate Best Fit musp **************************************
void MainWindow::on_pushButton_bestFit_clicked()
{
    if(_otherPlotsFlag)
    {
        if (_para->nWavel >1)
        {
            MainWindowSupport support;
            utilities util;
            calculate cal;
            PlotData plot;
            double fRay, bMie;

            //Widget settings
            support.DisableWidgetsDuringSimulation(ui, _para, true);
            ui->label_currentMse->setText("<font color=\"red\">Wait!...    </font>");
            ui->label_currentA->setText("");
            util.Delay();

            //Calculate fRay and bMie
            if (_para->fittingComplex)
                cal.CalculatePowerLawAutoFitComplex(_para);
            else
                cal.CalculatePowerLawAutoFitSimple(_para);
            plot.AssignValuesMuspPowerLawPlots(ui,_para);
            fRay = _para->fRay;
            bMie = _para->bMie;
            ui->doubleSpinBox_powerLaw_f->setValue(fRay);
            ui->doubleSpinBox_powerLaw_b->setValue(bMie);
            UpdateMuspFitErrorDisplay();

            //Enable Widgets
            support.DisableWidgetsDuringSimulation(ui, _para, false);
        }
        else
        {
            QMessageBox msgBox;
            msgBox.setWindowTitle("Error");
            msgBox.setText("Use two or more wavelengths.");
            msgBox.exec();
        }
    }
}

//********************************** Close GUI ********************************************
void MainWindow::on_pushButton_close_clicked()
{
   this->close();
}

//********************************** Radio Button Selection ********************************

//radioButton_LinearYAxis_clicked: Select Linear Y-Axis
void MainWindow::on_radioButton_linearYaxis_clicked()
{
    PlotData plot;
    if (_otherPlotsFlag)
    {
        MainWindowSupport support;
        support.DisableEnableRealImagButtons(ui);
        plot.AssignValuesS1S2Plot(ui, _para);
        plot.AssignValuesPhaseFunctionLinearPlot(ui, _para);
        plot.AssignValuesAllOtherPlots(ui, _para);
    }
    if (_distPlotFlag)
        plot.AssignValuesDistributionPlot(ui, _para);
}

//radioButton_logYaxis_clicked: Select Log Y-Axis
void MainWindow::on_radioButton_logYaxis_clicked()
{
    PlotData plot;
    if (_otherPlotsFlag)
    {
        MainWindowSupport support;
        support.DisableEnableRealImagButtons(ui);
        plot.AssignValuesS1S2Plot(ui, _para);
        plot.AssignValuesPhaseFunctionLinearPlot(ui, _para);
        plot.AssignValuesAllOtherPlots(ui, _para);
    }
    if (_distPlotFlag)
        plot.AssignValuesDistributionPlot(ui, _para);
}

//radioButton_LinearXAxis_clicked: Select Linear X-Axis
void MainWindow::on_radioButton_linearXaxis_clicked()
{
    if (_distPlotFlag)
    {
        PlotData plot;
        plot.AssignValuesDistributionPlot(ui, _para);
    }
}

//radioButton_LogXAxis_clicked: Select Log X-Axis
void MainWindow::on_radioButton_logXaxis_clicked()
{
    if (_distPlotFlag)
    {
        PlotData plot;
        plot.AssignValuesDistributionPlot(ui, _para);
    }
}

//radioButtonMonoDisperse_clicked: Select Mono disperse distribution
void MainWindow::on_radioButton_monoDisperse_clicked()
{
    MainWindowSupport support;
    ui->comboBox_distribution->setCurrentIndex(0);
    support.SetWidgets(ui);
    PlotData plot;
    plot.ClearPlots(ui,_para);
    ui->tabWidget_pFunction->setTabText(0,"S1/S2");
    ui->tabWidget_pFunction->setTabText(1,"Phase Function (Polar)");
    ui->tabWidget_pFunction->setTabText(2,"Phase Function (Linear)");
    ui->tabWidget_scat->setTabText(1,"Csca");
    ui->tabWidget_scat->setTabText(2,"Cext");
    ui->tabWidget_scat->setTabText(3,"Cback");
    ui->tabWidget_ns->setTabText(0,"Concentration (Number Density: Ns)");
    ui->tabWidget_ns->setTabEnabled(1,true);
    ui->tabWidget_ns->setTabText(1,"Size Parameter");
    _otherPlotsFlag = false;
    _distPlotFlag = false;
}

//radioButtonPolyDisperse_clicked: Select Poly disperse distribution
void MainWindow::on_radioButton_polyDisperse_clicked()
{
    MainWindowSupport support;
    support.SetWidgets(ui);
    PlotData plot;
    plot.ClearPlots(ui,_para);
    ui->tabWidget_pFunction->setTabText(0,"Ave. S1/S2");
    ui->tabWidget_pFunction->setTabText(1,"Ave. Phase Function (Polar)");
    ui->tabWidget_pFunction->setTabText(2,"Ave. Phase Function (Linear)");
    ui->tabWidget_scat->setTabText(1,"Ave. Csca");
    ui->tabWidget_scat->setTabText(2,"Ave. Cext");
    ui->tabWidget_scat->setTabText(3,"Ave. Cback");
        ui->tabWidget_ns->setTabText(0,"Scatterer Distribution and Concentration (Ns)");
    ui->tabWidget_ns->setTabEnabled(1,false);
    ui->tabWidget_ns->setTabText(1,"");
    _otherPlotsFlag = false;
    _distPlotFlag = false;
}

//radioButton_VolFrac_clicked: Select Volume Fraction (sphere volume/cubic mm)
void MainWindow::on_radioButton_volFrac_clicked()
{
    _para->volFraction = ui->lineEdit_volFrac->text().toDouble();
    ui->lineEdit_conc_mm3->setText(QString::number(_para->sphNumDensity));

    MainWindowSupport support;
    support.SetWidgets(ui);
    PlotData plot;
    plot.ClearPlots(ui,_para);
    _otherPlotsFlag = false;
    _distPlotFlag = false;
}

//radioButton_Conc_mm3_clicked: Select Concentraction (number of spheres in 1 cubic mm)
void MainWindow::on_radioButton_conc_mm3_clicked()
{
    ui->lineEdit_volFrac->setText(QString::number(_para->volFraction));
    _para->sphNumDensity = ui->lineEdit_conc_mm3->text().toDouble();

    MainWindowSupport support;
    support.SetWidgets(ui);
    PlotData plot;
    plot.ClearPlots(ui,_para);
    _otherPlotsFlag = false;
    _distPlotFlag = false;
}

//radioButtonS1_clicked: Select S1 plot
void MainWindow::on_radioButton_s1_amp_clicked()
{
    UpdateS1S2Plot();
}

//radioButtonS2_clicked: Select S2 plot
void MainWindow::on_radioButton_s2_amp_clicked()
{
    UpdateS1S2Plot();
}

//radioButton_s1_ampS2_Abs_clicked: Select absolute value
void MainWindow::on_radioButton_s1s2_abs_clicked()
{
    UpdateS1S2Plot();
}

//radioButton_s1_ampS2_Real_clicked: Select Real component
void MainWindow::on_radioButton_s1s2_real_clicked()
{
    UpdateS1S2Plot();
}

//radioButton_s1_ampS2_Imag_clicked: Select Imaginary component
void MainWindow::on_radioButton_s1s2_imag_clicked()
{
    UpdateS1S2Plot();
}

//radioButton_PhaseLinear_clicked: Select Linear Scale
void MainWindow::on_radioButton_polarPlotScale_linear_clicked()
{
    UpdatePhaseFunctionPolarPlot();
}

//radioButton_PhaseLog_clicked: Select Log Scale
void MainWindow::on_radioButton_polarPlotScale_log_clicked()
{
    UpdatePhaseFunctionPolarPlot();
}

//radioButton_PhaseAverage_clicked: Select average
void MainWindow::on_radioButton_pFunction_ave_clicked()
{
    UpdatePhaseFunctionPolarPlot();
}

//radioButton_PhasePara_clicked: Select Parallel component
void MainWindow::on_radioButton_pFunction_para_clicked()
{
    UpdatePhaseFunctionPolarPlot();
}

//radioButton_PhasePerp_clicked: Select Perpendicular component
void MainWindow::on_radioButton_pFunction_perp_clicked()
{
    UpdatePhaseFunctionPolarPlot();
}

//radioButton_dTheta0_1 is clicked: Simulation time warning
void MainWindow:: on_radioButton_dThetaStep0_1_clicked()
{
    QMessageBox msgBox;

    msgBox.setIcon(QMessageBox::Warning);
    msgBox.setWindowTitle("Warning");
    msgBox.setText("This selection will increase the simulation time.");
    msgBox.setInformativeText("Do you want to continue?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::No);
    int ret = msgBox.exec();

    switch (ret)
    {
      case QMessageBox::Yes:
          ui->radioButton_dThetaStep0_1->setChecked(true);
          break;
      case QMessageBox::No:
          ui->radioButton_dThetaStep0_5->setChecked(true);
          break;
    }
}

//********************************** Slider Selection **********************************

//sliderConcPercentChange_slidervalueChanged: Change Conc or VolFrac
void MainWindow::on_slider_concPercentChange_valueChanged(int position)
{
    double value;
    double margin = (1.0 + position /200.0);
    ui->label_actualConcPercent->setText(QString::number(position/2.0)+"%");

    if (ui->radioButton_conc_mm3->isChecked())
    {
        value =  _para->sphNumDensity * margin;
        ui->lineEdit_conc_mm3->setText(QString::number(value));
    }
    if (ui->radioButton_volFrac->isChecked())
    {
        value =  _para->volFraction * margin;
        ui->lineEdit_volFrac->setText(QString::number(value));
    }

    PlotData plot;
    if (_distPlotFlag)
        plot.AssignValuesDistributionPlot(ui, _para);
    if (_otherPlotsFlag)
    {   MainWindowSupport support;
        support.DisableEnableRealImagButtons(ui);
        plot.AssignValuesAllOtherPlots(ui, _para);
    }
}

//on_slider_WL_PFPolar_valueChanged: Change Phase Function Wavelenghth Slider
void MainWindow::on_slider_pFunctionPolar_wavel_valueChanged(int value)
{
   int wavel = static_cast<int>(_para->startWavel + value*_para->stepWavel);
   ui->label_pFunctionPolar_curWavel->setText(QString::number(wavel));
   UpdatePhaseFunctionPolarPlot();
}

//on_slider_WL_PFLinear_valueChanged: Change Phase Function Wavelenghth Slider
void MainWindow::on_slider_pFunctionLinear_wavel_valueChanged(int value)
{
   int wavel = static_cast<int>(_para->startWavel + value*_para->stepWavel);
   ui->label_pFunctionLinear_curWavel->setText(QString::number(wavel));
   UpdatePhaseFunctionLinearPlot();
}

//on_slider_WL_S1S2_valueChanged: Change S1/S2 Wavelenghth Slider
void MainWindow::on_slider_s1s2_wavel_valueChanged(int value)
{
   int wavel = static_cast<int>(_para->startWavel + value*_para->stepWavel);
   ui->label_s1s2_curWavel->setText(QString::number(wavel));
   UpdateS1S2Plot();
}

//SpinBoxF_valueChanged: Calculate fitting parameters and display
void MainWindow::on_doubleSpinBox_powerLaw_f_valueChanged(double arg1)
{
    ui->qwtslider_powerLaw_f->setValue(arg1);
    _para->fRay = arg1;
    _para->bMie = ui->qwtslider_powerLaw_b->value();
    UpdateMuspFitPlot();
    UpdateMuspFitErrorDisplay();
}

//SpinBoxB_valueChanged: Calculate fitting parameters and display
void MainWindow::on_doubleSpinBox_powerLaw_b_valueChanged(double arg1)
{
    ui->qwtslider_powerLaw_b->setValue(arg1);
    _para->bMie = arg1;
    if (_para->fittingComplex)
        _para->fRay = ui->qwtslider_powerLaw_f->value();
    UpdateMuspFitPlot();
    UpdateMuspFitErrorDisplay();
}

//********************************** Update Functions **********************************
//Update Musp Fit
void MainWindow::UpdateMuspFitPlot()
{
    if (_otherPlotsFlag)
    {
        if (_para->nWavel >1)
        {
            PlotData plot;
            plot.AssignValuesMuspPowerLawPlots(ui,_para);
        }
    }
}

//Update S1S2 plot
void MainWindow::UpdateS1S2Plot()
{
    if (_otherPlotsFlag)
    {
        PlotData plot;
        plot.AssignValuesS1S2Plot(ui, _para);
    }
}

//Update Phase Function (Linear) plot
void MainWindow::UpdatePhaseFunctionLinearPlot()
{
    if (_otherPlotsFlag)
    {
        PlotData plot;
        plot.AssignValuesPhaseFunctionLinearPlot(ui, _para);
    }
}

//Update Phase Function (Polar) plot
void MainWindow::UpdatePhaseFunctionPolarPlot()
{
    if (_otherPlotsFlag)
    {
        PlotData plot;
        plot.AssignValuesPhaseFunctionPolarPlot(ui,_para);
    }
}

//Update Phase Function (Polar) plot
void MainWindow::UpdateMuspFitErrorDisplay()
{
    ui->label_currentMse->setText("<font color=\"red\">M.S. Error = </font>"+QString::number(_para->muspFittingError,'g',6));
    ui->label_currentA->setText("<font color=\"blue\">A = </font>"+
                                QString::number(_para->muspAtRefWavel,'g',6)+
                                " mm<sup>-1</sup>");
}


//********************************** When 'PolyDisperse:Custom' is selected **********************************

//comboBox_distribution_valueChanged, read data from a file
void MainWindow::on_comboBox_distribution_currentIndexChanged(int value)
{
    MainWindowSupport support;
    PlotData plot;
    plot.ClearPlots(ui,_para);
    _otherPlotsFlag = false;
    _distPlotFlag = false;
    // 0: Log Normal
    // 1: Gaussian
    // 2: Custom
    if (value == 2) //Custom Distribution
    {
        support.DisableWidgetsDuringCustomPolyDisperseData(ui, true);
        _loadCustomDataCheckFlag = true;
    }
    else
    {
        support.DisableWidgetsDuringCustomPolyDisperseData(ui, false);
        _loadCustomDataCheckFlag = false;
    }
}

//********************************** Mouse move over customplots **********************************

//mouseMove ScatteringCrossSectionPlot
void MainWindow::mouseMove_customPlot_csca(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_csca;
    QString strNameX = "WL";
    if (ui->radioButton_logYaxis->isChecked())
        strNameY = "Log(Csca)";
    else
        strNameY = "Csca";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}


//mouseMove ExtinctionCrossSectionPlot
void MainWindow::mouseMove_customPlot_cext(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_cext;
    QString strNameX = "WL";
    if (ui->radioButton_logYaxis->isChecked())
        strNameY = "Log(Cext)";
    else
        strNameY = "Cext";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//mouseMove BackscatteringCrossSectionPlot
void MainWindow::mouseMove_customPlot_cback(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_cback;
    QString strNameX = "WL";
    if (ui->radioButton_logYaxis->isChecked())
        strNameY = "Log(Cback)";
    else
        strNameY = "Cback";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//mouseMove SizeParameterPlot
void MainWindow::mouseMove_customPlot_sizePara(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_sizePara;
    QString strNameX = "WL";
    if (ui->radioButton_logYaxis->isChecked())
        strNameY = "Log(Size Parameter)";
    else
        strNameY = "Size Parameter";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//mouseMove MusPlot
void MainWindow::mouseMove_customPlot_mus(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_mus;
    QString strNameX = "WL";
    if (ui->radioButton_logYaxis->isChecked())
        strNameY = "Log(μs)";
    else
        strNameY = "μs";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//mouseMove MuspPlot
void MainWindow::mouseMove_customPlot_musp(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_musp;
    QString strNameX = "WL";
    if (ui->radioButton_logYaxis->isChecked())
        strNameY = "Log(μs')";
    else
        strNameY = "μs'";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//mouseMove GPlot
void MainWindow::mouseMove_customPlot_g(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_g;
    QString strNameX = "WL";
    if (ui->radioButton_logYaxis->isChecked())
        strNameY = "Log(g)";
    else
        strNameY = "g";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//mouseMove ForwardBackwardPlot
void MainWindow::mouseMove_customPlot_forwardBackward(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_forwardBackward;
    QString strNameX = "WL";
    if (ui->radioButton_logYaxis->isChecked())
        strNameY = "Log(%)";
    else
        strNameY = "%";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//mouseMove S1S2Plot
void MainWindow::mouseMove_customPlot_s1s2(QMouseEvent *event)
{
    QString strNameX = "Ang.";
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_s1s2;
    if (ui->radioButton_logYaxis->isChecked())
    {
        if (ui->radioButton_s1_amp->isChecked())
            strNameY = "Log(S1)";
        if (ui->radioButton_s2_amp->isChecked())
            strNameY = "Log(S2)";
    }
    else
    {
        if (ui->radioButton_s1_amp->isChecked())
            strNameY = "S1";
        if (ui->radioButton_s2_amp->isChecked())
            strNameY = "S2";
    }
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//mouseMove PhaseFunctionLinearPlot
void MainWindow::mouseMove_customPlot_pFunctionLinear(QMouseEvent *event)
{
    QString strNameX = "Ang.";
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_pFunctionLinear;
    if (ui->radioButton_logYaxis->isChecked())
        strNameY = "Log(Magnitude)";
    else
        strNameY = "Magnitude";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//mouseMove DistributionPlot
void MainWindow::mouseMove_customPlot_distribution(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_distribution;
    QString strNameX = "Dia.";
    if (ui->radioButton_logYaxis->isChecked())
        strNameY = "Log(Ns)";
    else
        strNameY = "Ns";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//on_mouseMove_ MuspPowerLawPlot
void MainWindow::mouseMove_customPlot_muspPowerLaw(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_muspPowerLaw;
    QString strNameX = "WL";
    if (ui->radioButton_logYaxis->isChecked())
        strNameY = "Log(μs')";
    else
        strNameY = "μs'";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//Display x and y values
void MainWindow::DisplayCurveData(QMouseEvent *event, QCustomPlot *curPlot,
                                  QString strNameX, QString strNameY)
{
    QCPAbstractPlottable *plottable = curPlot->plottableAt(event->localPos());
    if(plottable)
    {
        double x = curPlot->xAxis->pixelToCoord(event->localPos().x());
        double *dx;

        QCPGraph *graph =  qobject_cast<QCPGraph*>(plottable);
        if (graph)
        {
            double key = 0;
            double value = 0;
            bool ok = false;
            double maxx = std::numeric_limits<double>::max();

            QCPDataRange dataRange = graph->data()->dataRange();
            QCPGraphDataContainer::const_iterator begin = graph->data()->at(dataRange.begin()); // get range begin iterator from index
            QCPGraphDataContainer::const_iterator end = graph->data()->at(dataRange.end()); // get range end iterator from index

            unsigned int n = static_cast<unsigned int>(end-begin);
            if (n>0)
            {
                dx = new double[n];
                int index =0;
                for (QCPGraphDataContainer::const_iterator it=begin; it<end; it++)
                {
                    dx[index] = qAbs(x - it->key);
                    if ((dx[index] < maxx) )
                    {
                        key = it->key;
                        value = it->value;
                        ok = true;
                        maxx = dx[index];
                    }
                    index++;
                }
                delete[] dx;


                if (ok)
                {
                    QToolTip::showText(event->globalPos(),
                    tr("<table>"
                         "<tr>"
                           "<td>%L1:</td>" "<td>%L2</td>"
                         "</tr>"
                         "<tr>"
                           "<td>%L3:</td>" "<td>%L4</td>"
                         "</tr>"
                      "</table>").arg(strNameX).arg(key).arg(strNameY).arg(value),curPlot, curPlot->rect());
                }
            }
        }
    }
    else
        QToolTip::hideText();
}
