/**********************************************************************
** All MainWindow functions are listed in this file.
**********************************************************************/

#include "mainwindow.h"

/***********************Starting Main window and initializations ******************************/
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent), ui(new Ui::MainWindow)
{
    //Initialize UI
    ui->setupUi(this);

   //Intialize Window
    Initialize();

    //Connect signals and slots
    connect(ui->customPlot_Csca, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotScatteringCrossSection(QMouseEvent*)));
    connect(ui->customPlot_Cext, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotExtinctionCrossSection(QMouseEvent*)));
    connect(ui->customPlot_Cback, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotBackscatteringCrossSection(QMouseEvent*)));
    connect(ui->customPlot_SizePara, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotSizeParameter(QMouseEvent*)));
    connect(ui->customPlot_Mus, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotMus(QMouseEvent*)));
    connect(ui->customPlot_Musp, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotMusp(QMouseEvent*)));
    connect(ui->customPlot_G, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotG(QMouseEvent*)));
    connect(ui->customPlot_FB, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotForwardBackward(QMouseEvent*)));
    connect(ui->customPlot_S1S2, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotS1S2(QMouseEvent*)));
    connect(ui->customPlot_PhaseFunctionLinear, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotPhaseFunctionLinear(QMouseEvent*)));
    connect(ui->customPlot_Distribution, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotDistribution(QMouseEvent*)));
    connect(ui->customPlot_MuspPowerLaw, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotMuspPowerLaw(QMouseEvent*)));
}

MainWindow::~MainWindow()
{
    delete ui;
}

/***************************** Intialize *********************************************/
void MainWindow::Initialize()
{
    //Initialize paramaters
    mOtherPlotsFlag = false;
    mDistPlotFlag = false;
    mArrayFlag = false;
    mLoadCustomNoGoodFlag = false;
    mPara = new parameters();
    mPara->polarCurve = new QwtPolarCurve();

    MainWindowSupport support;
    support.InitializeGUI(ui);    
    support.SetWidgets(ui);

    PlotData plot;
    plot.InitializeDistributionPlot(ui);
    plot.InitializePhaseFunctionPolarPlot(ui,mPara);
    plot.InitializePhaseFunctionLinearPlot(ui);
    plot.InitializeAllOtherPlots(ui);
}

/********************************** Run Simulation ********************************************/
void MainWindow::on_pushButton_RunSimulation_clicked()
{    
    PlotData plot;
    MainWindowSupport support;

    //Initialize
    plot.ClearPlots(ui,mPara);
    support.SetWidgets(ui);
    support.LoadInputData(ui,mPara);
    support.SetWavelengthSliders(ui);              //set slider position

    if  (mLoadCustomNoGoodFlag)
    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Error");
        msgBox.setText("No Data available.");
        msgBox.setInformativeText("Click 'Load Custom Data'");
        msgBox.exec();
    }
    else
    {

        if (!support.CheckInputParameters(ui,mPara))    //sanity check
        {
            //Disable widgets
            support.DisableWidgetsDuringSimulation(ui, mPara, true);
            //If mArrayFlag is TRUE, delete dynamic arrays before next simulation
            if (mArrayFlag)
                support.DeleteArrays(mPara, &mArrayFlag);
            //Initialize dynamic arrays
            support.InitializeArrays(ui,mPara, &mArrayFlag);

            //Mono disperse
            if (ui->radioButton_MonoDisperse->isChecked())
            {
                support.ProcessDistribution(ui, mPara, 3);
                support.ProcessMonoDisperse(ui,mPara);
                ui->label_Progress->setText("<font color=\"green\"> Completed!</font>");
                ui->slider_ConcPercentChange->setValue(0);
                mOtherPlotsFlag = true;
                mDistPlotFlag = true;
            }

            //Poly disperse
            if (ui->radioButton_PolyDisperse->isChecked())
            {
                if (!support.CheckDistribution(ui,mPara))   //sanity check
                {
                    support.ProcessDistribution(ui, mPara, static_cast<unsigned int>(ui->comboBox_Distribution->currentIndex()));
                    support.ProcessPolyDisperse(ui,mPara);
                    ui->label_Progress->setText("<font color=\"green\"> Completed!</font>");
                    ui->slider_ConcPercentChange->setValue(0);
                    mOtherPlotsFlag = true;
                    mDistPlotFlag = true;
                }
                if (ui->comboBox_Distribution->currentIndex() ==2)
                    support.DisableWidgetsDuringCustomPolyDisperseData(ui, true);
                else
                    support.DisableWidgetsDuringCustomPolyDisperseData(ui, false);
            }
            //Enable widgets
            support.DisableWidgetsDuringSimulation(ui, mPara, false);
        }        
    }
}

/********************************** Show distribution ********************************************/
void MainWindow::on_pushButton_ShowDistributionAndCustom_clicked()
{
    MainWindowSupport support;
    bool dataValidFlag = false;

    if (ui->comboBox_Distribution->currentIndex() == 2)
    {
        MainWindowSupport support;
        QString fileName = QFileDialog::getOpenFileName
          (this,tr("Read Data"),"",tr("Text File (*.txt);;All Files (*)"));

        if (fileName.isEmpty())
            return;
        else
            support.ReadCustomData(mPara, fileName, &dataValidFlag);
    }

    if (dataValidFlag)
    {
        PlotData plot;
        plot.ClearPlots(ui,mPara);
        mOtherPlotsFlag = false;
        mDistPlotFlag = false;
        mLoadCustomNoGoodFlag = true;
    }
    else
    {
        support.LoadInputData(ui, mPara);
        if (!support.CheckInputParameters(ui,mPara))    //sanity check
        {
            if (!support.CheckDistribution(ui,mPara))   //sanity check
            {
                support.SetWidgets(ui);
                support.ProcessDistribution(ui, mPara, static_cast<unsigned int>(ui->comboBox_Distribution->currentIndex()));
                mDistPlotFlag = true;
            }
        }
        mLoadCustomNoGoodFlag = false;
    }
}

/**********************************  Save Data ********************************************/
void MainWindow::on_pushButton_SaveData_clicked()
{
    if (mOtherPlotsFlag)
    {
        OptionsDialog saveOptions;
        saveOptions.SaveData(ui, mPara);
    }
    else
    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Error");
        msgBox.setText("No Data available.");
        msgBox.exec();
    }
}

/********************************** Display Data ********************************************/
void MainWindow::on_pushButton_DisplayData_clicked()
{
    if (mOtherPlotsFlag)
    {
        DisplayDialog *display;
        display = new DisplayDialog();
        display->DisplayData(ui,  mPara);
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

/********************************** Calculate Best Fit musp ********************************************/
void MainWindow::on_pushButton_BestFit_clicked()
{    
    if(mOtherPlotsFlag)
    {
        if (mPara->nWavel >1)
        {
            MainWindowSupport support;
            utilities util;
            calculate cal;
            PlotData plot;
            double fRay, bMie;

            //Widget settings
            support.DisableWidgetsDuringSimulation(ui, mPara, true);
            ui->label_CurrentMSE->setText("<font color=\"red\">Wait!...    </font>");
            ui->label_CurrentA->setText("");
            util.Delay();

            //Calculate fRay and bMie
            if (mPara->fittingComplex)
                cal.CalculatePowerLawAutoFitComplex(mPara);
            else
                cal.CalculatePowerLawAutoFitSimple(mPara);
            plot.AssignValuesMuspPowerLawPlots(ui,mPara);
            fRay = mPara->fRay;
            bMie = mPara->bMie;
            ui->doubleSpinBox_F->setValue(fRay);
            ui->doubleSpinBox_B->setValue(bMie);
            UpdateMuspFitErrorDisplay();

            //Enable Widgets
            support.DisableWidgetsDuringSimulation(ui, mPara, false);
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

/********************************** Close GUI ********************************************/
void MainWindow::on_pushButton_Close_clicked()
{
   this->close();
}

/********************************** Radio Button Selection **********************************/

//radioButton_LinearYAxis_clicked: Select Linear Y-Axis
void MainWindow::on_radioButton_LinearYAxis_clicked()
{
    PlotData plot;
    if (mOtherPlotsFlag)
    {
        MainWindowSupport support;
        support.DisableEnableRealImagButtons(ui);
        plot.AssignValuesS1S2Plot(ui, mPara);
        plot.AssignValuesPhaseFunctionLinearPlot(ui, mPara);
        plot.AssignValuesAllOtherPlots(ui, mPara);
    }
    if (mDistPlotFlag)
        plot.AssignValuesDistributionPlot(ui, mPara);
}

//radioButton_LogYAxis_clicked: Select Log Y-Axis
void MainWindow::on_radioButton_LogYAxis_clicked()
{
    PlotData plot;
    if (mOtherPlotsFlag)
    {
        MainWindowSupport support;
        support.DisableEnableRealImagButtons(ui);        
        plot.AssignValuesS1S2Plot(ui, mPara);
        plot.AssignValuesPhaseFunctionLinearPlot(ui, mPara);
        plot.AssignValuesAllOtherPlots(ui, mPara);
    }
    if (mDistPlotFlag)
        plot.AssignValuesDistributionPlot(ui, mPara);
}

//radioButton_LinearXAxis_clicked: Select Linear X-Axis
void MainWindow::on_radioButton_LinearXAxis_clicked()
{    
    if (mDistPlotFlag)
    {
        PlotData plot;
        plot.AssignValuesDistributionPlot(ui, mPara);
    }
}

//radioButton_LogXAxis_clicked: Select Log X-Axis
void MainWindow::on_radioButton_LogXAxis_clicked()
{
    if (mDistPlotFlag)
    {
        PlotData plot;
        plot.AssignValuesDistributionPlot(ui, mPara);
    }
}

//radioButtonMonoDisperse_clicked: Select Mono disperse distribution
void MainWindow::on_radioButton_MonoDisperse_clicked()
{    
    MainWindowSupport support;
    ui->comboBox_Distribution->setCurrentIndex(0);
    support.SetWidgets(ui);
    PlotData plot;
    plot.ClearPlots(ui,mPara);
    ui->tabWidget_PhaseFunction->setTabText(0,"S1/S2");
    ui->tabWidget_PhaseFunction->setTabText(1,"Phase Function (Polar)");
    ui->tabWidget_PhaseFunction->setTabText(2,"Phase Function (Linear)");
    ui->tabWidget_Scat->setTabText(1,"Csca");
    ui->tabWidget_Scat->setTabText(2,"Cext");
    ui->tabWidget_Scat->setTabText(3,"Cback");           
    ui->tabWidget_Ns->setTabText(0,"Concentration (Number Density: Ns)");
    ui->tabWidget_Ns->setTabEnabled(1,true);
    ui->tabWidget_Ns->setTabText(1,"Size Parameter");
    mOtherPlotsFlag = false;
    mDistPlotFlag = false;
}

//radioButtonPolyDisperse_clicked: Select Poly disperse distribution
void MainWindow::on_radioButton_PolyDisperse_clicked()
{
    MainWindowSupport support;
    support.SetWidgets(ui);    
    PlotData plot;
    plot.ClearPlots(ui,mPara);    
    ui->tabWidget_PhaseFunction->setTabText(0,"Ave. S1/S2");
    ui->tabWidget_PhaseFunction->setTabText(1,"Ave. Phase Function (Polar)");
    ui->tabWidget_PhaseFunction->setTabText(2,"Ave. Phase Function (Linear)");
    ui->tabWidget_Scat->setTabText(1,"Ave. Csca");
    ui->tabWidget_Scat->setTabText(2,"Ave. Cext");
    ui->tabWidget_Scat->setTabText(3,"Ave. Cback");
        ui->tabWidget_Ns->setTabText(0,"Scatterer Distribution and Concentration (Ns)");
    ui->tabWidget_Ns->setTabEnabled(1,false);    
    ui->tabWidget_Ns->setTabText(1,"");
    mOtherPlotsFlag = false;
    mDistPlotFlag = false;
}

//radioButton_VolFrac_clicked: Select Volume Fraction (sphere volume/cubic mm)
void MainWindow::on_radioButton_VolFrac_clicked()
{
    mPara->volFraction = ui->lineEdit_VolFrac->text().toDouble();
    ui->lineEdit_Conc_mm3->setText(QString::number(mPara->sphNumDensity));

    MainWindowSupport support;
    support.SetWidgets(ui);
    PlotData plot;
    plot.ClearPlots(ui,mPara);
    mOtherPlotsFlag = false;
    mDistPlotFlag = false;
}

//radioButton_Conc_mm3_clicked: Select Concentraction (number of spheres in 1 cubic mm)
void MainWindow::on_radioButton_Conc_mm3_clicked()
{
    mPara->sphNumDensity = ui->lineEdit_Conc_mm3->text().toDouble();
    ui->lineEdit_VolFrac->setText(QString::number(mPara->volFraction));

    MainWindowSupport support;
    support.SetWidgets(ui);
    PlotData plot;
    plot.ClearPlots(ui,mPara);
    mOtherPlotsFlag = false;
    mDistPlotFlag = false;
}

//radioButtonS1_clicked: Select S1 plot
void MainWindow::on_radioButton_S1_clicked()
{
    UpdateS1S2Plot();
}

//radioButtonS2_clicked: Select S2 plot
void MainWindow::on_radioButton_S2_clicked()
{
    UpdateS1S2Plot();
}

//radioButton_S1S2_Abs_clicked: Select absolute value
void MainWindow::on_radioButton_S1S2_Abs_clicked()
{
    UpdateS1S2Plot();
}

//radioButton_S1S2_Real_clicked: Select Real component
void MainWindow::on_radioButton_S1S2_Real_clicked()
{
    UpdateS1S2Plot();
}

//radioButton_S1S2_Imag_clicked: Select Imaginary component
void MainWindow::on_radioButton_S1S2_Imag_clicked()
{
    UpdateS1S2Plot();
}

//radioButton_PhaseLinear_clicked: Select Linear Scale
void MainWindow::on_radioButton_PhaseLinear_clicked()
{
    UpdatePhaseFunctionPolarPlot();
}

//radioButton_PhaseLog_clicked: Select Log Scale
void MainWindow::on_radioButton_PhaseLog_clicked()
{
    UpdatePhaseFunctionPolarPlot();
}

//radioButton_PhaseAverage_clicked: Select average
void MainWindow::on_radioButton_PhaseAverage_clicked()
{
    UpdatePhaseFunctionPolarPlot();
}

//radioButton_PhasePara_clicked: Select Parallel component
void MainWindow::on_radioButton_PhasePara_clicked()
{
    UpdatePhaseFunctionPolarPlot();
}

//radioButton_PhasePerp_clicked: Select Perpendicular component
void MainWindow::on_radioButton_PhasePerp_clicked()
{
    UpdatePhaseFunctionPolarPlot();
}

//radioButton_dTheta0_1 is clicked: Simulation time warning
void MainWindow::on_radioButton_Phase_DTheta0_1_clicked()
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
          ui->radioButton_Phase_DTheta0_1->setChecked(true);
          break;
      case QMessageBox::No:
          ui->radioButton_Phase_DTheta0_5->setChecked(true);
          break;
    }
}

/********************************** Slider Selection **********************************/

//sliderConcPercentChange_slidervalueChanged: Change Conc or VolFrac
void MainWindow::on_slider_ConcPercentChange_valueChanged(int position)
{
    double value;
    double margin = (1.0 + position /200.0);
    ui->label_ActualConcPercent->setText(QString::number(position/2.0)+"%");

    if (ui->radioButton_Conc_mm3->isChecked())
    {
        value =  mPara->sphNumDensity * margin;
        ui->lineEdit_Conc_mm3->setText(QString::number(value));
    }
    if (ui->radioButton_VolFrac->isChecked())
    {
        value =  mPara->volFraction * margin;
        ui->lineEdit_VolFrac->setText(QString::number(value));
    }

    PlotData plot;
    if (mDistPlotFlag)
        plot.AssignValuesDistributionPlot(ui, mPara);
    if (mOtherPlotsFlag)
    {   MainWindowSupport support;
        support.DisableEnableRealImagButtons(ui);
        plot.AssignValuesAllOtherPlots(ui, mPara);
    }
}

//on_slider_WL_PFPolar_valueChanged: Change Phase Function Wavelenghth Slider
void MainWindow::on_slider_WL_PFPolar_valueChanged(int value)
{
   int wavel = static_cast<int>(mPara->startWavel + value*mPara->stepWavel);
   ui->label_CurrentWL_PFPolar->setText(QString::number(wavel));
   UpdatePhaseFunctionPolarPlot();
}

//on_slider_WL_PFLinear_valueChanged: Change Phase Function Wavelenghth Slider
void MainWindow::on_slider_WL_PFLinear_valueChanged(int value)
{
   int wavel = static_cast<int>(mPara->startWavel + value*mPara->stepWavel);
   ui->label_CurrentWL_PFLinear->setText(QString::number(wavel));
   UpdatePhaseFunctionLinearPlot();
}

//on_slider_WL_S1S2_valueChanged: Change S1/S2 Wavelenghth Slider
void MainWindow::on_slider_WL_S1S2_valueChanged(int value)
{
   int wavel = static_cast<int>(mPara->startWavel + value*mPara->stepWavel);
   ui->label_CurrentWL_S1S2->setText(QString::number(wavel));
   UpdateS1S2Plot();
}

//SpinBoxF_valueChanged: Calculate fitting parameters and display
void MainWindow::on_doubleSpinBox_F_valueChanged(double arg1)
{
    ui->qwtslider_F->setValue(arg1);
    mPara->fRay = arg1;
    mPara->bMie = ui->qwtslider_B->value();
    UpdateMuspFitPlot();
    UpdateMuspFitErrorDisplay();
}

//SpinBoxB_valueChanged: Calculate fitting parameters and display
void MainWindow::on_doubleSpinBox_B_valueChanged(double arg1)
{
    ui->qwtslider_B->setValue(arg1);
    mPara->bMie = arg1;
    if (mPara->fittingComplex)
        mPara->fRay = ui->qwtslider_F->value();
    UpdateMuspFitPlot();
    UpdateMuspFitErrorDisplay();
}

//radioButton_FittingSimple_clicked: Select Fitting Simple
void MainWindow::on_radioButton_FittingSimple_clicked()
{
    mPara->fittingComplex = false;
    ui->qwtslider_F->setDisabled(true);
    ui->qwtslider_F->setValue(0);
    ui->doubleSpinBox_F->setDisabled(true);
    QString labelB = "<P><b><FONT COLOR='#aa0000' FONT SIZE = 4>";
    labelB .append("b");
    labelB .append("</b></P></br>");
    ui->label_BLabel->setText(labelB);
    ui->label_FLabel->setText("");
}

//radioButton_FittingComplex_clicked: Select Fitting Complex
void MainWindow::on_radioButton_FittingComplex_clicked()
{
    mPara->fittingComplex = true;
    ui->qwtslider_F->setDisabled(false);
    ui->doubleSpinBox_F->setDisabled(false);
    QString labelB = "<P><b><FONT COLOR='#aa0000' FONT SIZE = 4>";
    labelB .append("b<sub>Mie</sub>");
    labelB .append("</b></P></br>");
    ui->label_BLabel->setText(labelB);
    QString labelF = "<P><b><FONT COLOR='#005500' FONT SIZE = 4>";
    labelF .append("f<sub>Ray</sub>");
    labelF .append("</b></P></br>");
    ui->label_FLabel->setText(labelF);
}

//radioButton_RefWavel500_clicked: Set RefWavel = 500
void MainWindow::on_radioButton_RefWavel500_clicked()
{
    mPara->refWavel = 500.0;
    mPara->refWavelIdx = 0;
}

//radioButton_RefWavel500_clicked: Set RefWavel = 500
void MainWindow::on_radioButton_RefWavel600_clicked()
{
    mPara->refWavel = 600.0;
    mPara->refWavelIdx = 1;
}

//radioButton_RefWavel500_clicked: Set RefWavel = 500
void MainWindow::on_radioButton_RefWavel700_clicked()
{
    mPara->refWavel = 700.0;
    mPara->refWavelIdx = 2;
}

//radioButton_RefWavel500_clicked: Set RefWavel = 500
void MainWindow::on_radioButton_RefWavel800_clicked()
{
    mPara->refWavel = 800.0;
    mPara->refWavelIdx = 3;
}

//radioButton_RefWavel500_clicked: Set RefWavel = 500
void MainWindow::on_radioButton_RefWavel900_clicked()
{
    mPara->refWavel = 900.0;
    mPara->refWavelIdx = 4;
}

//radioButton_RefWavel1000_clicked: Set RefWavel = 1000
void MainWindow::on_radioButton_RefWavel1000_clicked()
{
    mPara->refWavel = 1000.0;
    mPara->refWavelIdx = 5;
}

/********************************** Update Functions **********************************/
//Update Musp Fit
void MainWindow::UpdateMuspFitPlot()
{
    if (mOtherPlotsFlag)
    {
        if (mPara->nWavel >1)
        {
            PlotData plot;
            plot.AssignValuesMuspPowerLawPlots(ui,mPara);
        }
    }
}

//Update S1S2 plot
void MainWindow::UpdateS1S2Plot()
{
    if (mOtherPlotsFlag)
    {
        PlotData plot;
        plot.AssignValuesS1S2Plot(ui, mPara);
    }
}

//Update Phase Function (Linear) plot
void MainWindow::UpdatePhaseFunctionLinearPlot()
{
    if (mOtherPlotsFlag)
    {
        PlotData plot;
        plot.AssignValuesPhaseFunctionLinearPlot(ui, mPara);
    }
}

//Update Phase Function (Polar) plot
void MainWindow::UpdatePhaseFunctionPolarPlot()
{
    if (mOtherPlotsFlag)
    {
        PlotData plot;
        plot.AssignValuesPhaseFunctionPolarPlot(ui,mPara);
    }
}

//Update Phase Function (Polar) plot
void MainWindow::UpdateMuspFitErrorDisplay()
{
    ui->label_CurrentMSE->setText("<font color=\"red\">M.S. Error = </font>"+QString::number(mPara->muspFittingError,'g',6));
    ui->label_CurrentA->setText("<font color=\"blue\">a = </font>"+
                                QString::number(mPara->muspAtRefWavel[mPara->refWavelIdx],'g',6)+
                                " mm<sup>-1</sup>");
}


/********************************** When 'PolyDisperse:Custom' is selected **********************************/

//comboBox_Distribution_valueChanged, read data from a file
void MainWindow::on_comboBox_Distribution_currentIndexChanged(int value)
{
    MainWindowSupport support;
    PlotData plot;
    plot.ClearPlots(ui,mPara);
    mOtherPlotsFlag = false;
    mDistPlotFlag = false;
    // 0: Log Normal
    // 1: Gaussian
    // 2: Custom
    if (value == 2) //Custom Distribution
    {
        support.DisableWidgetsDuringCustomPolyDisperseData(ui, true);
        mLoadCustomNoGoodFlag = true;
    }
    else
    {
        support.DisableWidgetsDuringCustomPolyDisperseData(ui, false);
        mLoadCustomNoGoodFlag = false;
    }
}

/********************************** Mouse over customplots **********************************/

//MouseOver PlotScatteringCrossSection
void MainWindow::MouseOverPlotScatteringCrossSection(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_Csca;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(Csca)";
    else
        strNameY = "Csca";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}


//MouseOver PlotExtinctionCrossSection
void MainWindow::MouseOverPlotExtinctionCrossSection(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_Cext;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(Cext)";
    else
        strNameY = "Cext";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//MouseOver PlotBackscatteringCrossSection
void MainWindow::MouseOverPlotBackscatteringCrossSection(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_Cback;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(Cback)";
    else
        strNameY = "Cback";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//MouseOver PlotSizeParameter
void MainWindow::MouseOverPlotSizeParameter(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_SizePara;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(Size Parameter)";
    else
        strNameY = "Size Parameter";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//MouseOver PlotMus
void MainWindow::MouseOverPlotMus(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_Mus;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(μs)";
    else
        strNameY = "μs";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//MouseOver PlotMusp
void MainWindow::MouseOverPlotMusp(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_Musp;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(μs')";
    else
        strNameY = "μs'";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//MouseOver PlotG
void MainWindow::MouseOverPlotG(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_G;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(g)";
    else
        strNameY = "g";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//MouseOver PlotForwardBackward
void MainWindow::MouseOverPlotForwardBackward(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_FB;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(%)";
    else
        strNameY = "%";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//MouseOver PlotS1S2
void MainWindow::MouseOverPlotS1S2(QMouseEvent *event)
{    
    QString strNameX = "Ang.";
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_S1S2;
    if (ui->radioButton_LogYAxis->isChecked())
    {
        if (ui->radioButton_S1->isChecked())
            strNameY = "Log(S1)";
        if (ui->radioButton_S2->isChecked())
            strNameY = "Log(S2)";
    }
    else
    {
        if (ui->radioButton_S1->isChecked())
            strNameY = "S1";
        if (ui->radioButton_S2->isChecked())
            strNameY = "S2";
    }
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//MouseOver PlotPhaseFunctionLinear
void MainWindow::MouseOverPlotPhaseFunctionLinear(QMouseEvent *event)
{
    QString strNameX = "Ang.";
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_PhaseFunctionLinear;
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(Magnitude)";
    else
        strNameY = "Magnitude";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//MouseOver PlotDistribution
void MainWindow::MouseOverPlotDistribution(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_Distribution;
    QString strNameX = "Dia.";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(Ns)";
    else
        strNameY = "Ns";
    DisplayCurveData(event, curPlot, strNameX, strNameY);
}

//MouseOver Musp power Law
void MainWindow::MouseOverPlotMuspPowerLaw(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_MuspPowerLaw;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
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
