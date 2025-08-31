/**********************************************************************
** All MainWindow functions are listed in this file.
**********************************************************************/

#include "mainwindow.h"
#include "calc/calculate.h"
#include "calc/utilities.h"
#include "dialog/plotdata.h"
#include "dialog/optionsdialog.h"
#include "dialog/displaydialog.h"
#include "dialog/mainwindowsupport.h"


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
    connect(ui->customPlot_PhaseFunctionPolar, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotPhaseFunctionPolar(QMouseEvent*)));
    connect(ui->customPlot_Distribution, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotDistribution(QMouseEvent*)));
    connect(ui->customPlot_MuspPowerLaw, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotMuspPowerLaw(QMouseEvent*)));
    connect(ui->slider_F, &QSlider::valueChanged, this, &::MainWindow::Slider_F_valueChanged);
    connect(ui->slider_B, &QSlider::valueChanged, this, &::MainWindow::Slider_B_valueChanged);
    connect(ui->doubleSpinBox_F, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &::MainWindow::DoubleSpinBox_F_valueChanged);
    connect(ui->doubleSpinBox_B, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &::MainWindow::DoubleSpinBox_B_valueChanged);
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
    mPara = new Parameters();

    MainWindowSupport support;
    support.InitializeGUI(ui, mPara);
    support.SetWidgets(ui, mPara);

    PlotData plot;
    plot.InitialSetupDistributionPlot(ui);
    plot.InitialSetupPhaseFunctionPolarPlot(ui);
    plot.InitialSetupPhaseFunctionLinearPlot(ui);
    plot.InitialSetupS1S2Plot(ui);
    plot.InitialSetupMuspPowerLawFit(ui);
    plot.InitialSetupOtherPlots(ui);
}

/********************************** Update Functions **********************************/
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
        plot.SetupPolarPlotForData(ui, mPara);
        plot.AssignValuesPhaseFunctionPolarPlot(ui,mPara);
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

//Update Musp Fit plot
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

//Update Musp error display
void MainWindow::UpdateMuspFitErrorDisplay()
{
    ui->label_CurrentMSE->setText("<font color=\"red\">M.S. Error = "
                                  "</font>"+QString::number(mPara->muspFittingError,'g',6));
    ui->label_CurrentA->setText("A = "+
                                QString::number(mPara->muspAtRefWavel[mPara->refWavelIdx],'g',6)+
                                " mm<sup>-1</sup>");
}

/********************************** Run Simulation ********************************************/
void MainWindow::on_pushButton_RunSimulation_clicked()
{    
    PlotData plot;
    MainWindowSupport support;

    //Initialize
    plot.ClearPlots(ui);
    support.SetWidgets(ui, mPara);
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
        if (mPara->CheckCommonParameters(ui->radioButton_MonoDisperse,
                                          ui->radioButton_NumDen,
                                          ui->radioButton_VolFrac))
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
                support.ProcessDistribution(ui, mPara, mPara->MonoDisperse);  // index = 3: mono disperse
                support.ProcessMonoDisperse(ui,mPara);
                ui->label_Progress->setText("<font color=\"green\"> Completed!</font>");
                ui->slider_ConcPercentChange->setValue(0);
                mOtherPlotsFlag = true;
                mDistPlotFlag = true;
            }

            //Poly disperse
            if (ui->radioButton_PolyDisperse->isChecked())
            {
                if (mPara->CheckDistributionParameters(ui->comboBox_Distribution))   //sanity check
                {
                    support.ProcessDistribution(ui, mPara, static_cast<unsigned int>(ui->comboBox_Distribution->currentIndex()));
                    support.ProcessPolyDisperse(ui,mPara);
                    ui->label_Progress->setText("<font color=\"green\"> Completed!</font>");
                    ui->slider_ConcPercentChange->setValue(0);
                    mOtherPlotsFlag = true;
                    mDistPlotFlag = true;
                }
                if (ui->comboBox_Distribution->currentIndex() == mPara->Custom)
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

    if (ui->comboBox_Distribution->currentIndex() == mPara->Custom)
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
        plot.ClearPlots(ui);
        mOtherPlotsFlag = false;
        mDistPlotFlag = false;
        mLoadCustomNoGoodFlag = true;
    }
    else
    {
        support.LoadInputData(ui, mPara);
        if (mPara->CheckCommonParameters(ui->radioButton_MonoDisperse,
                                          ui->radioButton_NumDen,
                                          ui->radioButton_VolFrac))    //sanity check
        {
            if (mPara->CheckDistributionParameters(ui->comboBox_Distribution))   //sanity check
            {
                support.SetWidgets(ui, mPara);
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
        saveOptions.SaveData(ui->radioButton_MonoDisperse,
                             ui->radioButton_PolyDisperse,
                             ui->radioButton_NumDen,
                             ui->radioButton_VolFrac,
                             ui->comboBox_Distribution,
                             ui->slider_ConcPercentChange,
                             ui->radioButton_PhaseAverage,
                             ui->radioButton_PhasePara,
                             ui->radioButton_PhasePerp,mPara);
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
        display->DisplayData(ui->radioButton_MonoDisperse,
                             ui->radioButton_PolyDisperse,
                             ui->radioButton_NumDen,
                             ui->radioButton_VolFrac,
                             ui->comboBox_Distribution,
                             ui->slider_ConcPercentChange,
                             ui->slider_WL_PFPolar,
                             ui->radioButton_PhaseAverage,
                             ui->radioButton_PhasePara,
                             ui->radioButton_PhasePerp,  mPara);
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
            Utilities util;
            Calculate cal;
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

/********************************** Save Plots ********************************************/

//SavePlot_Mus
void MainWindow::on_pushButton_SavePlot_Mus_clicked()
{
    SavePlot(ui->customPlot_Mus, "Mus_plot");
}

//SavePlot_Csca
void MainWindow::on_pushButton_SavePlot_Csca_clicked()
{
    SavePlot(ui->customPlot_Csca, "Csca_plot");
}

//SavePlot_Cext
void MainWindow::on_pushButton_SavePlot_Cext_clicked()
{
    SavePlot(ui->customPlot_Cext, "Cext_plot");
}

//SavePlot_Cback
void MainWindow::on_pushButton_SavePlot_Cback_clicked()
{
    SavePlot(ui->customPlot_Cback, "Cback_plot");
}

//SavePlot_Musp
void MainWindow::on_pushButton_SavePlot_Musp_clicked()
{
    SavePlot(ui->customPlot_Musp, "Musp_plot");
}

//SavePlot_MuspPowerLaw
void MainWindow::on_pushButton_SavePlot_MuspPowerLaw_clicked()
{
    SavePlot(ui->customPlot_MuspPowerLaw, "MuspPowerLaw_plot");
}

//SavePlot_Distribution
void MainWindow::on_pushButton_SavePlot_Distribution_clicked()
{
    SavePlot(ui->customPlot_Distribution, "Distribution_plot");
}

//SavePlot_SizePara
void MainWindow::on_pushButton_SavePlot_SizePara_clicked()
{
    SavePlot(ui->customPlot_SizePara, "SizePara_plot");
}

//SavePlot_S1S2
void MainWindow::on_pushButton_SavePlot_S1S2_clicked()
{
    SavePlot(ui->customPlot_S1S2, "S1S2_plot");
}

//SavePlot_PhaseFunctionPolar
void MainWindow::on_pushButton_SavePlot_PhaseFunctionPolar_clicked()
{
    SavePlot(ui->customPlot_PhaseFunctionPolar, "PhaseFunctionPolar_plot");
}

//SavePlot_PhaseFunctionLinear
void MainWindow::on_pushButton_SavePlot_PhaseFunctionLinear_clicked()
{
    SavePlot(ui->customPlot_PhaseFunctionLinear, "PhaseFunctionLinear_plot");
}

//SavePlot_G
void MainWindow::on_pushButton_SavePlot_G_clicked()
{
    SavePlot(ui->customPlot_G, "G_plot");
}

//SavePlot_ForwardBackwardScat
void MainWindow::on_pushButton_SavePlot_FB_clicked()
{
    SavePlot(ui->customPlot_FB, "ForwardBackwardScat_plot");
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
        plot.AssignValuesOtherPlots(ui, mPara);
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
        plot.AssignValuesOtherPlots(ui, mPara);
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
    ui->comboBox_Distribution->setCurrentIndex(mPara->LogNormal);
    support.SetWidgets(ui, mPara);
    PlotData plot;
    plot.ClearPlots(ui);
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
    support.SetWidgets(ui, mPara);
    PlotData plot;
    plot.ClearPlots(ui);
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

//radioButton_NumDen_clicked: Select Concentraction (number of spheres in 1 cubic mm)
void MainWindow::on_radioButton_NumDen_clicked()
{
    mPara->sphNumDensity = ui->lineEdit_NumDen->text().toDouble();
    ui->lineEdit_VolFrac->setText(QString::number(mPara->volFraction));

    MainWindowSupport support;
    support.SetWidgets(ui, mPara);
    PlotData plot;
    plot.ClearPlots(ui);
    mOtherPlotsFlag = false;
    mDistPlotFlag = false;
}

//radioButton_VolFrac_clicked: Select Volume Fraction (sphere volume/cubic mm)
void MainWindow::on_radioButton_VolFrac_clicked()
{
    mPara->volFraction = ui->lineEdit_VolFrac->text().toDouble();
    ui->lineEdit_NumDen->setText(QString::number(mPara->sphNumDensity));

    MainWindowSupport support;
    support.SetWidgets(ui, mPara);
    PlotData plot;
    plot.ClearPlots(ui);
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

//radioButtonS2_clicked: Select S1/S2 plot
void MainWindow::on_radioButton_S1S2_clicked()
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

//radioButton_PhasePerp_clicked: Select both Parallel and Perpendicular components
void MainWindow::on_radioButton_PhaseAll_clicked()
{
    UpdatePhaseFunctionPolarPlot();
}

//radioButton_dTheta0_1 is clicked: Simulation time warning
void MainWindow::on_radioButton_Phase_DTheta0_1_clicked()
{
    QMessageBox msgBox;

    msgBox.setIcon(QMessageBox::Warning);
    msgBox.setMinimumSize(320,180);
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

//radioButton_FittingSimple_clicked: Select Fitting Simple
void MainWindow::on_radioButton_FittingSimple_clicked()
{
    mPara->fittingComplex = false;
    ui->slider_F->setDisabled(true);
    ui->slider_F->setValue(0);
    ui->labelF_0_0->setDisabled(true);
    ui->labelF_0_2->setDisabled(true);
    ui->labelF_0_4->setDisabled(true);
    ui->labelF_0_6->setDisabled(true);
    ui->labelF_0_8->setDisabled(true);
    ui->labelF_1_0->setDisabled(true);
    ui->doubleSpinBox_F->setDisabled(true);
    QString labelB = "<P><b><FONT COLOR='#ff5500'>";
    labelB .append("b");
    labelB .append("</b></P></br>");
    ui->label_BLabel->setText(labelB);
    ui->label_FLabel->setText("");
}

//radioButton_FittingComplex_clicked: Select Fitting Complex
void MainWindow::on_radioButton_FittingComplex_clicked()
{
    mPara->fittingComplex = true;
    ui->slider_F->setDisabled(false);
    ui->doubleSpinBox_F->setDisabled(false);
    ui->labelF_0_0->setDisabled(false);
    ui->labelF_0_2->setDisabled(false);
    ui->labelF_0_4->setDisabled(false);
    ui->labelF_0_6->setDisabled(false);
    ui->labelF_0_8->setDisabled(false);
    ui->labelF_1_0->setDisabled(false);
    QString labelB = "<P><b><FONT COLOR='#ff5500'>";
    labelB .append("b<sub>Mie</sub>");
    labelB .append("</b></P></br>");
    ui->label_BLabel->setText(labelB);
    QString labelF = "<P><b><FONT COLOR='#00c800'>";
    labelF .append("f<sub>Ray</sub>");
    labelF .append("</b></P></br>");
    ui->label_FLabel->setText(labelF);
}

//radioButton_RefWavel500_clicked: Set RefWavel = 500
void MainWindow::on_radioButton_RefWavel500_clicked()
{
    mPara->refWavel = 500.0;
    mPara->refWavelIdx = mPara->wavel500;
}

//radioButton_RefWavel500_clicked: Set RefWavel = 600
void MainWindow::on_radioButton_RefWavel600_clicked()
{
    mPara->refWavel = 600.0;
    mPara->refWavelIdx = mPara->wavel600;
}

//radioButton_RefWavel500_clicked: Set RefWavel = 700
void MainWindow::on_radioButton_RefWavel700_clicked()
{
    mPara->refWavel = 700.0;
    mPara->refWavelIdx = mPara->wavel700;
}

//radioButton_RefWavel500_clicked: Set RefWavel = 800
void MainWindow::on_radioButton_RefWavel800_clicked()
{
    mPara->refWavel = 800.0;
    mPara->refWavelIdx = mPara->wavel800;
}

//radioButton_RefWavel500_clicked: Set RefWavel = 900
void MainWindow::on_radioButton_RefWavel900_clicked()
{
    mPara->refWavel = 900.0;
    mPara->refWavelIdx = mPara->wavel900;
}

//radioButton_RefWavel1000_clicked: Set RefWavel = 1000
void MainWindow::on_radioButton_RefWavel1000_clicked()
{
    mPara->refWavel = 1000.0;
    mPara->refWavelIdx = mPara->wavel1000;
}

/********************************** Slider Selection **********************************/

//sliderConcPercentChange_slidervalueChanged: Change Conc or VolFrac
void MainWindow::on_slider_ConcPercentChange_valueChanged(int position)
{
    double value;
    double margin = (1.0 + position /200.0);
    ui->label_ActualConcPercent->setText(QString::number(position/2.0)+"%");

    if (ui->radioButton_NumDen->isChecked())
    {
        value =  mPara->sphNumDensity * margin;
        ui->lineEdit_NumDen->setText(QString::number(value));
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
        plot.AssignValuesOtherPlots(ui, mPara);
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

/********************************** When 'PolyDisperse:Custom' is selected ***********************/

//comboBox_Distribution_valueChanged, read data from a file
void MainWindow::on_comboBox_Distribution_currentIndexChanged(int value)
{
    MainWindowSupport support;
    PlotData plot;
    plot.ClearPlots(ui);
    mOtherPlotsFlag = false;
    mDistPlotFlag = false;
    // 0: Log Normal
    // 1: Gaussian
    // 2: Custom
    if (value == mPara->Custom) //Custom Distribution
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
    QCustomPlot *customPlot = ui->customPlot_Csca;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(Csca)";
    else
        strNameY = "Csca";
    DisplayGraphData(event, customPlot, strNameX, strNameY);
}


//MouseOver PlotExtinctionCrossSection
void MainWindow::MouseOverPlotExtinctionCrossSection(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *customPlot = ui->customPlot_Cext;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(Cext)";
    else
        strNameY = "Cext";
    DisplayGraphData(event, customPlot, strNameX, strNameY);
}

//MouseOver PlotBackscatteringCrossSection
void MainWindow::MouseOverPlotBackscatteringCrossSection(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *customPlot = ui->customPlot_Cback;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(Cback)";
    else
        strNameY = "Cback";
    DisplayGraphData(event, customPlot, strNameX, strNameY);
}

//MouseOver PlotSizeParameter
void MainWindow::MouseOverPlotSizeParameter(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *customPlot = ui->customPlot_SizePara;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(Size Parameter)";
    else
        strNameY = "Size Parameter";
    DisplayGraphData(event, customPlot, strNameX, strNameY);
}

//MouseOver PlotMus
void MainWindow::MouseOverPlotMus(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *customPlot = ui->customPlot_Mus;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(μs)";
    else
        strNameY = "μs";
    DisplayGraphData(event, customPlot, strNameX, strNameY);
}

//MouseOver PlotMusp
void MainWindow::MouseOverPlotMusp(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *customPlot = ui->customPlot_Musp;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(μs')";
    else
        strNameY = "μs'";
    DisplayGraphData(event, customPlot, strNameX, strNameY);
}

//MouseOver PlotG
void MainWindow::MouseOverPlotG(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *customPlot = ui->customPlot_G;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(g)";
    else
        strNameY = "g";
    DisplayGraphData(event, customPlot, strNameX, strNameY);
}

//MouseOver PlotForwardBackward
void MainWindow::MouseOverPlotForwardBackward(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *customPlot = ui->customPlot_FB;
    customPlot->legend->setVisible(false);
    customPlot->replot();
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(%)";
    else
        strNameY = "%";
    DisplayGraphData(event, customPlot, strNameX, strNameY);
}

//MouseOver PlotS1S2
void MainWindow::MouseOverPlotS1S2(QMouseEvent *event)
{    
    QString strNameX = "Ang.";
    QString strNameY;
    QCustomPlot *customPlot = ui->customPlot_S1S2;
    customPlot->legend->setVisible(false);
    customPlot->replot();
    if (ui->radioButton_LogYAxis->isChecked())
    {
        if (ui->radioButton_S1->isChecked())
            strNameY = "Log(S1)";
        if (ui->radioButton_S2->isChecked())
            strNameY = "Log(S2)";
        if (ui->radioButton_S1S2->isChecked())
            strNameY = "Log(S)";
    }
    else
    {
        if (ui->radioButton_S1->isChecked())
            strNameY = "S1";
        if (ui->radioButton_S2->isChecked())
            strNameY = "S2";
        if (ui->radioButton_S1S2->isChecked())
            strNameY = "S";
    }
    DisplayGraphData(event, customPlot, strNameX, strNameY);
}

//MouseOver PlotPhaseFunctionLinear
void MainWindow::MouseOverPlotPhaseFunctionLinear(QMouseEvent *event)
{
    QString strNameX = "Ang.";
    QString strNameY;
    QCustomPlot *customPlot = ui->customPlot_PhaseFunctionLinear;
    customPlot->legend->setVisible(false);
    customPlot->replot();
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(Magnitude)";
    else
        strNameY = "Magnitude";
    DisplayGraphData(event, customPlot, strNameX, strNameY);
}

void MainWindow::MouseOverPlotPhaseFunctionPolar(QMouseEvent *event)
{
    QString strNameX = "Ang.";
    QString strNameY;
    QCustomPlot *customPlot = ui->customPlot_PhaseFunctionPolar;
    customPlot->legend->setVisible(false);
    customPlot->replot();
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(Magnitude)";
    else
        strNameY = "Magnitude";
    PlotData plot;
    DisplayPolarCurveData(event, customPlot, strNameX, strNameY);
}

//MouseOver PlotDistribution
void MainWindow::MouseOverPlotDistribution(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *customPlot = ui->customPlot_Distribution;
    QString strNameX = "Dia.";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(Ns)";
    else
        strNameY = "Ns";
    DisplayGraphData(event, customPlot, strNameX, strNameY);
}

//MouseOver Musp power Law
void MainWindow::MouseOverPlotMuspPowerLaw(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *customPlot = ui->customPlot_MuspPowerLaw;
    customPlot->legend->setVisible(false);
    customPlot->replot();
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(μs')";
    else
        strNameY = "μs'";
    DisplayGraphData(event, customPlot, strNameX, strNameY);
}

//Slot: SpinBoxF_valueChanged
void MainWindow::Slider_F_valueChanged(int value)
{
    ui->doubleSpinBox_F->blockSignals(true);
    ui->doubleSpinBox_F->setValue(0.01 * value);
    ui->doubleSpinBox_F->blockSignals(false);
    CalculateFittingAndDisplay();
}

//Slot: SpinBoxB_valueChanged
void MainWindow::Slider_B_valueChanged(int value)
{
    ui->doubleSpinBox_B->blockSignals(true);
    ui->doubleSpinBox_B->setValue(0.01 * value);
    ui->doubleSpinBox_B->blockSignals(false);
    CalculateFittingAndDisplay();
}

//Slot: DoubleSpinBox_F_valueChanged
void MainWindow::DoubleSpinBox_F_valueChanged(double value)
{
    ui->slider_F->blockSignals(true);
    ui->slider_F->setValue(100 * value);
    ui->slider_F->blockSignals(false);
    CalculateFittingAndDisplay();
}

//Slot: DoubleSpinBox_B_valueChanged
void MainWindow::DoubleSpinBox_B_valueChanged(double value)
{
    ui->slider_B->blockSignals(true);
    ui->slider_B->setValue(100 * value);
    ui->slider_B->blockSignals(false);
    CalculateFittingAndDisplay();
}

/******************************** Helper methods ***************************************/

//Update musp fitting plot
void MainWindow::CalculateFittingAndDisplay()
{
    mPara->bMie = ui->doubleSpinBox_B->value();
    mPara->fRay = ui->doubleSpinBox_F->value();
    UpdateMuspFitPlot();
    UpdateMuspFitErrorDisplay();
}

// Save QCustomPlot Plots
void MainWindow::SavePlot(QCustomPlot *customPlot, QString fileName)
{
    //Save Data
    QString filesTypes = tr("PNG (*.png);;JPEG (*.jpg);;BMP (*.bmp);;PDF (*.pdf)");

    QString name = QFileDialog::
        getSaveFileName(this, tr("Save Plot"),fileName,
                        filesTypes);
    if (name.isEmpty())
        return;
    else
    {
        if (name.endsWith(".png", Qt::CaseInsensitive))
            customPlot->savePng(name,QPrinter::HighResolution);
        if (name.endsWith(".jpg", Qt::CaseInsensitive))
            customPlot->saveJpg(name);
        if (name.endsWith(".bmp", Qt::CaseInsensitive))
            customPlot->saveBmp(name);
        if (name.endsWith(".pdf", Qt::CaseInsensitive))
            customPlot->savePdf(name);
        RememberLastDirectory(name);
    }
}

//Remember last directory
void MainWindow::RememberLastDirectory(QString fileName)
{
    int pos = fileName.lastIndexOf('/');
    QDir::setCurrent(fileName.left(pos));
}

//Display x and y values
void MainWindow::DisplayGraphData(QMouseEvent *event, QCustomPlot *customPlot,
                                  QString strNameX, QString strNameY)
{
    QCPAbstractPlottable *plottable = customPlot->plottableAt(event->position());
    if(plottable)
    {
        double x = customPlot->xAxis->pixelToCoord(event->position().x());
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
                    QToolTip::showText(event->globalPosition().toPoint(),
                                       tr("<table>"
                                          "<tr>"
                                          "<td>%L1:</td>" "<td>%L2</td>"
                                          "</tr>"
                                          "<tr>"
                                          "<td>%L3:</td>" "<td>%L4</td>"
                                          "</tr>"
                                          "</table>").arg(strNameX).arg(key).arg(strNameY).arg(value),
                                       customPlot, customPlot->rect());
                }
            }
        }
    }
    else
        QToolTip::hideText();
}

//Display x and y values
void MainWindow::DisplayPolarCurveData(QMouseEvent *event, QCustomPlot *customPlot,
                                       QString strNameX, QString strNameY)
{
    PlotData plot;

    QCPAbstractPlottable *plottable = customPlot->plottableAt(event->position());
    if(plottable)
    {
        double x = customPlot->xAxis->pixelToCoord(event->position().x());
        double *dx;

        QCPCurve *curve =  qobject_cast<QCPCurve*>(plottable);
        if (curve)
        {
            double key = 0;
            double value = 0;
            double rho = 0;
            double angle = 0;
            bool ok = false;
            double maxx = std::numeric_limits<double>::max();

            QCPCurveDataContainer::const_iterator begin = curve->data()->constBegin();
            QCPCurveDataContainer::const_iterator end = curve->data()->constEnd();

            unsigned int n = static_cast<unsigned int>(end-begin);
            if (n>0)
            {
                dx = new double[n];
                int index =0;
                for (QCPCurveDataContainer::const_iterator it=begin; it<end; it++)
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
                //For polar plot
                rho = sqrt(key*key + value*value);
                if (ui->radioButton_PhaseLog->isChecked())
                {
                    double logMin = plot.FindMinLogPolarPlot(mPara);
                    rho = pow(10, rho+logMin);
                }
                angle = 180*atan2(value,key)/M_PI;
                if (angle<0)
                    angle = angle + 360;

                if (ok)
                {
                    QToolTip::showText(event->globalPosition().toPoint(),
                                       tr("<table>"
                                          "<tr>"
                                          "<td>%L1:</td>" "<td>%L2</td>"
                                          "</tr>"
                                          "<tr>"
                                          "<td>%L3:</td>" "<td>%L4</td>"
                                          "</tr>"
                                          "</table>").arg(strNameX).arg(angle).arg(strNameY).arg(rho),
                                       customPlot, customPlot->rect());
                }
            }
        }
    }
    else
        QToolTip::hideText();
}
