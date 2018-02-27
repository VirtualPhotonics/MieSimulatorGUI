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
    connect(ui->customPlot_ScatCross, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotScatCross(QMouseEvent*)));
    connect(ui->customPlot_Mus, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotMus(QMouseEvent*)));
    connect(ui->customPlot_Musp, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotMusp(QMouseEvent*)));
    connect(ui->customPlot_G, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotG(QMouseEvent*)));
    connect(ui->customPlot_FB, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotForwardBackward(QMouseEvent*)));
    connect(ui->customPlot_S1S2, SIGNAL(mouseMove(QMouseEvent*)), SLOT(MouseOverPlotS1S2(QMouseEvent*)));
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
    plot.InitializePolarPlot(ui,mPara);
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
            support.DisableWidgetsDuringSimulation(ui, true);
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
                    support.ProcessDistribution(ui, mPara, ui->comboBox_Distribution->currentIndex());
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
            support.DisableWidgetsDuringSimulation(ui, false);
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
                support.ProcessDistribution(ui, mPara, ui->comboBox_Distribution->currentIndex());
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
        display->show();
        display->DisplayData(ui,  mPara);
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
            support.DisableWidgetsDuringSimulation(ui, true);
            ui->label_CurrentMSE->setText("<font color=\"red\">Wait!...    </font>");
            ui->label_CurrentA->setText("");
            util.Delay();

            //Calculation
            cal.CalculatePowerLawAutoFit(mPara);
            plot.AssignValuesMuspPowerLawPlots(ui,mPara);
            fRay = mPara->fRay;
            bMie = mPara->bMie;
            ui->doubleSpinBox_F->setValue(fRay);
            ui->doubleSpinBox_B->setValue(bMie);
            ui->label_CurrentMSE->setText("<font color=\"red\">M.S. Error = </font>"+QString::number(mPara->muspFittingError,'g',6));
            ui->label_CurrentA->setText("<font color=\"blue\">A = </font>"+QString::number(mPara->fittedA,'g',6));

            //Enable Widgets
            support.DisableWidgetsDuringSimulation(ui, false);
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
        plot.AssignValuesAllOtherPlots(ui, mPara);
        plot.AssignValuesS1S2Plot(ui, mPara);
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
        plot.AssignValuesAllOtherPlots(ui, mPara);
        plot.AssignValuesS1S2Plot(ui, mPara);
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
    ui->tabWidget_PhaseFunction->setTabText(1,"S1 and S2");
    ui->tabWidget_ScatCross->setTabText(1,"Scattering Cross Section");
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
    ui->tabWidget_PhaseFunction->setTabText(1,"Ave. S1 and S2");
    ui->tabWidget_ScatCross->setTabText(1,"Ave. Scattering Cross Section");
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
    UpdatePhaseFunctionPlot();
}

//radioButton_PhaseLog_clicked: Select Log Scale
void MainWindow::on_radioButton_PhaseLog_clicked()
{
    UpdatePhaseFunctionPlot();
}

//radioButton_PhaseAverage_clicked: Select average
void MainWindow::on_radioButton_PhaseAverage_clicked()
{
    UpdatePhaseFunctionPlot();
}

//radioButton_PhasePara_clicked: Select Parallel component
void MainWindow::on_radioButton_PhasePara_clicked()
{
    UpdatePhaseFunctionPlot();
}

//radioButton_PhasePerp_clicked: Select Perpendicular component
void MainWindow::on_radioButton_PhasePerp_clicked()
{
    UpdatePhaseFunctionPlot();
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

//on_sliderWL_PF_valueChanged: Change Phase Function Wavelenghth Slider
void MainWindow::on_slider_PF_WL_valueChanged(int value)
{
   int wavel = mPara->startWavel + value*mPara->stepWavel;
   ui->label_CurrentWL_PF->setText(QString::number(wavel));
   UpdatePhaseFunctionPlot();
}

//on_sliderWL_S1S2_valueChanged: Change S1/S2 Wavelenghth Slider
void MainWindow::on_slider_S1S2_WL_valueChanged(int value)
{
   int wavel = mPara->startWavel + value*mPara->stepWavel;
   ui->label_CurrentWL_S1S2->setText(QString::number(wavel));
   UpdateS1S2Plot();
}

//SpinBoxF_valueChanged: Calculate fitting parameters and display
void MainWindow::on_doubleSpinBox_F_valueChanged(double arg1)
{
    ui->qwtslider_F->setValue(arg1);
    mPara->fRay = arg1;
    UpdateMuspFitPlot();
}

//SpinBoxB_valueChanged: Calculate fitting parameters and display
void MainWindow::on_doubleSpinBox_B_valueChanged(double arg1)
{
    ui->qwtslider_B->setValue(arg1);
    mPara->bMie = arg1;
    UpdateMuspFitPlot();
}


/********************************** Update Functions **********************************/

//Process Musp Fit
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

//Process S1S2 plot
void MainWindow::UpdateS1S2Plot()
{
    if (mOtherPlotsFlag)
    {
        PlotData plot;
        plot.AssignValuesS1S2Plot(ui, mPara);
    }
}

//Process phase function plot
void MainWindow::UpdatePhaseFunctionPlot()
{
    if (mOtherPlotsFlag)
    {
        PlotData plot;
        plot.AssignValuesPolarPlot(ui,mPara);
    }
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

//MouseOver PlotScatCross
void MainWindow::MouseOverPlotScatCross(QMouseEvent *event)
{
    QString strNameY;
    QCustomPlot *curPlot = ui->customPlot_ScatCross;
    QString strNameX = "WL";
    if (ui->radioButton_LogYAxis->isChecked())
        strNameY = "Log(Csca)";
    else
        strNameY = "Csca";
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
            strNameY = "Log(S2)";\
    }
    else
    {
        if (ui->radioButton_S1->isChecked())
            strNameY = "S1";
        if (ui->radioButton_S2->isChecked())
            strNameY = "S2";
    } \
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

            int n = end-begin;
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
