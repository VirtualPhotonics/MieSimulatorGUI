/**********************************************************************
** All supporting functions to carry out MainWindow actions are listed
** in this file.
**********************************************************************/

#include "dialog/mainwindowsupport.h"
#include "dialog/plotdata.h"
#include <QMessageBox>

MainWindowSupport::MainWindowSupport(void)
{
}

// Initialize GUI
void MainWindowSupport::InitializeGUI(Ui_MainWindow *ui, Parameters *para)
{
    //Set variable limits
    ui->lineEdit_StartWL->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_EndWL->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_StepWL->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_Diameter->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_NumDen->setValidator( new QDoubleValidator(1e-50, 1e50, 12) );
    ui->lineEdit_VolFrac->setValidator( new QDoubleValidator(1e-50, 1e50, 12) );
    ui->lineEdit_MeanDiameter->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_StdDev->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_ScatRefReal->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_ScatRefImag->setValidator( new QDoubleValidator(-1e15, 1e15, 12) );
    ui->lineEdit_MedRef->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_NSphere->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );

    //Initial setting
    ui->label_ScatRefImag->setText("<font color=\"brown\">Sph.  Im.</font>");
    ui->comboBox_Distribution->blockSignals(true);
    ui->comboBox_Distribution->addItem("Log Normal");
    ui->comboBox_Distribution->addItem("Gaussian");
    ui->comboBox_Distribution->addItem("Custom");
    ui->comboBox_Distribution->setDisabled(true);
    ui->comboBox_Distribution->blockSignals(false);

    ui->lineEdit_StartWL->setText(QString::number(para->startWavel));
    ui->lineEdit_EndWL->setText(QString::number(para->endWavel));
    ui->lineEdit_StepWL->setText(QString::number(para->stepWavel));

    ui->lineEdit_Diameter->setText(QString::number(para->meanRadius));
    ui->lineEdit_ScatRefReal->setText(QString::number(para->scatRefReal));
    ui->lineEdit_ScatRefImag->setText(QString::number(para->scatRefImag));
    ui->lineEdit_MedRef->setText(QString::number(para->medRef));

    ui->lineEdit_MeanDiameter->setText(QString::number(para->meanRadius));
    ui->lineEdit_StdDev->setText(QString::number(para->stdDev));
    ui->lineEdit_NSphere->setText(QString::number(para->nRadius));
    ui->lineEdit_NumDen->setText(QString::number(para->sphNumDensity));
    ui->lineEdit_VolFrac->setText(QString::number(para->volFraction));
}

//Reset widgets
void MainWindowSupport::SetWidgets(Ui_MainWindow *ui, Parameters *para)
{
    bool falseFlag, trueFlag;

    if (ui->radioButton_MonoDisperse->isChecked())
    {
        falseFlag = false;
        trueFlag = true;
    }
	else
    {
        falseFlag = true;
        trueFlag = false;
    }

    //Mono Dispere Parameters
    ui->lineEdit_Diameter->setDisabled(falseFlag);
    ui->label_Diameter->setDisabled(falseFlag);

    //Poly Dispere Parameters
    ui->comboBox_Distribution->setDisabled(trueFlag);
    ui->pushButton_ShowDistributionAndCustom->setDisabled(trueFlag);
    ui->radioButton_LinearXAxis->setDisabled(trueFlag);
    ui->radioButton_LogXAxis->setDisabled(trueFlag);
    ui->label_DistXAxis->setDisabled(trueFlag);

    if (ui->comboBox_Distribution->currentIndex() != para->Custom)  // if not "Custom"""""""""
    {
        ui->lineEdit_MeanDiameter->setDisabled(trueFlag);
        ui->label_MeanDiameter->setDisabled(trueFlag);
        ui->lineEdit_StdDev->setDisabled(trueFlag);
        ui->label_StdDev->setDisabled(trueFlag);
        ui->lineEdit_NSphere->setDisabled(trueFlag);
        ui->label_NSphere->setDisabled(trueFlag);

        //Enable or disable according to the selection
        if (ui->radioButton_NumDen->isChecked())
        {
            ui->lineEdit_NumDen->setEnabled(true);
            ui->lineEdit_VolFrac->setDisabled(true);
        }

        if (ui->radioButton_VolFrac->isChecked())
        {
            ui->lineEdit_NumDen->setDisabled(true);
            ui->lineEdit_VolFrac->setEnabled(true);
        }
    }

    //Disable following widgets after initialization or different selection
    ui->slider_ConcPercentChange->setValue(0);
    ui->slider_ConcPercentChange->setDisabled(true);
    ui->label_ConcPercent->setDisabled(true);
    ui->label_ActualConcPercent->setDisabled(true);
    ui->label_CurrentTotNumDen->setVisible(false);
    ui->label_CurrentA->setText("");
    ui->label_CurrentMSE->setText("");
    ui->label_Progress->setText("");
}

//Copy input data to 'para' variable
void MainWindowSupport::LoadInputData(Ui_MainWindow *ui, Parameters *para)
{    
    para->startWavel = ui->lineEdit_StartWL->text().toDouble();
    para->endWavel = ui->lineEdit_EndWL->text().toDouble();
    para->stepWavel = ui->lineEdit_StepWL->text().toDouble();
    para->scatRefReal = ui->lineEdit_ScatRefReal->text().toDouble();
    para->scatRefImag = ui->lineEdit_ScatRefImag->text().toDouble();
    para->medRef = ui->lineEdit_MedRef->text().toDouble();

    if (ui->radioButton_NumDen->isChecked())
    {
        para->sphNumDensity = ui->lineEdit_NumDen->text().toDouble();
    }
    if (ui->radioButton_VolFrac->isChecked())
    {
        para->volFraction = ui->lineEdit_VolFrac->text().toDouble();
    }

    if (ui->radioButton_MonoDisperse->isChecked())
    {
        para->meanRadius = 0.5 * ui->lineEdit_Diameter->text().toDouble();
        para->nRadius = 1;
    }
    if ((ui->radioButton_PolyDisperse->isChecked()) &&
        (ui->comboBox_Distribution->currentIndex() != para->Custom))
    {
        para->meanRadius = 0.5 * ui->lineEdit_MeanDiameter->text().toDouble();
        para->stdDev = ui->lineEdit_StdDev->text().toDouble();
        para->nRadius = ui->lineEdit_NSphere->text().toUInt();
    }
    para->fRay = ui->doubleSpinBox_F->text().toDouble();
    para->bMie = ui->doubleSpinBox_B->text().toDouble();

    if (ui->radioButton_FittingComplex->isChecked())
    {
        para->fittingComplex = true;
    }
    if (ui->radioButton_FittingSimple->isChecked())
    {
        para->fittingComplex = false;
    }
    if (ui->radioButton_RefWavel500->isChecked())
    {
        para->refWavel = 500.0;
        para->refWavelIdx = para->wavel500;
    }
    if (ui->radioButton_RefWavel600->isChecked())
    {
        para->refWavel = 600.0;
        para->refWavelIdx = para->wavel600;
    }
    if (ui->radioButton_RefWavel700->isChecked())
    {
        para->refWavel = 700.0;
        para->refWavelIdx = para->wavel700;
    }
    if (ui->radioButton_RefWavel800->isChecked())
    {
        para->refWavel = 800.0;
        para->refWavelIdx = para->wavel800;
    }
    if (ui->radioButton_RefWavel900->isChecked())
    {
        para->refWavel = 900.0;
        para->refWavelIdx = para->wavel900;
    }
    if (ui->radioButton_RefWavel1000->isChecked())
    {
        para->refWavel = 1000.0;
        para->refWavelIdx = para->wavel1000;
    }
    SetWavelengthSliders(ui);
}

//Initialize dynamic arrays.
void MainWindowSupport::InitializeArrays(Ui_MainWindow *ui, Parameters *para, bool *arrayFlag)
{
    if (ui->radioButton_Phase_DTheta0_1->isChecked())
        para->nTheta = 1801;    //180/(1801-1) = 0.1degree step
    if (ui->radioButton_Phase_DTheta0_5->isChecked())
        para->nTheta = 361;    //180/(361-1) = 0.5degree step

    para->nWavel = static_cast<unsigned int>(floor(para->endWavel - para->startWavel) / para->stepWavel) + 1;
    para->wavelArray = new double [para->nWavel];
    for (unsigned int i=0; i<para->nWavel; ++i)
    {
        para->wavelArray[i] = para->startWavel + i * para->stepWavel;
    }

    para->phaseFunctionAve = new double *[para->nWavel];
    for (unsigned int i=0; i<para->nWavel; i++)
    {
        para->phaseFunctionAve[i] = new double [para->nTheta];
    }

    para->phaseFunctionPara = new double *[para->nWavel];
    for (unsigned int i=0; i<para->nWavel; i++)
    {
        para->phaseFunctionPara[i] = new double [para->nTheta];
    }

    para->phaseFunctionPerp = new double *[para->nWavel];
    for (unsigned int i=0; i<para->nWavel; i++)
    {
        para->phaseFunctionPerp[i] = new double [para->nTheta];
    }

    para->S1 = new std::complex<double> *[para->nWavel];
    for (unsigned int i=0; i<para->nWavel; i++)
    {
        para->S1[i] = new std::complex<double> [para->nTheta];
    }

    para->S2 = new std::complex<double> *[para->nWavel];
    for (unsigned int i=0; i<para->nWavel; i++)
    {
        para->S2[i] = new std::complex<double> [para->nTheta];
    }

    para->cSca = new double [para->nWavel];
    para->cExt = new double [para->nWavel];
    para->cBack = new double [para->nWavel];
    para->sizePara = new double [para->nWavel];
    para->mus = new double [para->nWavel];
    para->g = new double [para->nWavel];
    para->forward = new double [para->nWavel];
    para->backward = new double [para->nWavel];

    *arrayFlag = true;
}

//Delete dynamic arrays
void MainWindowSupport::DeleteArrays(Parameters *para, bool *arrayFlag)
{
    delete[] para->wavelArray;
    delete[] para->cSca;
    delete[] para->cExt;
    delete[] para->cBack;
    delete[] para->sizePara;
    delete[] para->mus;
    delete[] para->g;
    delete[] para->forward;
    delete[] para->backward;

    for (unsigned int i=0; i<para->nWavel; i++)
    {
        delete [] para->phaseFunctionAve[i];
    }
    for (unsigned int i=0; i<para->nWavel; i++)
    {
        delete [] para->phaseFunctionPara[i];
    }
    for (unsigned int i=0; i<para->nWavel; i++)
    {
        delete [] para->phaseFunctionPerp[i];
    }
    for (unsigned int i=0; i<para->nWavel; i++)
    {
        delete [] para->S1[i];
    }
    for (unsigned int i=0; i<para->nWavel; i++)
    {
        delete [] para->S2[i];
    }

    *arrayFlag = false;
}

//Set wavelenght sliders
void MainWindowSupport::SetWavelengthSliders(Ui_MainWindow *ui)
{
    double startWL = ui->lineEdit_StartWL->text().toDouble();
    double endWL = ui->lineEdit_EndWL->text().toDouble();
    double stepWL = ui->lineEdit_StepWL->text().toDouble();
    int nWL = static_cast<int>(floor(endWL - startWL)/stepWL) + 1;

    ui->slider_WL_PFPolar->setMinimum(0);
    ui->slider_WL_PFPolar->setMaximum(nWL-1);
    ui->slider_WL_PFPolar->setSingleStep(1);
    ui->label_CurrentWL_PFPolar->setText(QString::number(startWL));

    ui->slider_WL_PFLinear->setMinimum(0);
    ui->slider_WL_PFLinear->setMaximum(nWL-1);
    ui->slider_WL_PFLinear->setSingleStep(1);
    ui->label_CurrentWL_PFLinear->setText(QString::number(startWL));

    ui->slider_WL_S1S2->setMinimum(0);
    ui->slider_WL_S1S2->setMaximum(nWL-1);
    ui->slider_WL_S1S2->setSingleStep(1);
    ui->label_CurrentWL_S1S2->setText(QString::number(startWL));
}

// Run Mono disperse distribution
void MainWindowSupport::ProcessMonoDisperse(Ui_MainWindow *ui, Parameters *para)
{
    mCalc = new Calculate();
    PlotData plot;

    //Run Mie simulation for radius para->radArray[i]
    mCalc->DoSimulation(ui->label_Progress, para);

    //Get musp at reference wavelengths
    mCalc->ComputeMuspAtRefWavel(para);

    //Assign values and plot
    plot.SetupPolarPlotForData(ui, para);
    plot.AssignValuesPhaseFunctionPolarPlot(ui,para);
    plot.AssignValuesPhaseFunctionLinearPlot(ui,para);
    plot.AssignValuesS1S2Plot(ui, para);
    plot.AssignValuesOtherPlots(ui, para);

    //Check independent/dependent scattering
    double clearanceToWavelength, sizeParameter, volFraction, criticalWavelength;
    QString strRegime;
    bool flagVolFrac = ui->radioButton_VolFrac->isChecked();
    if(mCalc->CheckIndependentScattering(para, clearanceToWavelength, sizeParameter, volFraction,
                                          criticalWavelength, strRegime, flagVolFrac))
    {
        DisplayScatteringRegimeWarning(clearanceToWavelength, sizeParameter, volFraction, criticalWavelength, strRegime);
    }
}

// Run Poly disperse distribution
void MainWindowSupport::ProcessPolyDisperse(Ui_MainWindow *ui, Parameters *para)
{
    mCalc = new Calculate();
    PlotData plot;

    //Run Mie simulation for radius para->radArray[i]
    mCalc->DoSimulation(ui->label_Progress, para);

    //Get musp at reference wavelengths
    mCalc->ComputeMuspAtRefWavel(para);

    //Assign values and plot
    plot.SetupPolarPlotForData(ui, para);
    plot.AssignValuesPhaseFunctionPolarPlot(ui,para);
    plot.AssignValuesPhaseFunctionLinearPlot(ui,para);
    plot.AssignValuesS1S2Plot(ui, para);
    plot.AssignValuesOtherPlots(ui, para);

    //Check independent/dependent scattering
    double clearanceToWavelength, sizeParameter, volFraction, criticalWavelength;
    QString strRegime;
    bool flagVolFrac = ui->radioButton_VolFrac->isChecked();
    if(mCalc->CheckIndependentScattering(para, clearanceToWavelength, sizeParameter, volFraction,
                                          criticalWavelength, strRegime, flagVolFrac))
    {
        DisplayScatteringRegimeWarning(clearanceToWavelength, sizeParameter, volFraction, criticalWavelength, strRegime);
    }
}

//Sphere distribution in polydisperse
void MainWindowSupport::ProcessDistribution(Ui_MainWindow *ui, Parameters *para, unsigned int distIndex)
{
    PlotData plot;
    if (distIndex != 2)
    {
        bool flagVolOrConc = true;
        para->radArray = new double [para->nRadius];
        para->numDensityArray = new double [para->nRadius];
        para->scatRefRealArray = new double [para->nRadius];
        para->scatRefImagArray = new double [para->nRadius];
        para->medRefArray = new double [para->nRadius];

        //Find size of spheres
        if (ui->radioButton_VolFrac->isChecked())
        {
            flagVolOrConc = true;
        }
        if (ui->radioButton_NumDen->isChecked())
        {
            flagVolOrConc = false;
        }
        mCalc->DiameterRangeSetting(para, distIndex);
        mCalc->SetSphereRadiusAndRefIndex(para, distIndex, flagVolOrConc);
    }
    plot.AssignValuesDistributionPlot(ui, para);
}

//Disable Enable Real and Imaginary buttons
void MainWindowSupport::DisableEnableRealImagButtons(Ui_MainWindow *ui)
{
    if (ui->radioButton_LinearYAxis->isChecked())
    {
        ui->radioButton_S1S2_Real->setDisabled(false);
        ui->radioButton_S1S2_Imag->setDisabled(false);
    }
    if (ui->radioButton_LogYAxis->isChecked())
    {
        ui->radioButton_S1S2_Real->setDisabled(true);
        ui->radioButton_S1S2_Imag->setDisabled(true);
        ui->radioButton_S1S2_Abs->setChecked(true);
    }
}

// Disable widgets during simulation
void MainWindowSupport::DisableWidgetsDuringSimulation(Ui_MainWindow *ui, Parameters *para, bool flag)
{
    //disable widgets during simulation
    ui->pushButton_RunSimulation->setDisabled(flag);
    ui->radioButton_MonoDisperse->setDisabled(flag);
    ui->radioButton_PolyDisperse->setDisabled(flag);
    ui->radioButton_NumDen->setDisabled(flag);
    ui->radioButton_VolFrac->setDisabled(flag);
    ui->radioButton_LinearYAxis->setDisabled(flag);
    ui->radioButton_LogYAxis->setDisabled(flag);
    ui->pushButton_DisplayData->setDisabled(flag);
    ui->pushButton_SaveData->setDisabled(flag);
    ui->pushButton_BestFit->setDisabled(flag);
    ui->label_ConcPercent->setDisabled(flag);
    ui->label_ActualConcPercent->setDisabled(flag);
    ui->slider_ConcPercentChange->setDisabled(flag);
    ui->slider_B->setDisabled(flag);
    ui->doubleSpinBox_B->setDisabled(flag);
    if (!para->fittingComplex)
    {
        ui->slider_F->setDisabled(true);
        ui->doubleSpinBox_F->setDisabled(true);
    }
    if (para->fittingComplex)
    {
        ui->slider_F->setDisabled(flag);
        ui->doubleSpinBox_F->setDisabled(flag);
    }
}

// Disable/Enable widgets during "Custom" Polydisperse data
void MainWindowSupport::DisableWidgetsDuringCustomPolyDisperseData(Ui_MainWindow *ui, bool flag)
{    
    //disable widgets during simulation
    ui->lineEdit_MeanDiameter->setDisabled(flag);
    ui->lineEdit_StdDev->setDisabled(flag);
    ui->lineEdit_NSphere->setDisabled(flag);
    ui->lineEdit_ScatRefReal->setDisabled(flag);
    ui->lineEdit_ScatRefImag->setDisabled(flag);
    if (!flag)
    {
        if (ui->radioButton_NumDen->isChecked())
        {
            ui->lineEdit_NumDen->setDisabled(flag);
            ui->lineEdit_VolFrac->setEnabled(flag);
        }

        if (ui->radioButton_VolFrac->isChecked())
        {
            ui->lineEdit_NumDen->setEnabled(flag);
            ui->lineEdit_VolFrac->setDisabled(flag);
        }
    }
    else
    {
        ui->lineEdit_NumDen->setDisabled(flag);
        ui->lineEdit_VolFrac->setDisabled(flag);
    }
    ui->label_MeanDiameter->setDisabled(flag);
    ui->label_StdDev->setDisabled(flag);
    ui->label_NSphere->setDisabled(flag);
    ui->label_ScatRefReal->setDisabled(flag);
    ui->label_ScatRefImag->setDisabled(flag);
    ui->radioButton_NumDen->setDisabled(flag);
    ui->radioButton_VolFrac->setDisabled(flag);
    if (flag)
    {
        ui->pushButton_ShowDistributionAndCustom->setText("Load Custom Data");
    }
    else
    {
        ui->pushButton_ShowDistributionAndCustom->setText("Show Distribution");
    }
}

//Read "Custom" (PolyDisperse) data from a file
void MainWindowSupport::ReadCustomData(Parameters *para, QString fileName, bool *dataValidFlag)
{
    int count = -1;
    int badLines = -1;

    QFile file(fileName);
    QString line;
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Error");
        msgBox.setIcon(QMessageBox::Critical);
        msgBox.setText("Unable to open the file");
        msgBox.exec();
        return;
    }
    QTextStream in(&file);
    do
    {
        line = in.readLine();
        QStringList list = line.split(QRegularExpression(",|;|\t"));

        if (!(list.size() == 4 || list.size() == 5))
        {
            badLines++;
        }
        else
        {
            if (list.size()==4)
            {
                bool check1, check2, check3, check4;
                list.at(0).toDouble(&check1);
                list.at(1).toDouble(&check2);
                list.at(2).toDouble(&check3);
                list.at(3).toDouble(&check4);
                if (!check1 || !check2 || !check3 || !check4)
                {
                    badLines++;
                }
            }
            if (list.size()==5)
            {
                bool check1, check2, check3, check4, check5;
                list.at(0).toDouble(&check1);
                list.at(1).toDouble(&check2);
                list.at(2).toDouble(&check3);
                list.at(3).toDouble(&check4);
                list.at(4).toDouble(&check5);
                if (!check1 || !check2 || !check3 || !check4 || !check5)
                {
                    badLines++;
                }
            }
        }
        count++;
    }while (!line.isNull());

    file.close();

    if ((badLines > 0) || (count < 1))
    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Error");
        msgBox.setIcon(QMessageBox::Critical);
        msgBox.setText("Invalid data format!");
        msgBox.exec();
        *dataValidFlag = true;
    }
    else
    {
        para->nRadius = static_cast<unsigned int>(count);
        para->radArray = new double [para->nRadius];
        para->numDensityArray = new double [para->nRadius];
        para->scatRefRealArray = new double [para->nRadius];
        para->scatRefImagArray = new double [para->nRadius];
        para->medRefArray = new double [para->nRadius];

        file.open(QIODevice::ReadOnly | QIODevice::Text);
        QTextStream in(&file);
        int idx =0;
        double sumRad =0;
        double minRad = 1e100;
        double maxRad = 0;
        do
        {
            line = in.readLine();
            QStringList list = line.split(QRegularExpression(",|;|\t"));
            if (list.size() == 4)
            {
                 para->radArray[idx] = list.at(0).toDouble()/2.0;
                 para->numDensityArray[idx] = list.at(1).toDouble();
                 para->scatRefRealArray[idx] = list.at(2).toDouble();
                 para->scatRefImagArray[idx] = list.at(3).toDouble();
                 para->medRefArray[idx] = para->medRef; //use default
                 sumRad += para->radArray[idx];
                 if (para->radArray[idx] >maxRad)
                 {
                     maxRad = para->radArray[idx];
                 }
                 if (para->radArray[idx] < minRad)
                 {
                     minRad = para->radArray[idx];
                 }
                 idx++;
            }
            if (list.size() == 5)
            {
                 para->radArray[idx] = list.at(0).toDouble()/2.0;
                 para->numDensityArray[idx] = list.at(1).toDouble();
                 para->scatRefRealArray[idx] = list.at(2).toDouble();
                 para->scatRefImagArray[idx] = list.at(3).toDouble();
                 para->medRefArray[idx] = list.at(4).toDouble();
                 sumRad += para->radArray[idx];
                 if (para->radArray[idx] >maxRad)
                 {
                     maxRad = para->radArray[idx];
                 }
                 if (para->radArray[idx] < minRad)
                 {
                     minRad = para->radArray[idx];
                 }
                 idx++;
            }
        }while (!line.isNull());

        file.close();

        para->meanRadius = sumRad/para->nRadius;
        para->minRadius = minRad;
        para->maxRadius = maxRad;
        //Assign this values to pass the sanity check
        para->stdDev = 1.0;
    }
}

void MainWindowSupport::DisplayScatteringRegimeWarning(double clearanceToWavelength, double sizeParameter,
                                                       double volFraction, double criticalWavelength,
                                                       QString strRegime)
{
    QString strTienDorlen = (clearanceToWavelength > 0.5)
                                ? "<font color='green'>Independent</font>"
                                : "<font color='red'>Dependent</font>";

    QString strGaly = "<font color='red'>Dependent</font>";

    QString msg = QString(
                      "<b>Dependent Scattering Warning</b><br>"
                      "The current configuration deviates from the independent scattering regime. "
                      "Since MieSimulatorGUI is best suited for dilute mixtures, results at this "
                      "concentration should be interpreted with caution.<br><br>"
                      "<b>Independent or Dependent:</b><br>"
                      "• Per <b>Tien and Drolen (1987):</b> %1<br>"
                      "• Per <b>Galy et al. (2020):</b> %2<br><br>"
                      "<b>Details:</b><br>"
                      "• <b>Regime:</b> %3<br>"
                      "• <b>Volume Fraction:</b> %4<br>"
                      "• <b>Critical Wavelength:</b> %5<br>"
                      "• <b>Size Parameter:</b> %6<br>"
                      "• <b>Clearance/Wavelength Ratio:</b> %7"
                      )
                      .arg(strTienDorlen)
                      .arg(strGaly)
                      .arg(strRegime)
                      .arg(volFraction, 0, 'g', 6)
                      .arg(criticalWavelength, 0, 'g', 6)
                      .arg(sizeParameter, 0, 'g', 6)
                      .arg(clearanceToWavelength, 0, 'g', 6);

    DisplayWarning(msg);
}

// Display warning message
void MainWindowSupport::DisplayWarning(QString warningMessage)
{
    QMessageBox msgBoxWarning;
    msgBoxWarning.setWindowTitle("Warning");
    msgBoxWarning.setIcon(QMessageBox::Warning);
    msgBoxWarning.setText(warningMessage);
    msgBoxWarning.exec();
}

