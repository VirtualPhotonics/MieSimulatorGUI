//**********************************************************************
//** All supporting functions to carry out MainWindow actions are listed
//** in this file.
//**********************************************************************

#include "mainwindowsupport.h"


MainWindowSupport::MainWindowSupport(void)
{
}

// Initialize GUI
void MainWindowSupport::InitializeGUI(Ui_MainWindow *ui)
{
    //Set variable limits
    ui->lineEdit_startWavel->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_endWavel->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_stepWavel->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_diameter->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_conc_mm3->setValidator( new QDoubleValidator(1e-50, 1e50, 12) );
    ui->lineEdit_volFrac->setValidator( new QDoubleValidator(1e-50, 1e50, 12) );
    ui->lineEdit_meanDiameter->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_stdDev->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_scatRefReal->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_scatRefImag->setValidator( new QDoubleValidator(-1e15, 1e15, 12) );
    ui->lineEdit_medRef->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );
    ui->lineEdit_nSphere->setValidator( new QDoubleValidator(1e-15, 1e15, 12) );

    //Initial setting
    ui->label_scatRefImag->setText("<font color=\"brown\">Sph.  Im.</font>");
    ui->comboBox_distribution->addItem("Log Normal");
    ui->comboBox_distribution->addItem("Gaussian");
    ui->comboBox_distribution->addItem("Custom");
    ui->comboBox_distribution->setDisabled(true);

    ui->lineEdit_startWavel->setText("600");
    ui->lineEdit_endWavel->setText("1000");
    ui->lineEdit_stepWavel->setText("10");

    ui->lineEdit_diameter->setText("0.5");
    ui->lineEdit_scatRefReal->setText("1.377");
    ui->lineEdit_scatRefImag->setText("0.0");
    ui->lineEdit_medRef->setText("1.333");

    ui->lineEdit_meanDiameter->setText("0.5");
    ui->lineEdit_stdDev->setText("0.25");
    ui->lineEdit_nSphere->setText("11");
    ui->lineEdit_conc_mm3->setText("1e8");
    ui->lineEdit_volFrac->setText("0.1");
}

//Reset widgets
void MainWindowSupport::SetWidgets(Ui_MainWindow *ui)
{
    bool falseFlag, trueFlag;

    if (ui->radioButton_monoDisperse->isChecked())
    {
        falseFlag = false;
        trueFlag = true;
    }
	else
    {
        falseFlag = true;
        trueFlag = false;
    }
    //Mono Dispere parameters
    ui->lineEdit_diameter->setDisabled(falseFlag);
    ui->label_diameter->setDisabled(falseFlag);
    //Poly Dispere parameters
    ui->comboBox_distribution->setDisabled(trueFlag);
    ui->pushButton_showDistributionAndCustom->setDisabled(trueFlag);
    ui->radioButton_linearXaxis->setDisabled(trueFlag);
    ui->radioButton_logXaxis->setDisabled(trueFlag);
    ui->label_distXAxis->setDisabled(trueFlag);

    if (ui->comboBox_distribution->currentIndex() !=2)  // if not "Custom"
    {
        ui->lineEdit_meanDiameter->setDisabled(trueFlag);
        ui->label_meanDiameter->setDisabled(trueFlag);
        ui->lineEdit_stdDev->setDisabled(trueFlag);
        ui->label_stdDev->setDisabled(trueFlag);
        ui->lineEdit_nSphere->setDisabled(trueFlag);
        ui->label_nSphere->setDisabled(trueFlag);

        //Enable or disable according to the selection
        if (ui->radioButton_conc_mm3->isChecked())
        {
            ui->lineEdit_conc_mm3->setDisabled(false);
            ui->lineEdit_volFrac->setDisabled(true);
        }

        if (ui->radioButton_volFrac->isChecked())
        {
            ui->lineEdit_conc_mm3->setDisabled(true);
            ui->lineEdit_volFrac->setDisabled(false);
        }
    }

    //Disable following widgets after initialization or different selection
    ui->slider_concPercentChange->setValue(0);
    ui->slider_concPercentChange->setDisabled(true);
    ui->label_concPercent->setDisabled(true);
    ui->label_actualConcPercent->setDisabled(true);
    ui->label_currentTotNumDen->setVisible(false);
    ui->label_currentA->setText("");
    ui->label_currentMse->setText("");
    ui->label_progress->setText("");
}

//Copy input data to 'para' variable
void MainWindowSupport::LoadInputData(Ui_MainWindow *ui, parameters *para)
{    
    para->startWavel = ui->lineEdit_startWavel->text().toDouble();
    para->endWavel = ui->lineEdit_endWavel->text().toDouble();
    para->stepWavel = ui->lineEdit_stepWavel->text().toDouble();
    para->scatRefReal = ui->lineEdit_scatRefReal->text().toDouble();
    para->scatRefImag = ui->lineEdit_scatRefImag->text().toDouble();
    para->medRef = ui->lineEdit_medRef->text().toDouble();
    if (ui->radioButton_conc_mm3->isChecked())
        para->sphNumDensity = ui->lineEdit_conc_mm3->text().toDouble();
    if (ui->radioButton_volFrac->isChecked())
        para->volFraction = ui->lineEdit_volFrac->text().toDouble();

    if (ui->radioButton_monoDisperse->isChecked())
    {
        para->meanRadius = 0.5 * ui->lineEdit_diameter->text().toDouble();
        para->nRadius = 1;
        if (ui->radioButton_volFrac->isChecked())
        {
            double volume = 4.0 * M_PI *para->meanRadius * para->meanRadius * para->meanRadius / 3.0;
            para->sphNumDensity = para->volFraction * 1e9 /volume ;
        }
    }
    if ((ui->radioButton_polyDisperse->isChecked()) && (ui->comboBox_distribution->currentIndex() != 2))
    {
        para->meanRadius = 0.5 * ui->lineEdit_meanDiameter->text().toDouble();
        para->stdDev = ui->lineEdit_stdDev->text().toDouble();
        para->nRadius = ui->lineEdit_nSphere->text().toUInt();
    }
    para->fRay = ui->doubleSpinBox_powerLaw_f->text().toDouble();
    para->bMie = ui->doubleSpinBox_powerLaw_b->text().toDouble();

    para->fittingComplex = true;
    para->refWavel = 1000.0;

    SetWavelengthSliders(ui);
}

//Initialize dynamic arrays.
void MainWindowSupport::InitializeArrays(Ui_MainWindow *ui, parameters *para, bool *arrayFlag)
{
    para->minTheta = 0;
    para->maxTheta = M_PI;
    if (ui->radioButton_dThetaStep0_1->isChecked())
        para->nTheta = 1801;    //180/(1801-1) = 0.1degree step
    if (ui->radioButton_dThetaStep0_5->isChecked())
        para->nTheta = 361;    //180/(361-1) = 0.5degree step

    para->stepTheta = (para->maxTheta - para->minTheta)/static_cast<double>(para->nTheta - 1);
    para->nWavel = static_cast<unsigned int>(floor(para->endWavel - para->startWavel) / para->stepWavel) + 1;
    para->wavelArray = new double [para->nWavel];
    for (unsigned int i=0; i<para->nWavel; ++i)
        para->wavelArray[i] = para->startWavel + i * para->stepWavel;

    para->phaseFunctionAve = new double *[para->nWavel];
    for (unsigned int i=0; i<para->nWavel; i++)
        para->phaseFunctionAve[i] = new double [para->nTheta];

    para->phaseFunctionPara = new double *[para->nWavel];
    for (unsigned int i=0; i<para->nWavel; i++)
        para->phaseFunctionPara[i] = new double [para->nTheta];

    para->phaseFunctionPerp = new double *[para->nWavel];
    for (unsigned int i=0; i<para->nWavel; i++)
        para->phaseFunctionPerp[i] = new double [para->nTheta];

    para->S1 = new std::complex<double> *[para->nWavel];
    for (unsigned int i=0; i<para->nWavel; i++)
        para->S1[i] = new std::complex<double> [para->nTheta];

    para->S2 = new std::complex<double> *[para->nWavel];
    for (unsigned int i=0; i<para->nWavel; i++)
        para->S2[i] = new std::complex<double> [para->nTheta];

    para->cSca = new double [para->nWavel];
    para->cExt = new double [para->nWavel];
    para->cBack = new double [para->nWavel];
    para->SizePara = new double [para->nWavel];
    para->mus = new double [para->nWavel];
    para->g = new double [para->nWavel];
    para->forward = new double [para->nWavel];
    para->backward = new double [para->nWavel];

    *arrayFlag = true;
}

//Delete dynamic arrays
void MainWindowSupport::DeleteArrays(parameters *para, bool *arrayFlag)
{
    delete[] para->wavelArray;
    delete[] para->cSca;
    delete[] para->cExt;
    delete[] para->cBack;
    delete[] para->SizePara;
    delete[] para->mus;
    delete[] para->g;
    delete[] para->forward;
    delete[] para->backward;

    for (unsigned int i=0; i<para->nWavel; i++)
        delete [] para->phaseFunctionAve[i];
    for (unsigned int i=0; i<para->nWavel; i++)
        delete [] para->phaseFunctionPara[i];
    for (unsigned int i=0; i<para->nWavel; i++)
        delete [] para->phaseFunctionPerp[i];
    for (unsigned int i=0; i<para->nWavel; i++)
        delete [] para->S1[i];
    for (unsigned int i=0; i<para->nWavel; i++)
        delete [] para->S2[i];    

    *arrayFlag = false;
}

//Set wavelenght sliders
void MainWindowSupport::SetWavelengthSliders(Ui_MainWindow *ui)
{
    double startWL = ui->lineEdit_startWavel->text().toDouble();
    double endWL = ui->lineEdit_endWavel->text().toDouble();
    double stepWL = ui->lineEdit_stepWavel->text().toDouble();
    int nWL = static_cast<int>(floor(endWL - startWL)/stepWL) + 1;

    ui->slider_pFunctionPolar_wavel->setMinimum(0);
    ui->slider_pFunctionPolar_wavel->setMaximum(nWL-1);
    ui->slider_pFunctionPolar_wavel->setSingleStep(1);
    ui->label_pFunctionPolar_curWavel->setText(QString::number(startWL));

    ui->slider_pFunctionLinear_wavel->setMinimum(0);
    ui->slider_pFunctionLinear_wavel->setMaximum(nWL-1);
    ui->slider_pFunctionLinear_wavel->setSingleStep(1);
    ui->label_pFunctionLinear_curWavel->setText(QString::number(startWL));

    ui->slider_s1s2_wavel->setMinimum(0);
    ui->slider_s1s2_wavel->setMaximum(nWL-1);
    ui->slider_s1s2_wavel->setSingleStep(1);
    ui->label_s1s2_curWavel->setText(QString::number(startWL));
}

// Run Mono disperse distribution
void MainWindowSupport::ProcessMonoDisperse(Ui_MainWindow *ui, parameters *para)
{
    mCalc = new calculate();
    PlotData plot;

    //Run Mie simulation for radius para->radArray[i]
    mCalc->DoSimulation(ui, para);

    //Get musp at reference wavelengths
    mCalc->ComputeMuspAtRefWavel(para);

    //Assign values and plot    
    plot.AssignValuesPhaseFunctionPolarPlot(ui,para);
    plot.AssignValuesPhaseFunctionLinearPlot(ui,para);
    plot.AssignValuesS1S2Plot(ui, para);
    plot.AssignValuesAllOtherPlots(ui, para);
}

// Run Poly disperse distribution
void MainWindowSupport::ProcessPolyDisperse(Ui_MainWindow *ui, parameters *para)
{
    mCalc = new calculate();
    PlotData plot;

    //Run Mie simulation for radius para->radArray[i]
    mCalc->DoSimulation(ui, para);

    //Get musp at reference wavelengths
    mCalc->ComputeMuspAtRefWavel(para);

    //Assign values and plot
    plot.AssignValuesPhaseFunctionPolarPlot(ui,para);
    plot.AssignValuesPhaseFunctionLinearPlot(ui,para);
    plot.AssignValuesS1S2Plot(ui, para);
    plot.AssignValuesAllOtherPlots(ui, para);
}

//Sphere distribution in poly disperse
void MainWindowSupport::ProcessDistribution(Ui_MainWindow *ui, parameters *para, unsigned int distIndex)
{
    PlotData plot;
    if (distIndex !=2)
    {
        bool flagVolOrConc = true;
        para->radArray = new double [para->nRadius];
        para->numDensityArray = new double [para->nRadius];
        para->scatRefRealArray = new double [para->nRadius];
        para->scatRefImagArray = new double [para->nRadius];

        //Find size of spheres
        if (ui->radioButton_volFrac->isChecked())
            flagVolOrConc = true;
        if (ui->radioButton_conc_mm3->isChecked())
            flagVolOrConc = false;
        mCalc->DiameterRangeSetting(para, distIndex);
        mCalc->SetSphereRadiusAndRefIndex(para, distIndex, flagVolOrConc);
    }
    plot.AssignValuesDistributionPlot(ui, para);
}

//Disable Enable Real and Imaginary buttons
void MainWindowSupport::DisableEnableRealImagButtons(Ui_MainWindow *ui)
{
    if (ui->radioButton_linearYaxis->isChecked())
    {
        ui->radioButton_s1s2_real->setDisabled(false);
        ui->radioButton_s1s2_imag->setDisabled(false);
    }
    if (ui->radioButton_logYaxis->isChecked())
    {
        ui->radioButton_s1s2_real->setDisabled(true);
        ui->radioButton_s1s2_imag->setDisabled(true);
        ui->radioButton_s1s2_abs->setChecked(true);
    }
}

// Disable widgets during simulation
void MainWindowSupport::DisableWidgetsDuringSimulation(Ui_MainWindow *ui, parameters *para, bool flag)
{
    //disable widgets during simulation
    ui->pushButton_runSimulation->setDisabled(flag);
    ui->radioButton_monoDisperse->setDisabled(flag);
    ui->radioButton_polyDisperse->setDisabled(flag);
    ui->radioButton_conc_mm3->setDisabled(flag);
    ui->radioButton_volFrac->setDisabled(flag);
    ui->radioButton_linearYaxis->setDisabled(flag);
    ui->radioButton_logYaxis->setDisabled(flag);
    ui->pushButton_displayData->setDisabled(flag);
    ui->pushButton_saveData->setDisabled(flag);
    ui->pushButton_bestFit->setDisabled(flag);
    ui->label_concPercent->setDisabled(flag);
    ui->label_actualConcPercent->setDisabled(flag);
    ui->slider_concPercentChange->setDisabled(flag);
    ui->qwtslider_powerLaw_f->setDisabled(flag);
    ui->qwtslider_powerLaw_b->setDisabled(flag);
    if (!para->fittingComplex)
    {
        ui->qwtslider_powerLaw_f->setDisabled(true);
        ui->doubleSpinBox_powerLaw_f->setDisabled(true);
    }
    else
    {
        ui->qwtslider_powerLaw_f->setDisabled(flag);
        ui->doubleSpinBox_powerLaw_f->setDisabled(flag);
    }
}

// Disable/Enable widgets during "Custom" Polydisperse data
void MainWindowSupport::DisableWidgetsDuringCustomPolyDisperseData(Ui_MainWindow *ui, bool flag)
{
    //disable widgets during simulation
    ui->lineEdit_meanDiameter->setDisabled(flag);
    ui->lineEdit_stdDev->setDisabled(flag);
    ui->lineEdit_nSphere->setDisabled(flag);
    ui->lineEdit_scatRefReal->setDisabled(flag);
    ui->lineEdit_scatRefImag->setDisabled(flag);
    ui->lineEdit_conc_mm3->setDisabled(flag);
    ui->lineEdit_volFrac->setDisabled(flag);
    ui->label_meanDiameter->setDisabled(flag);
    ui->label_stdDev->setDisabled(flag);
    ui->label_nSphere->setDisabled(flag);
    ui->label_scatRefReal->setDisabled(flag);
    ui->label_scatRefImag->setDisabled(flag);
    ui->radioButton_conc_mm3->setDisabled(flag);
    ui->radioButton_volFrac->setDisabled(flag);
    if (flag)
        ui->pushButton_showDistributionAndCustom->setText("Load Custom Data");
    else
        ui->pushButton_showDistributionAndCustom->setText("Show Distribution");
}

//Read "Custom" (PolyDisperse) data from a file
void MainWindowSupport::ReadCustomData(parameters *para, QString fileName, bool *dataValidFlag)
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
        QStringList list = line.split(QRegExp(",|;|\t"));
        if (list.size()!=4)
            badLines++;
        count++;
    }while (!line.isNull());
    file.close();

    if ((badLines>0) || (count < 1))
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


        file.open(QIODevice::ReadOnly | QIODevice::Text);
        QTextStream in(&file);
        int idx =0;
        double sumRad =0;
        double minRad = 1e100;
        double maxRad = 0;
        do
        {
            line = in.readLine();
            QStringList list = line.split(QRegExp(",|;|\t"));
            if (list.size()==4)
            {
                 para->radArray[idx] = list.at(0).toDouble()/2.0;
                 para->numDensityArray[idx] = list.at(1).toDouble();
                 para->scatRefRealArray[idx] = list.at(2).toDouble();
                 para->scatRefImagArray[idx] = list.at(3).toDouble();
                 sumRad += para->radArray[idx];
                 if (para->radArray[idx] >maxRad)
                     maxRad = para->radArray[idx];
                 if (para->radArray[idx] < minRad)
                     minRad = para->radArray[idx];
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

//*************************  Sanity check  *********************************************
//Input parameter check
bool MainWindowSupport::CheckInputParameters(Ui_MainWindow *ui, parameters *para)
{
    QMessageBox msgBox, msgBoxWarn;
    msgBox.setWindowTitle("Error");
    msgBox.setIcon(QMessageBox::Critical);
    msgBoxWarn.setWindowTitle("Warning");
    msgBoxWarn.setIcon(QMessageBox::Warning);

    if ((para->scatRefReal <= 0.0) || (para->medRef <= 0.0))
    {
        msgBox.setText("Refractive index cannot be zero");
        msgBox.exec();
        return 1;
    }
    if ((para->scatRefReal/para->medRef == 1.0))
    {
        msgBox.setText("Relative refractive index cannot be 1.0");
        msgBox.exec();
        return 1;
    }
    double m = para->scatRefReal / para->medRef;
    if ((m < 0.05) || (m > 5.0))
    {
        msgBoxWarn.setText("Unrealistic relative refractive index.");
        msgBoxWarn.setInformativeText("Check sphere and medium refractive index values");
        msgBoxWarn.exec();
        return 1;
    }
    if (para->scatRefImag >= 5.0)
    {
        msgBox.setText("Imaginary refractive index must be negative or less than 5.0.");
        msgBox.exec();
        return 1;
    }
    //Avoid zeros
    if ((para->startWavel <= 0.0) || (para->endWavel <= 0.0))
    {
        msgBox.setText("The starting or ending wavelength cannot be zero");
        msgBox.exec();
        return 1;
    }
    //Avoid zeros
    if (para->stepWavel <= 0.0)
    {
        msgBox.setText("Wavelength step cannot be zero");
        msgBox.exec();
        return 1;
    }
    if (para->startWavel < 50.0)
    {
        msgBox.setText("Minimum wavlength is 50nm");
        msgBox.exec();
        return 1;
    }
    if (para->endWavel > 3000.0)
    {
        msgBox.setText("Maximum wavlength is 3000nm");
        msgBox.exec();
        return 1;
    }
    if (para->endWavel - para->startWavel < 0.0)
    {
        msgBox.setText("The starting wavelength is greater than the ending wavelength");
        msgBox.exec();
        return 1;
    }
    if (ui->radioButton_conc_mm3->isChecked())
    {
        if (para->sphNumDensity <= 0.0)
        {
            msgBox.setText("Sphere concentration cannot be zero");
            msgBox.exec();
            return 1;
        }
    }
    if (ui->radioButton_volFrac->isChecked())
    {
        if (para->volFraction <= 0.0)
        {
            msgBox.setText("Volume Fraction cannot be zero");
            msgBox.exec();
            return 1;
        }

        if (para->volFraction >= 1.0)
        {
            msgBox.setText("Volume Fraction cannot exceed 1.0");
            msgBox.exec();
            return 1;
        }
    }
    if ((para->meanRadius < 0.00005) ||(para->meanRadius >150))
    {
        msgBox.setText("Diameter is out of range!");
        msgBox.setInformativeText("Enter a value between 0.0001μm and 300μm");
        msgBox.exec();
        return 1;
    }
    if (ui->radioButton_monoDisperse->isChecked())
    {
        if (ui->radioButton_conc_mm3->isChecked())
        {
            double volume = 4.0 * M_PI *para->meanRadius * para->meanRadius * para->meanRadius / 3.0;
            if (para->sphNumDensity*volume >= 1e9)
            {
                msgBoxWarn.setText("Concentration x Sphere Volume exceeds 1mm³.");
                msgBoxWarn.setInformativeText("Reduce Concentration.");
                msgBoxWarn.exec();
                return 1;
            }
        }
    }
    return 0;
}

//Poly disperse distribution parameter check
bool MainWindowSupport::CheckDistribution(Ui_MainWindow *ui, parameters *para)
{
    QMessageBox msgBox;
    msgBox.setWindowTitle("Error");

    if (para->stdDev == 0.0)
    {
        msgBox.setText("Standard Deviation is zero.");
        msgBox.setInformativeText("Use 'Mono Disperse'.");
        msgBox.exec();
        return 1;
    }
    if(ui->comboBox_distribution->currentIndex() == 0)
    {
        if (para->stdDev > 3.0)
        {
            msgBox.setIcon(QMessageBox::Critical);
            msgBox.setText("Large standard deviation provides an abnormal Log Normal distribution.");
            msgBox.setInformativeText("The limit was set to 3.0μm.");
            msgBox.exec();
            return 1;
        }
        if (para->stdDev < 1e-5)
        {
            msgBox.setIcon(QMessageBox::Critical);
            msgBox.setText("The standard deviation is too small.");
            msgBox.setInformativeText("Use 'Mono Disperse'.");
            msgBox.exec();
            return 1;
        }
    }

    if(ui->comboBox_distribution->currentIndex() == 1)
    {
        if (para->stdDev > 50.0)
        {
            msgBox.setIcon(QMessageBox::Critical);
            msgBox.setText("The standard deviation is too large.");
            msgBox.setInformativeText("The limit was set to 50.0μm.");
            msgBox.exec();
            return 1;
        }
        if (para->stdDev < 1e-8)
        {
            msgBox.setIcon(QMessageBox::Critical);
            msgBox.setText("The standard deviation is too small.");
            msgBox.setInformativeText("Use 'Mono Disperse'.");
            msgBox.exec();
            return 1;
        }
    }
    if (para->nRadius == 1)
    {
        msgBox.setText("Discrete sphere size is 1. ");
        msgBox.setInformativeText("Use 'Mono Disperse'.");
        msgBox.exec();
        return 1;
    }
    if ((para->nRadius < 2.0) ||(para->nRadius >101.0))
    {
        msgBox.setText("Number of sphere sizes is out of range.");
        msgBox.setInformativeText("Enter a value between 2 and 101.");
        msgBox.exec();
        return 1;
    }
    if ((para->meanRadius < 0.0005) ||(para->meanRadius >25))
    {
        msgBox.setText("Diameter is out of range.");
        msgBox.setInformativeText("Enter a value between 0.001μm and 50μm.");
        msgBox.exec();
        return 1;
    }
    if (para->stdDev/para->meanRadius < 1.999e-5)
    {
        msgBox.setIcon(QMessageBox::Critical);
        msgBox.setText("Standard deviation to mean diameter ratio is smaller than 1e-5.");
        msgBox.setInformativeText("Use 'Mono Disperse'");
        msgBox.exec();
        return 1;
    }
    return 0;
}
