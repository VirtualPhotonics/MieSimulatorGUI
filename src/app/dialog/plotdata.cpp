/**********************************************************************
** All QCustomPlot asignments and plotting actions are listed here.
**********************************************************************/

#include "plotdata.h"
#include "calc/utilities.h"
#include "lib/qcustomplot.h"

PlotData::PlotData(void)
{
}

PlotData::~PlotData(void)
{
}

//Clear and Initialize Plots
void PlotData::ClearPlots(Ui_MainWindow *ui)
{
    InitialSetupDistributionPlot(ui);
    InitialSetupPhaseFunctionPolarPlot(ui);
    InitialSetupPhaseFunctionLinearPlot(ui);
    InitialSetupS1S2Plot(ui);
    InitialSetupMuspPowerLawFit(ui);
    InitialSetupOtherPlots(ui);
}

//Initialize distribution plot
void PlotData::InitialSetupDistributionPlot(Ui_MainWindow *ui)
{
    const double BUFFER_PERCENT = 0.05;
    double minX, maxX;
    minX = (1 - BUFFER_PERCENT) * (ui->lineEdit_Diameter->text().toDouble());
    maxX = (1 + BUFFER_PERCENT) * (ui->lineEdit_Diameter->text().toDouble());

    auto customPlot = ui->customPlot_Distribution;

    RemovePlotables(customPlot);
    InitialSetupPlot(customPlot, "Diameter (μm)", "Ns  (/vol. of 1mm³)", minX, maxX);
}

//Initialize phase function linear plot
void PlotData::InitialSetupPhaseFunctionLinearPlot(Ui_MainWindow *ui)
{    
    double minX = -180;
    double maxX = 180;

    auto customPlot = ui->customPlot_PhaseFunctionLinear;
    InitialSetupPlot(customPlot, "Angle (deg.)", "Magnitude", minX, maxX);

    QVector<double> ticks{-180, -135, -90, -45, 0, 45, 90, 135, 180};
    QVector<QString> labels;
    for (double tick : ticks)
        labels << QString::number(tick);
    QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
    textTicker->addTicks(ticks, labels);
    customPlot->xAxis->setTicker(textTicker);
    customPlot->replot();
}

//Initialize phase function polar plot
void PlotData::InitialSetupPhaseFunctionPolarPlot(Ui_MainWindow *ui)
{
    auto customPlot = ui->customPlot_PhaseFunctionPolar;
    mPolarMinRadius = 0.0;
    mPolarMaxRadius = 1.0;
    DrawPolarPlotGrid(customPlot, true); //true for linear
    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
}

//Initialize S1,S2
void PlotData::InitialSetupS1S2Plot(Ui_MainWindow *ui)
{
    auto customPlot = ui->customPlot_S1S2;
    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);

    // remove previous graphs (if any)
    RemoveGraphs(customPlot);

    customPlot->xAxis->setLabel("Angle (deg.)");
    customPlot->xAxis->setRange(0,180);
    QVector<double> ticks;
    QVector<QString> labels;
    ticks << 0 << 45 << 90 << 135 << 180;
    labels << "0" << "45" << "90" << "135" << "180";
    QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
    textTicker->addTicks(ticks, labels);
    customPlot->xAxis->setTicker(textTicker);
    customPlot->yAxis->setLabel("S1");
    customPlot->yAxis->setRange(0, 5);
    customPlot->replot();
}

//Initialize all Other (Scattering cross section, g and Musp) plots
void PlotData::InitialSetupOtherPlots(Ui_MainWindow *ui)
{
    double minX = ui->lineEdit_StartWL->text().toDouble();
    double maxX = ui->lineEdit_EndWL->text().toDouble();

    InitialSetupPlot(ui->customPlot_Csca, "Wavelength (nm)", "Scattering Cross Section (μm²)", minX, maxX);
    InitialSetupPlot(ui->customPlot_Cext, "Wavelength (nm)", "Extinction Cross Section (μm²)", minX, maxX);
    InitialSetupPlot(ui->customPlot_Cback, "Wavelength (nm)", "Backscattering Cross Section (μm²)", minX, maxX);
    InitialSetupPlot(ui->customPlot_SizePara, "Wavelength (nm)", "Size Parameter", minX, maxX);
    InitialSetupPlot(ui->customPlot_Mus, "Wavelength (nm)", "μs (mmˉˡ)", minX, maxX);
    InitialSetupPlot(ui->customPlot_G, "Wavelength (nm)", "g (Average Cosine of phase function)", minX, maxX);
    InitialSetupPlot(ui->customPlot_FB, "Wavelength (nm)", "Forward & Backward Scattering %", minX, maxX);
    InitialSetupPlot(ui->customPlot_Musp, "Wavelength (nm)", "μs' (mmˉˡ)", minX, maxX);
}

//Plot Musp for fitting
void PlotData::InitialSetupMuspPowerLawFit(Ui_MainWindow *ui)
{   
    double minX = ui->lineEdit_StartWL->text().toDouble();
    double maxX = ui->lineEdit_EndWL->text().toDouble();

    InitialSetupPlot(ui->customPlot_MuspPowerLaw, "Wavelength (nm)", "μs' (mmˉˡ)", minX, maxX);

    ui->doubleSpinBox_B->setValue(4);
    ui->doubleSpinBox_F->setValue(0);
}

// Prepare polar plot based on data
void PlotData::SetupPolarPlotForData(Ui_MainWindow *ui, Parameters *para )
{
    auto customPlot = ui->customPlot_PhaseFunctionPolar;

    // Set initial boundaries.
    double maxPolarRadial = -1e100;
    double minPolarRadial = 1e100;
    int nTheta = static_cast<int>(para->nTheta);
    int nWavel = static_cast<int>(para->nWavel);

    for (int i = 0; i < nWavel; i++)
    {
        for (int j = 0; j < nTheta; j++)
        {
            double currentPara = para->phaseFunctionPara[i][j];
            double currentPerp = para->phaseFunctionPerp[i][j];
            double currentAve = para->phaseFunctionAve[i][j];

            if (ui->radioButton_PhaseLog->isChecked())
            {
                currentPara = log10(currentPara);
                currentPerp = log10(currentPerp);
                currentAve = log10(currentAve);
            }

            maxPolarRadial = std::max({maxPolarRadial, currentPara, currentPerp, currentAve});
            if (ui->radioButton_PhaseLog->isChecked())
                minPolarRadial = std::min({minPolarRadial, currentPara, currentPerp, currentAve});
        }
    }

    if (ui->radioButton_PhaseLog->isChecked())
    {
        mPolarMinRadius = pow(10, floor(minPolarRadial));
        mPolarMaxRadius = pow(10, maxPolarRadial);
    } else {
        mPolarMinRadius = 0.0;
        mPolarMaxRadius = maxPolarRadial;
    }
    DrawPolarPlotGrid(customPlot, !ui->radioButton_PhaseLog->isChecked());
    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
}

//Assign values for distribution plot
void PlotData::AssignValuesDistributionPlot(Ui_MainWindow *ui, Parameters* para)
{
    QVector<double> xDist(static_cast<int>(para->nRadius));
    QVector<double> yDist(static_cast<int>(para->nRadius));
    double tiny = 1e-100;  //add a small number before log calculation to avoid NaN
    double margin = (1.0 + ui->slider_ConcPercentChange->value()/200.0);

    double totNumDen = 0.0;
    for (int i=0; i<static_cast<int>(para->nRadius); i++)
    {
        xDist[i] = 2.0 * para->radArray[i];
        yDist[i] = para->numDensityArray[i]*margin;
        if (ui->radioButton_LogYAxis->isChecked())
            yDist[i] = log10(yDist[i] + tiny);
        totNumDen += para->numDensityArray[i]*margin;
    }
    if (ui->radioButton_MonoDisperse->isChecked())
        ui->label_CurrentTotNumDen->setVisible(false);
    if (ui->radioButton_PolyDisperse->isChecked())
    {
        ui->label_CurrentTotNumDen->setVisible(true);
        ui->label_CurrentTotNumDen->setText("<font color=\"red\">Σ Ns: "+QString::number(totNumDen)+" in vol. of 1mm³");
    }
    PlotDistribution(ui, para, xDist,yDist);
}

//Assign values for Phase Function Linear plot
void PlotData::AssignValuesPhaseFunctionLinearPlot(Ui_MainWindow *ui, Parameters *para)
{
    int nTheta = static_cast<int>(para->nTheta);
    int totalSize = 2 * nTheta - 1;
    bool flagThetaNegPosOrPos = true;  //theta: -180 to +180
    QVector<double> phaseFuncPara(totalSize);
    QVector<double> phaseFuncPerp(totalSize);
    QVector<double> phaseFuncAve(totalSize);
    QVector<double> theta(totalSize);

    int indexWL = ui->slider_WL_PFLinear->value();
    ui->label_CurrentWL_PFLinear->setText(QString::number(para->startWavel + indexWL*para->stepWavel));

    RearrangePhaseFunctionData(para, theta,phaseFuncPara,phaseFuncPerp,phaseFuncAve,indexWL,
                               !ui->radioButton_LogYAxis->isChecked(), flagThetaNegPosOrPos);
    PlotPhaseFunctionLinear(ui, theta, phaseFuncPara, phaseFuncPerp, phaseFuncAve);
}

//Assign values for phase function polar plot
void PlotData::AssignValuesPhaseFunctionPolarPlot(Ui_MainWindow *ui, Parameters *para)
{
    int nTheta = static_cast<int>(para->nTheta);
    int totalSize = 2 * nTheta - 1;
    bool flagThetaNegPosOrPos = false;  //theta: 0 to 360
    QVector<double> phaseFuncPara(totalSize);
    QVector<double> phaseFuncPerp(totalSize);
    QVector<double> phaseFuncAve(totalSize);
    QVector<double> theta(totalSize);

    int indexWL = ui->slider_WL_PFPolar->value();
    ui->label_CurrentWL_PFPolar->setText(QString::number(para->startWavel + indexWL*para->stepWavel));

    RearrangePhaseFunctionData(para, theta,phaseFuncPara,phaseFuncPerp,phaseFuncAve,indexWL,
                               !ui->radioButton_PhaseLog->isChecked(), flagThetaNegPosOrPos);
    PlotPhaseFunctionPolar(ui, theta, phaseFuncPara, phaseFuncPerp, phaseFuncAve, totalSize);

    if (para->nWavel==1)
        ui->slider_WL_PFPolar->setDisabled(true);
    else
        ui->slider_WL_PFPolar->setDisabled(false);
}

//Assign values for S1,S2 plot
void PlotData::AssignValuesS1S2Plot(Ui_MainWindow *ui, Parameters *para)
{
    int nTheta = static_cast<int>(para->nTheta);
    QVector<double> S1(nTheta);
    QVector<double> S2(nTheta);
    QVector<double> theta(nTheta);

    int indexWL = ui->slider_WL_S1S2->value();
    ui->label_CurrentWL_S1S2->setText(QString::number(para->startWavel + indexWL*para->stepWavel));

    std::complex<double> *tempS1 = nullptr;
    std::complex<double> *tempS2 = nullptr;

    tempS1 = para->S1[indexWL];
    tempS2 = para->S2[indexWL];

    for (int i=0; i<nTheta; i++)
    {
        theta[i] = 180.0 * i /(para->nTheta-1);

        if (ui->radioButton_LogYAxis->isChecked())
        {
            S1[i] = log10(sqrt(tempS1[i].real()*tempS1[i].real() + tempS1[i].imag()*tempS1[i].imag()));
            S2[i] = log10(sqrt(tempS2[i].real()*tempS2[i].real() + tempS2[i].imag()*tempS2[i].imag()));
        }
        else
        {
            if (ui->radioButton_S1S2_Abs->isChecked())
            {
                S1[i] = sqrt(tempS1[i].real()*tempS1[i].real() + tempS1[i].imag()*tempS1[i].imag());
                S2[i] = sqrt(tempS2[i].real()*tempS2[i].real() + tempS2[i].imag()*tempS2[i].imag());
            }
            if (ui->radioButton_S1S2_Real->isChecked())
            {
                S1[i] = tempS1[i].real();
                S2[i] = tempS2[i].real();
            }
            if (ui->radioButton_S1S2_Imag->isChecked())
            {
                S1[i] = tempS1[i].imag();
                S2[i] = tempS2[i].imag();
            }
        }
    }
    PlotS1S2(ui, para, theta, S1, S2);

}

//Assign values for Musp power law plots
void PlotData::AssignValuesMuspPowerLawPlots(Ui_MainWindow *ui, Parameters* para)
{
    QVector<double> x(static_cast<int>(para->nWavel));
    QVector<double> yMusp(static_cast<int>(para->nWavel));
    QVector<double> yFit(static_cast<int>(para->nWavel));
    double tiny = 1e-100;  //add a small number before log calculation to avoid NaN
    double tempError = 0.0;
    double error;
    double mus, fitA;
    double margin = (1.0 + ui->slider_ConcPercentChange->value() /200.0);

    for (int i=0; i<static_cast<int>(para->nWavel); i++)
    {
        mus = para->mus[i]*margin;  // set mus according to conc slider value

        x[i] = para->wavelArray[i];
        yMusp[i] = mus * (1.0 - para->g[i]);          //reduced scattering coefficient

        //Steve L Jacques,"Optical properties of biological tissues: a review" Phys. Med & Bio. 58(2013) R37-R61.
        //wavelength λ is normalized by a reference wavelength, 500 nm or 1000nm
        fitA = para->muspAtRefWavel[para->refWavelIdx] *margin;
        if (para->fittingComplex)
            yFit[i] = fitA *(para->fRay*pow(x[i]/para->refWavel, -4.0) + (1-para->fRay)*pow(x[i]/para->refWavel, -para->bMie));
        else
            yFit[i] = fitA *pow(x[i]/para->refWavel, -para->bMie);   // A(lambda/lambdaRef)^-b

        error = yFit[i] - yMusp[i];
        tempError += error*error;

        if (ui->radioButton_LogYAxis->isChecked())
        {
            yFit[i] = log10(yFit[i]+tiny);
            yMusp[i] = log10((mus * (1.0 - para->g[i])) + tiny);    //reduced scattering coefficient
        }
    }
    para->muspFittingError = tempError/para->nWavel;
    PlotMuspPowerLaw(ui, x,yMusp,yFit);
}

//Assign values for all Other (Scattering cross section, g and Musp) plots
void PlotData::AssignValuesOtherPlots(Ui_MainWindow *ui, Parameters* para)
{
    QVector<double> x(static_cast<int>(para->nWavel));
    QVector<double> yCsca(static_cast<int>(para->nWavel));
    QVector<double> yCext(static_cast<int>(para->nWavel));
    QVector<double> yCback(static_cast<int>(para->nWavel));
    QVector<double> ySizePara(static_cast<int>(para->nWavel));
    QVector<double> yG(static_cast<int>(para->nWavel));
    QVector<double> yMus(static_cast<int>(para->nWavel));
    QVector<double> yMusp(static_cast<int>(para->nWavel));
    QVector<double> yF(static_cast<int>(para->nWavel));
    QVector<double> yB(static_cast<int>(para->nWavel));

    double tiny = 1e-100;  //add a small number before log calculation to avoid NaN
    double mus;
    double margin = (1.0 + ui->slider_ConcPercentChange->value() /200.0);

    int fbLegendCheckLocation = static_cast<int>(0.75*para->nWavel);
    bool fbLegnedFlag;
    double fbLimit =75;   //75%

    for (int i=0; i<static_cast<int>(para->nWavel); i++)
    {
        mus = para->mus[i]*margin;  // set mus according to conc slider value

        x[i] = para->wavelArray[i];
        yCsca[i] = para->cSca[i];    //scattering cross section
        yCext[i] = para->cExt[i];    //extinction cross section
        yCback[i] = para->cBack[i];  //backscattering cross section
        yG[i] = para->g[i];          //g - average cosine of phase function
        yMus[i] = mus;               //scattering coefficient
        yMusp[i] = mus * (1.0 - para->g[i]);          //reduced scattering coefficient
        yF[i] = para->forward[i];
        yB[i] = para->backward[i];
        if (ui->radioButton_MonoDisperse->isChecked())
            ySizePara[i] = para->SizePara[i];

        if (ui->radioButton_LogYAxis->isChecked())
        {
            yCsca[i] = log10(para->cSca[i] + tiny);    //scattering cross section
            yCext[i] = log10(para->cExt[i] + tiny);    //extinction cross section
            yCback[i] = log10(para->cBack[i] + tiny);  //backscattering cross section
            yG[i] = log10(para->g[i] +tiny);           //g - average cosine of phase function
            yMus[i] = log10(mus + tiny);              //scattering coefficient
            yMusp[i] = log10((mus * (1.0 - para->g[i])) + tiny);    //reduced scattering coefficient
            yF[i] = log10(para->forward[i] + tiny);
            yB[i] = log10(para->backward[i] + tiny);
            if (ui->radioButton_MonoDisperse->isChecked())
                ySizePara[i] = log10(para->SizePara[i] + tiny);
            fbLimit = 2;   //100%
        }
    }
    PlotG(ui, x, yG);
    if (yF[fbLegendCheckLocation]>fbLimit)
        fbLegnedFlag = true;
    else
        fbLegnedFlag = false;
    PlotForwardBackward(ui, x, yF, yB, fbLegnedFlag);

    PlotScatteringCrossSection(ui, x, yCsca);
    PlotExtinctionCrossSection(ui, x, yCext);
    PlotBackscatteringCrossSection(ui, x, yCback);
    PlotMus(ui, x, yMus);
    PlotMusp(ui, x, yMusp);
    PlotMuspCurveForPowerLaw(ui, x, yMusp);
    if (ui->radioButton_MonoDisperse->isChecked())
        PlotSizeParameter(ui, x, ySizePara);
}

//Find MinLog value for MouseOver Polar Plot
double PlotData::FindMinLogPolarPlot(Parameters *para)
{
    Utilities util;

    double minPolarRadial = 1e100;
    int nTheta = static_cast<int>(para->nTheta);
    int nWavel = static_cast<int>(para->nWavel);

    for (int i = 0; i < nWavel; i++)
    {
        QVector<double> yPara, yPerp, yAve;

        for (int j = 0; j < nTheta; j++)
        {
            yPara.append(log10(para->phaseFunctionPara[i][j]));
            yPerp.append(log10(para->phaseFunctionPerp[i][j]));
            yAve.append(log10(para->phaseFunctionAve[i][j]));
        }

        double currentMin = util.FindMinMax(yPara, yPerp, yAve, true);  //true for Min
        if (currentMin < minPolarRadial)
            minPolarRadial = currentMin;
    }
    return floor(minPolarRadial);
}

//plot distribution plot
void PlotData::PlotDistribution(Ui_MainWindow *ui, Parameters *para, QVector<double> x, QVector<double> yDist)
{    
    const double MAX_Y_BUFFER = 1.1;
    const double X_BUFFER_MIN = 0.92;
    const double X_BUFFER_MAX = 1.08;
    const double BAR_WIDTH_FACTOR = 0.05;
    const double SHIFT_SCALE = 3.0;
    const double BAR_WIDTH_SCALE = 0.1;

    double minX, maxX, minY, maxY;
    minY = 0;
    maxY = yDist[0];
    for(int i = 0; i<static_cast<int>(para->nRadius); i++)
    {
        if (maxY < yDist[i])
            maxY = yDist[i];
    }
    maxY *= MAX_Y_BUFFER;

    auto customPlot = ui->customPlot_Distribution;

    // Clear all plottables, which includes graphs and bars
    customPlot->clearPlottables();

    //Add new graph
    PlotSingleGraph(customPlot, x, yDist, Qt::red, "Distirbution", 0, 2);

    QCPBars *bar = new QCPBars(customPlot->xAxis, customPlot->yAxis);
    bar->setData(x, yDist);
    customPlot->yAxis->setRange(minY, maxY);
    bar->setPen(Qt::NoPen);
    bar->setBrush(QColor(0,0,0,255));

    if (ui->radioButton_MonoDisperse->isChecked())
    {
        minX = X_BUFFER_MIN * 2.0 * para->minRadius;
        maxX = X_BUFFER_MAX * 2.0 * para->maxRadius;
        bar->setWidth(BAR_WIDTH_FACTOR*(maxX -minX));
        customPlot->xAxis->setRange(minX, maxX);
        customPlot->xAxis->ticker()->setTickCount(3);
    }
    if (ui->radioButton_PolyDisperse->isChecked())
    {
        if (para->maxRadius == para->minRadius)      //Array of spherical scatterers with similar radius
        {
            minX = para->minRadius;
            maxX = SHIFT_SCALE * para->maxRadius;
            bar->setWidth(BAR_WIDTH_SCALE*para->maxRadius);
        }
        else
        {
            double shift = (para->maxRadius - para->minRadius)/(para->nRadius);
            minX = 2.0 * para->minRadius - shift;
            maxX = 2.0 * para->maxRadius + shift;
            bar->setWidth((para->maxRadius - para->minRadius)/para->nRadius);
            customPlot->xAxis->ticker()->setTickCount(5);
        }
        customPlot->xAxis->setRange(minX, maxX);
    }

    if (ui->radioButton_LinearYAxis->isChecked())
        customPlot->yAxis->setLabel("Concentration (Number Density: Ns)");
    if (ui->radioButton_LogYAxis->isChecked())
        customPlot->yAxis->setLabel("Log (Concentration (Number Density: Ns) )");
    if (ui->radioButton_LogXAxis->isChecked())
    {
        customPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);
        customPlot->xAxis->setLabel("Diameter (μm) {Axis in Log Scale}");
    }
    if (ui->radioButton_LinearXAxis->isChecked())
    {
        customPlot->xAxis->setScaleType(QCPAxis::stLinear);
        customPlot->xAxis->setLabel("Diameter (μm)");
    }
    customPlot->replot();
}

//Plot Phase function linear plot
void PlotData::PlotPhaseFunctionLinear(Ui_MainWindow *ui, QVector<double> x, QVector<double> yPara,
                                       QVector<double> yPerp, QVector<double> yAve)
{
    Utilities util;
    auto customPlot = ui->customPlot_PhaseFunctionLinear;
    RemoveGraphs(customPlot);
    RemoveLegends(customPlot);

    //Set Min and Max values for the y-axis
    double minY = util.FindMinMax(yPara, yPerp, yAve, true);   //true:Min
    double maxY = util.FindMinMax(yPara, yPerp, yAve, false);  //false: Max
    maxY = maxY + 0.02*(maxY-minY);
    minY = minY - 0.02*(maxY-minY);
    customPlot->yAxis->setRange(minY, maxY);    

    //Add new graph
    PlotSingleGraph(customPlot, x, yAve, Qt::red, "Ave.", 0, 2);
    PlotSingleGraph(customPlot, x, yPara, Qt::blue, "Para.", 1, 2);
    PlotSingleGraph(customPlot, x, yPerp, Qt::darkGreen, "Perp.", 2, 2);

    if (ui->radioButton_LinearYAxis->isChecked())
        customPlot->yAxis->setLabel("Magnitude");
    if (ui->radioButton_LogYAxis->isChecked())
        customPlot->yAxis->setLabel("Log (Magnitude)");

    //Legend    
    customPlot->legend->setBrush(QBrush(QColor(255, 255, 255, 0)));
    customPlot->legend->setBorderPen(QPen(Qt::NoPen));
    customPlot->legend->setVisible(true);
    customPlot->axisRect()->insetLayout()->setInsetAlignment(0,Qt::AlignTop|Qt::AlignRight);
    customPlot->replot();
}

//Plot phase function polar plot
void PlotData::PlotPhaseFunctionPolar(Ui_MainWindow *ui, QVector<double> theta, QVector<double> yPara,
                                      QVector<double> yPerp, QVector<double> yAve, int totalSize)
{
    double radialMin = 0;
    if (ui->radioButton_PhaseLog->isChecked())
        radialMin = log10(mPolarMinRadius);

    QVector<double> phaseData1 = yAve;  //Default
    QVector<double> phaseData2 = yPara;
    QVector<double> phaseData3 = yPerp;

    if (ui->radioButton_PhasePara->isChecked())
        phaseData1 = yPara;
    if (ui->radioButton_PhasePerp->isChecked())
        phaseData1 = yPerp;

    auto customPlot = ui->customPlot_PhaseFunctionPolar;

    if (ui->radioButton_PhaseAll->isChecked())
    {
        QVector<double> x1, y1;
        QVector<double> x2, y2;
        QVector<double> x3, y3;
        for (int i = 0; i < totalSize; i++)
        {
            double cost = cos(theta[i]);
            double sint = sin(theta[i]);

            double rho1 = phaseData1[i] - radialMin;
            x1.append(rho1 * cost);
            y1.append(rho1 * sint);

            double rho2 = phaseData2[i] - radialMin;
            x2.append(rho2 * cost);
            y2.append(rho2 * sint);

            double rho3 = phaseData3[i] - radialMin;
            x3.append(rho3 * cost);
            y3.append(rho3 * sint);
        }
        PlotSingleCurve(customPlot, x1, y1, Qt::red, "Ave.", totalSize);
        PlotSingleCurve(customPlot, x2, y2, Qt::blue, "Para.", totalSize);
        PlotSingleCurve(customPlot, x3, y3, Qt::darkGreen, "Perp.", totalSize);

        // Add legend
        customPlot->legend->setVisible(true);
        customPlot->legend->setBrush(QBrush(QColor(255, 255, 255, 0)));
        customPlot->legend->setBorderPen(QPen(Qt::NoPen));
        customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignBottom);
    }
    else
    {
        QVector<double> x1, y1;
        for (int i = 0; i < totalSize; i++)
        {
            double cost = cos(theta[i]);
            double sint = sin(theta[i]);
            double rho1 = phaseData1[i] - radialMin;
            x1.append(rho1 * cost);
            y1.append(rho1 * sint);
        }
        PlotSingleCurve(customPlot, x1, y1, Qt::red, "", totalSize);
    }
    customPlot->replot();
}

//Plot S1S2 plots
void PlotData::PlotS1S2(Ui_MainWindow *ui, Parameters *para, QVector<double> x, QVector<double> yS1, QVector<double> yS2)
{
    auto customPlot = ui->customPlot_S1S2;

    //remove previous graphs (if any)
    RemoveGraphs(customPlot);
    RemoveLegends(customPlot);

    //Add new graph
    if (ui->radioButton_S1->isChecked())
    {
        PlotSingleGraph(customPlot, x, yS1, Qt::red, "S1", 0, 2);
        customPlot->yAxis->setLabel("S1");
        if (ui->radioButton_LogYAxis->isChecked())
            customPlot->yAxis->setLabel("Log (S1)");
        customPlot->legend->setVisible(false);
    }
    if (ui->radioButton_S2->isChecked())
    {
        PlotSingleGraph(customPlot, x, yS2, Qt::red, "S2", 0, 2);
        customPlot->yAxis->setLabel("S2");
        if (ui->radioButton_LogYAxis->isChecked())
            customPlot->yAxis->setLabel("Log (S2)");
        customPlot->legend->setVisible(false);
    }
    if (ui->radioButton_S1S2->isChecked())
    {
        PlotSingleGraph(customPlot, x, yS1, Qt::red, "S1", 0, 2);
        PlotSingleGraph(customPlot, x, yS2, Qt::blue, "S2", 1, 2);
        customPlot->yAxis->setLabel("S1,S2");
        if (ui->radioButton_LogYAxis->isChecked())
            customPlot->yAxis->setLabel("Log (S1,S2)");
        customPlot->legend->setBrush(QBrush(QColor(255, 255, 255, 0)));
        customPlot->legend->setBorderPen(QPen(Qt::NoPen));
        customPlot->legend->setVisible(true);
    }
    customPlot->rescaleAxes(true);
    customPlot->replot();

    if (para->nWavel==1)
        ui->slider_WL_S1S2->setDisabled(true);
    else
        ui->slider_WL_S1S2->setDisabled(false);
}

//Plot Musp curve for power law
void PlotData::PlotMuspCurveForPowerLaw(Ui_MainWindow *ui, QVector<double> x, QVector<double> yMusp)
{
    auto customPlot = ui->customPlot_MuspPowerLaw;

    // remove previous graphs (if any)
    RemoveGraphs(customPlot);

    if (ui->radioButton_LinearYAxis->isChecked())
        customPlot->yAxis->setLabel("μs' (mmˉˡ)");
    if (ui->radioButton_LogYAxis->isChecked())
        customPlot->yAxis->setLabel("Log (μs' (mmˉˡ) )");

    //Add new graph
    PlotSingleGraph(customPlot, x, yMusp, Qt::red, "Musp", 0, 4);
    customPlot->graph(0)->rescaleAxes();
    customPlot->replot();
}

//Plot power law fit curve
void PlotData::PlotMuspPowerLaw(Ui_MainWindow *ui, QVector<double> x, QVector<double> yMusp, QVector<double> yFit)
{
    auto customPlot = ui->customPlot_MuspPowerLaw;

    // remove previous graphs (if any)
    RemoveGraphs(customPlot);    

    if (ui->radioButton_LinearYAxis->isChecked())
        customPlot->yAxis->setLabel("μs' (mmˉˡ)");
    if (ui->radioButton_LogYAxis->isChecked())
        customPlot->yAxis->setLabel("Log (μs' (mmˉˡ))");

    //Main plot
    PlotSingleGraph(customPlot, x, yMusp, Qt::red, "μs' Data", 0, 4);

    //Power law fit plot
    PlotSingleGraph(customPlot, x, yFit, Qt::blue, "Best Fit", 1, 4);

    customPlot->legend->setBrush(QBrush(QColor(255, 255, 255, 0)));
    customPlot->legend->setBorderPen(QPen(Qt::NoPen));
    customPlot->legend->setVisible(true);
    customPlot->rescaleAxes(true);
    customPlot->replot();
}

//Plot Scattering Cross Section plot
void PlotData::PlotScatteringCrossSection(Ui_MainWindow *ui, QVector<double> x, QVector<double> yCsca)
{
    auto customPlot = ui->customPlot_Csca;

    // remove previous graphs (if any)
    RemoveGraphs(customPlot);

    if (ui->radioButton_LinearYAxis->isChecked())
        customPlot->yAxis->setLabel("Scattering Cross Section (μm²)");
    if (ui->radioButton_LogYAxis->isChecked())
        customPlot->yAxis->setLabel("Log (Scattering Cross Section (μm²) )");

    //Add new graph
    PlotSingleGraph(customPlot, x, yCsca, Qt::red, "Csca", 0, 4);
    customPlot->graph(0)->rescaleAxes();
    customPlot->replot();
}

//Plot Extinction Cross Section plot
void PlotData::PlotExtinctionCrossSection(Ui_MainWindow *ui, QVector<double> x, QVector<double> yCext)
{
    auto customPlot = ui->customPlot_Cext;

    // remove previous graphs (if any)
    RemoveGraphs(customPlot);

    if (ui->radioButton_LinearYAxis->isChecked())
        customPlot->yAxis->setLabel("Extinction Cross Section (μm²)");
    if (ui->radioButton_LogYAxis->isChecked())
        customPlot->yAxis->setLabel("Log (Extinction Cross Section (μm²) )");

    //Add new graph
    PlotSingleGraph(customPlot, x, yCext, Qt::red, "Cext", 0, 4);
    customPlot->graph(0)->rescaleAxes();
    customPlot->replot();
}

//Plot Backscattering Cross Section plot
void PlotData::PlotBackscatteringCrossSection(Ui_MainWindow *ui, QVector<double> x, QVector<double> yCback)
{
    auto customPlot = ui->customPlot_Cback;

    // remove previous graphs (if any)
    RemoveGraphs(customPlot);

    if (ui->radioButton_LinearYAxis->isChecked())
        customPlot->yAxis->setLabel("Backscattering Cross Section (μm²)");
    if (ui->radioButton_LogYAxis->isChecked())
        customPlot->yAxis->setLabel("Log (Backscattering Cross Section (μm²) )");

    //Add new graph
    PlotSingleGraph(customPlot, x, yCback, Qt::red, "Cback", 0, 4);
    customPlot->graph(0)->rescaleAxes();
    customPlot->replot();
}

//Plot Size Parameter
void PlotData::PlotSizeParameter(Ui_MainWindow *ui, QVector<double> x, QVector<double> ySizePara)
{
    auto customPlot = ui->customPlot_SizePara;

    // remove previous graphs (if any)
    RemoveGraphs(customPlot);

    if (ui->radioButton_LinearYAxis->isChecked())
        customPlot->yAxis->setLabel("Size Parameter");
    if (ui->radioButton_LogYAxis->isChecked())
        customPlot->yAxis->setLabel("Log (Size Parameter)");

    //Add new graph
    PlotSingleGraph(customPlot, x, ySizePara, Qt::red, "Size Para", 0, 4);
    customPlot->graph(0)->rescaleAxes();
    customPlot->replot();
}

//Plot Mus plot
void PlotData::PlotMus(Ui_MainWindow *ui,QVector<double> x, QVector<double> yMus)
{
    auto customPlot = ui->customPlot_Mus;

    // remove previous graphs (if any)
    RemoveGraphs(customPlot);

    if (ui->radioButton_LinearYAxis->isChecked())
        customPlot->yAxis->setLabel("μs (mmˉˡ)");
    if (ui->radioButton_LogYAxis->isChecked())
        customPlot->yAxis->setLabel("Log (μs (mmˉˡ) )");

    //Add new graph
    PlotSingleGraph(customPlot, x, yMus, Qt::red, "Mus", 0, 4);
    customPlot->graph(0)->rescaleAxes();
    customPlot->replot();
}

//Plot g (average phase function) plot
void PlotData::PlotG(Ui_MainWindow *ui, QVector<double> x, QVector<double> yG)
{
    auto customPlot = ui->customPlot_G;

    // remove previous graphs (if any)
    RemoveGraphs(customPlot);

    if (ui->radioButton_LinearYAxis->isChecked())
        customPlot->yAxis->setLabel("g (Average Cosine of phase function)");
    if (ui->radioButton_LogYAxis->isChecked())
        customPlot->yAxis->setLabel("Log (g (Average Cosine of phase function))");

    //Add new graph
    PlotSingleGraph(customPlot, x, yG, Qt::red, "g", 0, 4);
    customPlot->graph(0)->rescaleAxes();
    customPlot->replot();
}

//Plot Musp plot
void PlotData::PlotMusp(Ui_MainWindow *ui, QVector<double> x, QVector<double> yMusp)
{
    auto customPlot = ui->customPlot_Musp;

    // remove previous graphs (if any)
    RemoveGraphs(customPlot);

    if (ui->radioButton_LinearYAxis->isChecked())
        customPlot->yAxis->setLabel("μs' (mmˉˡ)");
    if (ui->radioButton_LogYAxis->isChecked())
        customPlot->yAxis->setLabel("Log (μs' (mmˉˡ) )");

    //Add new graph
    PlotSingleGraph(customPlot, x, yMusp, Qt::red, "Musp", 0, 4);
    customPlot->graph(0)->rescaleAxes();
    customPlot->replot();
}

//Plot Forward /Backward scattering percentage plot
void PlotData::PlotForwardBackward(Ui_MainWindow *ui, QVector<double> x, QVector<double> yF, QVector<double> yB, bool legendFlag)
{
    auto customPlot = ui->customPlot_FB;

    // remove previous graphs and items (if any)
    RemoveGraphs(customPlot);
    RemoveItems(customPlot);
    RemoveLegends(customPlot);

    //Forward and backward Plot
    PlotSingleGraph(customPlot, x, yF, Qt::red, "Forward Scat. %", 0, 4);
    PlotSingleGraph(customPlot, x, yB, Qt::blue, "Backward Scat. %", 1, 4);

    if (ui->radioButton_LinearYAxis->isChecked())
    {
        customPlot->yAxis->setLabel("Forward & Backward Scattering %");
        customPlot->rescaleAxes(true);
        customPlot->yAxis->setRange(0,115);
    }
    if (ui->radioButton_LogYAxis->isChecked())
    {
        customPlot->yAxis->setLabel("Log (Forward & Backward Scattering %)");
        customPlot->rescaleAxes(true);
        customPlot->yAxis->setRange(-3, 4);
    }
    customPlot->legend->setVisible(true);
    customPlot->legend->setBrush(QBrush(QColor(255, 255, 255, 0)));
    customPlot->legend->setBorderPen(QPen(Qt::NoPen));

    if (legendFlag)
        customPlot->axisRect()->insetLayout()->setInsetAlignment(0,Qt::AlignCenter|Qt::AlignRight);
    else
        customPlot->axisRect()->insetLayout()->setInsetAlignment(0,Qt::AlignTop|Qt::AlignRight);
    customPlot->replot();
}

//Plot single graph
void PlotData::PlotSingleGraph(QCustomPlot* customPlot, const QVector<double> &x, const QVector<double> &y,
                     QColor color, const QString &name, int graphIndex, int sizeCircle)
{
    customPlot->addGraph();
    customPlot->graph(graphIndex)->setData(x, y);
    customPlot->graph(graphIndex)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, sizeCircle));
    customPlot->graph(graphIndex)->setPen(QPen(color, 1));
    customPlot->graph(graphIndex)->setName(name);
}

//Plot Single Curve
void PlotData::PlotSingleCurve(QCustomPlot* customPlot, const QVector<double>& xData, const QVector<double>& yData,
                               const QColor& color, const QString& name, int totalSize)
{
    if (xData.size() == totalSize && yData.size() == totalSize)
    {
        QCPCurve* curve = new QCPCurve(customPlot->xAxis, customPlot->yAxis);
        curve->setData(xData, yData);
        curve->setPen(QPen(color, 2));
        curve->setName(name);
    }
}

//Initial setup plot
void PlotData::InitialSetupPlot(QCustomPlot *customPlot, const QString &xLabel, const QString &yLabel, double minX, double maxX)
{
    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    customPlot->xAxis->setLabel(xLabel);
    customPlot->xAxis->setRange(minX, maxX);
    customPlot->yAxis->setLabel(yLabel);
    customPlot->yAxis->setRange(0, 5);
    customPlot->clearGraphs();
    customPlot->replot();
}

//Draw axes of polar plot
void PlotData::DrawPolarPlotGrid(QCustomPlot *customPlot, bool flagLinearLog)
{    
    RemoveLegends(customPlot);
    RemovePlotables(customPlot);
    RemoveEllipseLineTextItems(customPlot);
    customPlot->replot();

    //Add Radial and circular grids
    CreateRadialGrid(customPlot, flagLinearLog);
    CreateCircularGrid(customPlot, flagLinearLog);
    customPlot->replot();
}

//Create Circular Grid
void PlotData::CreateCircularGrid(QCustomPlot *customPlot, bool flagLinearLog)
{
    QPen majorGridPen(Qt::black, 1.0, Qt::SolidLine);
    QPen minorGridPen(Qt::gray, 0.5, Qt::SolidLine);
    QPen boundaryGridPen(Qt::black, 0.5, Qt::SolidLine);

    Utilities util;

    if (flagLinearLog)  // Linear spacing
    {
        // Start with a large number of circles and reduce until it's a "nice" number
        int numCircles = 10;
        double majorTickStep = util.NiceStep(mPolarMaxRadius, numCircles);
        numCircles = mPolarMaxRadius / majorTickStep;
        while (numCircles > 5)
        {
            numCircles--;
            majorTickStep = util.NiceStep(mPolarMaxRadius, numCircles);
        }

        double currentRadius = 0;
        while (currentRadius <= mPolarMaxRadius)
        {
            // Draw the major circle
            QCPItemEllipse *majorCircle = new QCPItemEllipse(customPlot);
            majorCircle->topLeft->setCoords(-currentRadius, currentRadius);
            majorCircle->bottomRight->setCoords(currentRadius, -currentRadius);
            majorCircle->setPen(majorGridPen);

            // Add label for the major circle
            QCPItemText *label = new QCPItemText(customPlot);
            label->setText(QString::number(currentRadius, 'g', 3));
            label->position->setCoords(0, currentRadius);
            label->setFont(QFont(customPlot->font().family(), 9));
            label->setPen(QPen(Qt::NoPen));
            label->setBrush(QBrush(Qt::white));

            // Draw minor circles between the current major and the next one
            double minorTickStep = majorTickStep / 5.0; // 5 minor circles between each major
            for (int i = 1; i < 5; ++i) {
                double minorRadius = currentRadius + minorTickStep * i;
                if (minorRadius > mPolarMaxRadius)
                    break; // Don't draw minor circles beyond the max radius

                QCPItemEllipse *minorCircle = new QCPItemEllipse(customPlot);
                minorCircle->topLeft->setCoords(-minorRadius, minorRadius);
                minorCircle->bottomRight->setCoords(minorRadius, -minorRadius);
                minorCircle->setPen(minorGridPen);
            }

            currentRadius += majorTickStep;
        }
        //Add boundary circle
        if (mPolarMaxRadius - floor(mPolarMaxRadius)>0)
        {
            QCPItemEllipse *boundaryCircle = new QCPItemEllipse(customPlot);
            boundaryCircle->topLeft->setCoords(-mPolarMaxRadius, mPolarMaxRadius);
            boundaryCircle->bottomRight->setCoords(mPolarMaxRadius, -mPolarMaxRadius);
            boundaryCircle->setPen(boundaryGridPen);
        }
    }
    else //Log spacing    
    {
        double logMin = log10(mPolarMinRadius);
        double logMax = log10(mPolarMaxRadius);
        double radialRange =logMax - logMin;

        // Determine the step for major circles and the minor radius multipliers
        int logStep = 1;
        QVector<int> minorMultipliers;

        if (radialRange > 10)
        {
            logStep = 3;
            minorMultipliers << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9;
            minorMultipliers << 10 << 20 << 30 << 40 << 50 << 60 << 70 << 80 << 90;
            minorMultipliers << 100 << 200 << 300 << 400 << 500 << 600 << 700 << 800 << 900;
        }
        else if (radialRange > 5)
        {
            logStep = 2;
            minorMultipliers << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9;
            minorMultipliers << 10 << 20 << 30 << 40 << 50 << 60 << 70 << 80 << 90;
        }
        else
        {
            logStep = 1;
            minorMultipliers << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9;
        }

        int startLog = floor(logMin);

        for (int i = startLog; i < logMax; i += logStep)
        {
            double majorRadius = pow(10, i);
            // Only draw if the radius is within the plot's range
            if (majorRadius >= mPolarMinRadius - 1e-9 && majorRadius <= mPolarMaxRadius + 1e-9)
            {
                QCPItemEllipse *circle = new QCPItemEllipse(customPlot);
                circle->topLeft->setCoords((-i + logMin), (i - logMin));
                circle->bottomRight->setCoords((i - logMin), (-i + logMin));
                circle->setPen(majorGridPen);

                QCPItemText *label = new QCPItemText(customPlot);
                label->setText(QString::number(majorRadius, 'g', 3));
                label->position->setCoords(0, (i - logMin));
                label->setFont(QFont(customPlot->font().family(), 9));
                label->setPen(QPen(Qt::NoPen));
                label->setBrush(QBrush(Qt::white));
            }

            // Draw minor circles between the current and the next major circle
            double nextMajorRadius = pow(10, i + logStep);
            if (nextMajorRadius > mPolarMaxRadius)
                nextMajorRadius = mPolarMaxRadius;

            // Use a const reference to avoid C++11 range-loop might detach Qt container warning
            const QVector<int>& multipliers = minorMultipliers;
            for (int factor : multipliers)
            {
                double minorRadius = pow(10, i) * factor;

                if (minorRadius > mPolarMaxRadius)
                    break;

                if (minorRadius >= mPolarMinRadius - 1e-9 && minorRadius <= nextMajorRadius + 1e-9)
                {
                    QCPItemEllipse *minorCircle = new QCPItemEllipse(customPlot);
                    minorCircle->topLeft->setCoords((-log10(minorRadius) + logMin), (log10(minorRadius) - logMin));
                    minorCircle->bottomRight->setCoords((log10(minorRadius) - logMin), (-log10(minorRadius) + logMin));
                    minorCircle->setPen(minorGridPen);
                }
            }
        }
        // Add boundary circle
        if (radialRange - floor(radialRange) > 0)
        {
            QCPItemEllipse *boundaryCircle = new QCPItemEllipse(customPlot);
            boundaryCircle->topLeft->setCoords(- radialRange, radialRange);
            boundaryCircle->bottomRight->setCoords(radialRange, -radialRange);
            boundaryCircle->setPen(boundaryGridPen);
        }
    }
    customPlot->replot();
}

//Create Radial Grid
void PlotData::CreateRadialGrid(QCustomPlot *customPlot, bool flagLinearLog)
{
    double radialMax;
    if (flagLinearLog)  // Linear spacing
        radialMax = mPolarMaxRadius;
    else //Log spacing
    {
        double logMin = log10(mPolarMinRadius);
        double logMax = log10(mPolarMaxRadius);
        radialMax = logMax - logMin;
    }

    QPen majorGridPen(Qt::black, 1.0, Qt::SolidLine);
    QPen minorGridPen(Qt::gray, 0.5, Qt::SolidLine);
    double labelRadius = radialMax * 1.17;
    int numLines = 36;
    double angleStep = 10;

    for (int i = 0; i < numLines; ++i)
    {
        QCPItemLine *radialLine = new QCPItemLine(customPlot);
        double angleRad = angleStep * i * M_PI / 180.0;

        radialLine->start->setCoords(0, 0);
        radialLine->end->setCoords(radialMax * cos(angleRad), radialMax * sin(angleRad));

        if (i % 3 == 0)  // Major line and label
        {
            radialLine->setPen(majorGridPen);
            QCPItemText *angleText = new QCPItemText(customPlot);
            angleText->position->setCoords(labelRadius * cos(angleRad), labelRadius * sin(angleRad));
            angleText->setText(QString::number(angleStep * i) + "°");
            angleText->setPen(QPen(Qt::NoPen));
            angleText->setBrush(QBrush(Qt::NoBrush));
        } else  // Minor line
            radialLine->setPen(minorGridPen);
    }
    customPlot->xAxis->setRange(-1.3*radialMax, 1.3*radialMax);
    customPlot->yAxis->setRange(-1.3*radialMax, 1.3*radialMax);
    HideCartesianAxes(customPlot);
}

//Hide Cartesian axes for polar plot
void PlotData::HideCartesianAxes(QCustomPlot *customPlot)
{
    customPlot->xAxis->grid()->setVisible(false);
    customPlot->yAxis->grid()->setVisible(false);
    customPlot->xAxis->setVisible(false);
    customPlot->yAxis->setVisible(false);
    customPlot->xAxis->setTickLabels(false);
    customPlot->yAxis->setTickLabels(false);
    customPlot->xAxis->setTicks(false);
    customPlot->yAxis->setTicks(false);
}

//Remove QCustomPlot graphs
void PlotData::RemoveGraphs(QCustomPlot *customPlot)
{
    if (customPlot->graphCount() > 0)
    {
        for (int i = customPlot->graphCount() - 1; i >= 0; --i)
            customPlot->removeGraph(i);
    }
}

void PlotData::RemoveLegends(QCustomPlot *customPlot)
{
    if (customPlot->legend->itemCount() >0)
    {
        for (int i = customPlot->itemCount() - 1; i >= 0; --i)
            customPlot->legend->removeItem(i);
    }
}

//Remove QCustomPlot items
void PlotData::RemoveItems(QCustomPlot *customPlot)
{
    if (customPlot->itemCount() >0)
    {
        for (int i = customPlot->itemCount() - 1; i >= 0; --i)
        {
            QCPAbstractItem *item = customPlot->item(i);
            customPlot->removeItem(item);
        }
    }
}

//Remove Ellipse and Lines
void PlotData::RemoveEllipseLineTextItems(QCustomPlot *customPlot)
{
    if (customPlot->itemCount() >0)
    {
        for (int i = customPlot->itemCount() - 1; i >= 0; --i)
        {
            QCPAbstractItem *item = customPlot->item(i);
            if (qobject_cast<QCPItemEllipse*>(item) || qobject_cast<QCPItemLine*>(item) ||
                qobject_cast<QCPItemText*>(item))
                customPlot->removeItem(item);
        }
    }
}

//Remove QCustomPlot plottables
void PlotData::RemovePlotables(QCustomPlot *customPlot)
{
    if (customPlot->plottableCount() >0)
    {
        for (int i = customPlot->plottableCount() - 1; i >= 0; --i)
            customPlot->removePlottable(i);
    }

}

//Rearrange Phase Function Data
void PlotData::RearrangePhaseFunctionData(Parameters *para, QVector<double> &theta,
                                          QVector<double> &phaseFuncPara,
                                          QVector<double> &phaseFuncPerp,
                                          QVector<double> &phaseFuncAve,
                                          int indexWL, bool flagLinearLog,
                                          bool flagThetaNegPosOrPos)
{
    double tiny = 1e-100;  //add a small number before log calculation to avoid NaN
    int nTheta = static_cast<int>(para->nTheta);

    //Assign "rho(r)" and "theta" values in polar plot
    for (int i=0; i<nTheta; i++)
    {
        double theta_val1 = -180 + (180.0 * i /(nTheta-1));
        double theta_val2 = 180 - (180.0 * i /(nTheta-1));

        //Theta is from -180 to 180 in degrees
        theta[i] = theta_val1;
        theta[2 * nTheta -2  -i] =  theta_val2;

        //Theta is from 0 to 360 in radians
        if (!flagThetaNegPosOrPos)
        {
            if (theta_val1 < 0)
                theta_val1 += 360;
            if (theta_val2 < 0)
                theta_val2 += 360;

            theta[i] = M_PI * theta_val1 /180;
            theta[2 * nTheta -2  -i] =  M_PI * theta_val2 /180;
        }

        double phaseFunctionPara = para->phaseFunctionPara[indexWL][i];
        double phaseFunctionPerp = para->phaseFunctionPerp[indexWL][i];
        double phaseFunctionAve = para->phaseFunctionAve[indexWL][i];

        int forwardIndex = nTheta - 1 + i;
        int backwardIndex = nTheta - 1 - i;

        if (flagLinearLog)  // Linear scale
        {
            phaseFuncPara[backwardIndex] = phaseFunctionPara;
            phaseFuncPara[forwardIndex] = phaseFunctionPara;

            phaseFuncPerp[backwardIndex] = phaseFunctionPerp;
            phaseFuncPerp[forwardIndex] = phaseFunctionPerp;

            phaseFuncAve[backwardIndex] = phaseFunctionAve;
            phaseFuncAve[forwardIndex] = phaseFunctionAve;
        }
        else  //log scale
        {
            phaseFuncPara[backwardIndex] = log10(phaseFunctionPara+tiny);
            phaseFuncPara[forwardIndex] = log10(phaseFunctionPara+tiny);

            phaseFuncPerp[backwardIndex] = log10(phaseFunctionPerp+tiny);
            phaseFuncPerp[forwardIndex] = log10(phaseFunctionPerp+tiny);

            phaseFuncAve[backwardIndex] = log10(phaseFunctionAve+tiny);
            phaseFuncAve[forwardIndex] = log10(phaseFunctionAve+tiny);
        }
    }
}



