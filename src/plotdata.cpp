//**********************************************************************
//** All QCustomPlot asignments and plotting actions are listed here
//**********************************************************************

#include "plotdata.h"

const QwtInterval radialInterval( 0.0, 10.0 );
const QwtInterval azimuthInterval( 0.0, 360.0 );

PlotData::PlotData(void)
{
}

PlotData::~PlotData(void)
{
}

//Initialize distribution plot
void PlotData::InitializeDistributionPlot(Ui_MainWindow *ui)
{
    double minX, maxX;
    minX = 0.95 * (ui->lineEdit_diameter->text().toDouble());
    maxX = 1.05 * (ui->lineEdit_diameter->text().toDouble());

    ui->customPlot_distribution->xAxis->setLabel("Diameter (μm)");
    ui->customPlot_distribution->xAxis->setRange(minX, maxX);
    ui->customPlot_distribution->yAxis->setLabel("Ns  (/vol. of 1mm³)");
    ui->customPlot_distribution->clearPlottables();
    ui->customPlot_distribution->replot();
}

//Initialize polar plot
void PlotData::InitializePhaseFunctionPolarPlot(Ui_MainWindow *ui, parameters *para )
{
//    SetSliderWL(ui);

    //Polar plot: phase function
    ui->qwtpolarplot_pFunctionPolar->setAutoReplot( false );
    ui->qwtpolarplot_pFunctionPolar->setPlotBackground( Qt::white );
    ui->qwtpolarplot_pFunctionPolar->setScale( QwtPolar::Azimuth,0, 360, 360 / 12 );
    ui->qwtpolarplot_pFunctionPolar->setScale( QwtPolar::Radius, 0, 1);
    mGrid = new QwtPolarGrid();
    mGrid->setPen( QPen( Qt::black ) );
    for ( int scaleId = 0; scaleId < QwtPolar::ScaleCount; scaleId++ )
    {
        mGrid->showGrid( scaleId );
        mGrid->showMinorGrid( scaleId );
        QPen minorPen( Qt::gray );
        #if 0
            minorPen.setStyle( Qt::DotLine );
        #endif
        mGrid->setMinorGridPen( scaleId, minorPen );
    }
    mGrid->setAxisPen( QwtPolar::AxisAzimuth, QPen( Qt::black ) );
    mGrid->showAxis( QwtPolar::AxisAzimuth, true );
    mGrid->showAxis( QwtPolar::AxisLeft, false );
    mGrid->showAxis( QwtPolar::AxisRight, false );
    mGrid->showAxis( QwtPolar::AxisTop, true );
    mGrid->showAxis( QwtPolar::AxisBottom, false );
    mGrid->showGrid( QwtPolar::Azimuth, true );
    mGrid->showGrid( QwtPolar::Radius, true );
    mGrid->attach( ui->qwtpolarplot_pFunctionPolar );
    para->polarCurve->detach();
    ui->qwtpolarplot_pFunctionPolar->replot();
}

//Initialize distribution plot
void PlotData::InitializePhaseFunctionLinearPlot(Ui_MainWindow *ui)
{
    ui->customPlot_pFunctionLinear->xAxis->setLabel("Angle (deg.)");
    ui->customPlot_pFunctionLinear->xAxis->setRange(-180,180);
    QVector<double> ticks;
    QVector<QString> labels;
    ticks << -180 << -135 << -90 << -45 << 0 << 45 << 90 << 135 << 180;
    labels << "-180" << "135" << "-90" << "-45" << "0" << "45" << "90" << "135" << "180";
    QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
    textTicker->addTicks(ticks, labels);
    ui->customPlot_pFunctionLinear->xAxis->setTicker(textTicker);
    ui->customPlot_pFunctionLinear->yAxis->setLabel("Magnitude");
    ui->customPlot_pFunctionLinear->clearGraphs();
    ui->customPlot_pFunctionLinear->replot();
    ui->customPlot_pFunctionLinear->addGraph(); //dummy
    ui->customPlot_pFunctionLinear->addGraph(); //dummy
    ui->customPlot_pFunctionLinear->addGraph(); //dummy
    ui->customPlot_pFunctionLinear->legend->setVisible(false);
}

//Initialize all other plots
void PlotData::InitializeAllOtherPlots(Ui_MainWindow *ui)
{
    double minX = ui->lineEdit_startWavel->text().toDouble();
    double maxX = ui->lineEdit_endWavel->text().toDouble();

    //Plot S1S2
    ui->customPlot_s1s2->xAxis->setLabel("Angle (deg.)");
    ui->customPlot_s1s2->xAxis->setRange(0,180);
    QVector<double> ticks;
    QVector<QString> labels;
    ticks << 0 << 45 << 90 << 135 << 180;
    labels << "0" << "45" << "90" << "135" << "180";
    QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
    textTicker->addTicks(ticks, labels);
    ui->customPlot_s1s2->xAxis->setTicker(textTicker);
    ui->customPlot_s1s2->yAxis->setLabel("S1");
    ui->customPlot_s1s2->clearGraphs();
    ui->customPlot_s1s2->replot();
    ui->customPlot_s1s2->addGraph();   //dummy

    //Plot Scattering Cross Section
    ui->customPlot_csca->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_csca->xAxis->setRange(minX, maxX);
    ui->customPlot_csca->yAxis->setLabel("Scattering Cross Section (μm²)");
    ui->customPlot_csca->clearGraphs();
    ui->customPlot_csca->replot();
    ui->customPlot_csca->addGraph();   //dummy

    //Plot Extinction Cross Section
    ui->customPlot_cext->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_cext->xAxis->setRange(minX, maxX);
    ui->customPlot_cext->yAxis->setLabel("Extinction Cross Section (μm²)");
    ui->customPlot_cext->clearGraphs();
    ui->customPlot_cext->replot();
    ui->customPlot_cext->addGraph();   //dummy

    //Plot Backscattering Cross Section
    ui->customPlot_cback->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_cback->xAxis->setRange(minX, maxX);
    ui->customPlot_cback->yAxis->setLabel("Backscattering Cross Section (μm²)");
    ui->customPlot_cback->clearGraphs();
    ui->customPlot_cback->replot();
    ui->customPlot_cback->addGraph();   //dummy

    //Plot Size Parameter
    ui->customPlot_sizePara->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_sizePara->xAxis->setRange(minX, maxX);
    ui->customPlot_sizePara->yAxis->setLabel("Size Parameter");
    ui->customPlot_sizePara->clearGraphs();
    ui->customPlot_sizePara->replot();
    ui->customPlot_sizePara->addGraph();   //dummy

    //Plot Mus
    ui->customPlot_mus->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_mus->xAxis->setRange(minX, maxX);
    ui->customPlot_mus->yAxis->setLabel("μs (mmˉˡ)");
    ui->customPlot_mus->clearGraphs();
    ui->customPlot_mus->replot();
    ui->customPlot_mus->addGraph();   //dummy

    //Plot G
    ui->customPlot_g->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_g->xAxis->setRange(minX, maxX);
    ui->customPlot_g->yAxis->setLabel("g (Average Cosine of phase function)");
    ui->customPlot_g->clearGraphs();
    ui->customPlot_g->replot();
    ui->customPlot_g->addGraph();   //dummy

    //Plot FB
    ui->customPlot_forwardBackward->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_forwardBackward->xAxis->setRange(minX, maxX);
    ui->customPlot_forwardBackward->yAxis->setLabel("Forward & Backward Scattering %");
    ui->customPlot_forwardBackward->clearGraphs();
    ui->customPlot_forwardBackward->replot();
    ui->customPlot_forwardBackward->addGraph(); //dummy
    ui->customPlot_forwardBackward->addGraph(); //dummy
    ui->customPlot_forwardBackward->legend->setVisible(false);

    //Plot Musp
    ui->customPlot_musp->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_musp->xAxis->setRange(minX, maxX);
    ui->customPlot_musp->yAxis->setLabel("μs' (mmˉˡ)");
    ui->customPlot_musp->clearGraphs();
    ui->customPlot_musp->replot();
    ui->customPlot_musp->addGraph();  //dummy

    //Plot Musp for fitting
    ui->customPlot_muspPowerLaw->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_muspPowerLaw->xAxis->setRange(minX, maxX);
    ui->customPlot_muspPowerLaw->yAxis->setLabel("μs' (mmˉˡ)");
    ui->customPlot_muspPowerLaw->clearGraphs();
    ui->customPlot_muspPowerLaw->replot();
    ui->customPlot_muspPowerLaw->addGraph();  //dummy
    ui->customPlot_muspPowerLaw->addGraph();  //dummy
    ui->customPlot_muspPowerLaw->legend->setVisible(false);
}

//Clear plots
void PlotData::ClearPlots(Ui_MainWindow *ui, parameters *para)
{
    InitializeDistributionPlot(ui);
    InitializePhaseFunctionPolarPlot(ui, para);
    InitializePhaseFunctionLinearPlot(ui);    
    InitializeAllOtherPlots(ui);
}

//Assign values for distribution plot
void PlotData::AssignValuesDistributionPlot(Ui_MainWindow *ui, parameters* para)
{
    QVector<double> xDist(static_cast<int>(para->nRadius));
    QVector<double> yDist(static_cast<int>(para->nRadius));
    double tiny = 1e-100;  //add a small number before log calculation to avoid NaN
    double margin = (1.0 + ui->slider_concPercentChange->value()/200.0);

    double totNumDen = 0.0;
    for (int i=0; i<static_cast<int>(para->nRadius); i++)
    {
        xDist[i] = 2.0 * para->radArray[i];
        yDist[i] = para->numDensityArray[i]*margin;
        if (ui->radioButton_logYaxis->isChecked())
            yDist[i] = log10(yDist[i] + tiny);
        totNumDen += para->numDensityArray[i]*margin;
    }
    if (ui->radioButton_monoDisperse->isChecked())
        ui->label_currentTotNumDen->setVisible(false);
    if (ui->radioButton_polyDisperse->isChecked())
    {
        ui->label_currentTotNumDen->setVisible(true);
        ui->label_currentTotNumDen->setText("<font color=\"red\">Σ Ns: "+QString::number(totNumDen)+" in vol. of 1mm³");
    }
    PlotDistribution(ui, para, xDist,yDist);
}

//plot distribution plot
void PlotData::PlotDistribution(Ui_MainWindow *ui, parameters *para, QVector<double> x, QVector<double> yDist)
{
    double minX, maxX, minY, maxY;

    minY = 0;
    maxY = yDist[0];
    for(int i = 0; i<static_cast<int>(para->nRadius); i++)
    {
        if (maxY < yDist[i])
               maxY = yDist[i];
    }
    maxY *= 1.1;

    ui->customPlot_distribution->removeGraph(0);
    ui->customPlot_distribution->clearPlottables();
    ui->customPlot_distribution->replot();
    ui->customPlot_distribution->addGraph();
    ui->customPlot_distribution->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);

    ui->customPlot_distribution->graph()->setData(x, yDist);
    ui->customPlot_distribution->graph()->setPen( QPen( Qt::red, 1 ) );
    ui->customPlot_distribution->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 3));

    QCPBars *bar = new QCPBars(ui->customPlot_distribution->xAxis, ui->customPlot_distribution->yAxis);
    bar->setData(x, yDist);
    ui->customPlot_distribution->yAxis->setRange(minY, maxY);
    bar->setPen(Qt::NoPen);
    bar->setBrush(QColor(0,0,0,255));

    if (ui->radioButton_monoDisperse->isChecked())
    {
        minX = 0.92 * 2.0 * para->minRadius;
        maxX = 1.08 * 2.0 * para->maxRadius;
        bar->setWidth(0.05*(maxX -minX));
        ui->customPlot_distribution->xAxis->setRange(minX, maxX);
        ui->customPlot_distribution->xAxis->ticker()->setTickCount(3);
    }
    if (ui->radioButton_polyDisperse->isChecked())
    {
        double shift = (para->maxRadius - para->minRadius)/(para->nRadius);
        if (shift == 0.0)      //Array of spherical scatterers with similar radius
        {
            minX = para->minRadius;
            maxX = 3.0 * para->maxRadius;
            bar->setWidth(0.1*para->maxRadius);
        }
        else
        {
            minX = 2.0 * para->minRadius - shift;
            maxX = 2.0 * para->maxRadius + shift;
            bar->setWidth((para->maxRadius - para->minRadius)/para->nRadius);
            ui->customPlot_distribution->xAxis->ticker()->setTickCount(5);
        }
        ui->customPlot_distribution->xAxis->setRange(minX, maxX);
    }

    if (ui->radioButton_linearYaxis->isChecked())
        ui->customPlot_distribution->yAxis->setLabel("Concentration (Number Density: Ns)");
    if (ui->radioButton_logYaxis->isChecked())
        ui->customPlot_distribution->yAxis->setLabel("Log (Concentration (Number Density: Ns) )");
    if (ui->radioButton_logXaxis->isChecked())
    {
        ui->customPlot_distribution->xAxis->setScaleType(QCPAxis::stLogarithmic);
        ui->customPlot_distribution->xAxis->setLabel("Diameter (μm) {Axis in Log Scale}");
    }
    if (ui->radioButton_linearXaxis->isChecked())
    {
        ui->customPlot_distribution->xAxis->setScaleType(QCPAxis::stLinear);
        ui->customPlot_distribution->xAxis->setLabel("Diameter (μm)");
    }
    ui->customPlot_distribution->replot();
}

//Assign values for phase function polar plot
void PlotData::AssignValuesPhaseFunctionPolarPlot(Ui_MainWindow *ui, parameters *para)
{
    QVector<double> phaseFunction(static_cast<int>(2*para->nTheta-1));
    QVector<double> theta(static_cast<int>(2*para->nTheta-1));
    int indexWL = ui->slider_pFunctionPolar_wavel->value();
    ui->label_pFunctionPolar_curWavel->setText(QString::number(para->startWavel + indexWL*para->stepWavel));

    //Set minimum and maximum values
    para->maxPolarPtheta = 0;
    para->minPolarPtheta = 1e100;
    for (int i = 0; i < static_cast<int>(para->nWavel); i++)
    {
        for (int j = 0; j < static_cast<int>(para->nTheta); j++)
        {
            //Get Max value
            if (para->phaseFunctionAve[i][j]>para->maxPolarPtheta)
                para->maxPolarPtheta = para->phaseFunctionAve[i][j];
            if (para->phaseFunctionPara[i][j]>para->maxPolarPtheta)
                para->maxPolarPtheta = para->phaseFunctionPara[i][j];
            if (para->phaseFunctionPerp[i][j]>para->maxPolarPtheta)
                para->maxPolarPtheta = para->phaseFunctionPerp[i][j];
            //Get Min value
            if (para->phaseFunctionAve[i][j]<para->minPolarPtheta)
                para->minPolarPtheta = para->phaseFunctionAve[i][j];
            if (para->phaseFunctionPara[i][j]<para->minPolarPtheta)
                para->minPolarPtheta = para->phaseFunctionPara[i][j];
            if (para->phaseFunctionPerp[i][j]<para->minPolarPtheta)
                para->minPolarPtheta = para->phaseFunctionPerp[i][j];
        }
    }

    //Assign "rho(r)" and "theta" values in polar plot
    for (int i=0; i<static_cast<int>(para->nTheta); i++)
    {
        theta[i] = 180.0 * i /(para->nTheta-1);
        theta[static_cast<int>(2*para->nTheta)-2 -i] =  360.0 - (180.0 * i /(para->nTheta-1));
        if (ui->radioButton_pFunction_ave->isChecked())
        {
            phaseFunction[i] = para->phaseFunctionAve[indexWL][i];
            phaseFunction[static_cast<int>(2*para->nTheta)-2 -i] = para->phaseFunctionAve[indexWL][i];
        }
        if (ui->radioButton_pFunction_para->isChecked())
        {
            phaseFunction[i] = para->phaseFunctionPara[indexWL][i];
            phaseFunction[static_cast<int>(2*para->nTheta)-2 -i] = para->phaseFunctionPara[indexWL][i];
        }
        if (ui->radioButton_pFunction_perp->isChecked())
        {
            phaseFunction[i] = para->phaseFunctionPerp[indexWL][i];
            phaseFunction[static_cast<int>(2*para->nTheta)-2 -i] = para->phaseFunctionPerp[indexWL][i];
        }
    }
    PlotPhaseFunctionPolar(ui, para, theta, phaseFunction);
}

//Plot phase function polar plot
void PlotData::PlotPhaseFunctionPolar(Ui_MainWindow *ui, parameters *para, QVector<double> theta, QVector<double> phaseFunction)
{
    //Polar plot: phase function
    QwtInterval radial( 0, para->maxPolarPtheta );
    QwtInterval azimuth( 0.0, 360.0 );
    ui->qwtpolarplot_pFunctionPolar->setScale(QwtPolar::Radius, 0, para->maxPolarPtheta);
    ui->qwtpolarplot_pFunctionPolar->setScale(QwtPolar::Azimuth, 0, 360, 360/12);

    para->polarCurve->detach();
    ui->qwtpolarplot_pFunctionPolar->replot();
    para->polarCurve->setStyle( QwtPolarCurve::Lines );
    para->polarCurve->setPen( QPen( Qt::red, 2 ) );
    para->polarCurve->setData(new Polar(radial,azimuth,2*para->nTheta-1,theta, phaseFunction));
    para->polarCurve->attach(ui->qwtpolarplot_pFunctionPolar);
    if (ui->radioButton_polarPlotScale_linear->isChecked())
        ui->qwtpolarplot_pFunctionPolar->setScaleEngine( QwtPolar::Radius, new QwtLinearScaleEngine() );
    if (ui->radioButton_polarPlotScale_log->isChecked())
    {
        ui->qwtpolarplot_pFunctionPolar->setScale(QwtPolar::Radius, para->minPolarPtheta, para->maxPolarPtheta);
        ui->qwtpolarplot_pFunctionPolar->setScaleEngine( QwtPolar::Radius, new QwtLogScaleEngine() );
    }
    ui->qwtpolarplot_pFunctionPolar->replot();
    if (para->nWavel==1)
        ui->slider_pFunctionPolar_wavel->setDisabled(true);
    else
        ui->slider_pFunctionPolar_wavel->setDisabled(false);
}

//Assign values for S1/S2 plot
void PlotData::AssignValuesPhaseFunctionLinearPlot(Ui_MainWindow *ui, parameters *para)
{
    double tiny = 1e-100;  //add a small number before log calculation to avoid NaN

    QVector<double> phaseFunctionPara(static_cast<int>(2*para->nTheta-1));
    QVector<double> phaseFunctionPerp(static_cast<int>(2*para->nTheta-1));
    QVector<double> phaseFunctionAve(static_cast<int>(2*para->nTheta-1));
    QVector<double> theta(static_cast<int>(2*para->nTheta-1));
    int indexWL = ui->slider_pFunctionLinear_wavel->value();
    ui->label_pFunctionLinear_curWavel->setText(QString::number(para->startWavel + indexWL*para->stepWavel));

    //Assign "rho(r)" and "theta" values in polar plot
    for (int i=0; i<static_cast<int>(para->nTheta); i++)
    {
        theta[i] = -180 + (180.0 * i /(para->nTheta-1));
        theta[static_cast<int>(2*para->nTheta)-2 -i] =  180 - (180.0 * i /(para->nTheta-1));

        phaseFunctionPara[static_cast<int>(para->nTheta)-1 -i] = para->phaseFunctionPara[indexWL][i];
        phaseFunctionPara[static_cast<int>(para->nTheta)-1 +i] = para->phaseFunctionPara[indexWL][i];

        phaseFunctionPerp[static_cast<int>(para->nTheta)-1 -i] = para->phaseFunctionPerp[indexWL][i];
        phaseFunctionPerp[static_cast<int>(para->nTheta)-1 +i] = para->phaseFunctionPerp[indexWL][i];

        phaseFunctionAve[static_cast<int>(para->nTheta)-1 -i] = para->phaseFunctionAve[indexWL][i];
        phaseFunctionAve[static_cast<int>(para->nTheta)-1 +i] = para->phaseFunctionAve[indexWL][i];

        if (ui->radioButton_logYaxis->isChecked())
        {
            phaseFunctionPara[static_cast<int>(para->nTheta)-1 -i] = log10(para->phaseFunctionPara[indexWL][i]+tiny);
            phaseFunctionPara[static_cast<int>(para->nTheta)-1 +i] = log10(para->phaseFunctionPara[indexWL][i]+tiny);

            phaseFunctionPerp[static_cast<int>(para->nTheta)-1 -i] = log10(para->phaseFunctionPerp[indexWL][i]+tiny);
            phaseFunctionPerp[static_cast<int>(para->nTheta)-1 +i] = log10(para->phaseFunctionPerp[indexWL][i]+tiny);

            phaseFunctionAve[static_cast<int>(para->nTheta)-1 -i] = log10(para->phaseFunctionAve[indexWL][i]+tiny);
            phaseFunctionAve[static_cast<int>(para->nTheta)-1 +i] = log10(para->phaseFunctionAve[indexWL][i]+tiny);
        }
    }
    PlotPhaseFunctionLinear(ui, theta, phaseFunctionPara, phaseFunctionPerp, phaseFunctionAve);
}

//Plot Forward /Backward scattering percentage plot
void PlotData::PlotPhaseFunctionLinear(Ui_MainWindow *ui, QVector<double> x, QVector<double> yPara,
                                       QVector<double> yPerp, QVector<double> yAve)
{
    //Set Min and Max values for the y-axis
    double minY, maxY;
    minY = 1e10;
    maxY = 0;
    for(int i = 0; i<x.size(); i++)
    {
        //Min search
        if (minY > yAve[i])
               minY = yAve[i];
        if (minY > yPara[i])
               minY = yPara[i];
        if (minY > yPerp[i])
               minY = yPerp[i];
        //Max search
        if (maxY < yAve[i])
               maxY = yAve[i];
        if (maxY < yPara[i])
               maxY = yPara[i];
        if (maxY < yPerp[i])
               maxY = yPerp[i];
    }
    maxY = maxY + 0.02*(maxY-minY);
    minY = minY - 0.02*(maxY-minY);

    //Clear previous graphs
    ui->customPlot_pFunctionLinear->removeGraph(2);
    ui->customPlot_pFunctionLinear->removeGraph(1);
    ui->customPlot_pFunctionLinear->removeGraph(0);
    ui->customPlot_pFunctionLinear->legend->clearItems();

    //Average Plot
    ui->customPlot_pFunctionLinear->addGraph();
    ui->customPlot_pFunctionLinear->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_pFunctionLinear->graph(0)->setData(x, yAve);
    ui->customPlot_pFunctionLinear->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    ui->customPlot_pFunctionLinear->graph(0)->setPen( QPen( Qt::red, 1 ) );
    ui->customPlot_pFunctionLinear->graph(0)->setName("Ave.");

    //Para plot
    ui->customPlot_pFunctionLinear->addGraph();
    ui->customPlot_pFunctionLinear->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_pFunctionLinear->graph(1)->setData(x, yPara);
    ui->customPlot_pFunctionLinear->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    ui->customPlot_pFunctionLinear->graph(1)->setPen( QPen( Qt::blue, 1 ) );
    ui->customPlot_pFunctionLinear->graph(1)->setName("Para.");

    //Perp plot
    ui->customPlot_pFunctionLinear->addGraph();
    ui->customPlot_pFunctionLinear->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_pFunctionLinear->graph(2)->setData(x, yPerp);
    ui->customPlot_pFunctionLinear->graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    ui->customPlot_pFunctionLinear->graph(2)->setPen( QPen( Qt::green, 1 ) );
    ui->customPlot_pFunctionLinear->graph(2)->setName("Perp.");

    if (ui->radioButton_linearYaxis->isChecked())
        ui->customPlot_pFunctionLinear->yAxis->setLabel("Magnitude");
    if (ui->radioButton_logYaxis->isChecked())
        ui->customPlot_pFunctionLinear->yAxis->setLabel("Log (Magnitude)");

    ui->customPlot_pFunctionLinear->yAxis->setRange(minY, maxY);
    ui->customPlot_pFunctionLinear->replot();

    //Legend    
    ui->customPlot_pFunctionLinear->legend->setVisible(true);
    ui->customPlot_pFunctionLinear->axisRect()->insetLayout()->setInsetAlignment(0,Qt::AlignTop|Qt::AlignRight);
    ui->customPlot_pFunctionLinear->replot();
}

//Assign values for S1/S2 plot
void PlotData::AssignValuesS1S2Plot(Ui_MainWindow *ui, parameters *para)
{
    QVector<double> S(static_cast<int>(para->nTheta));
    QVector<double> theta(static_cast<int>(para->nTheta));
    std::complex<double> *tempS;

    int indexWL = ui->slider_s1s2_wavel->value();
    ui->label_s1s2_curWavel->setText(QString::number(para->startWavel + indexWL*para->stepWavel));

    tempS = new std::complex<double> [para->nTheta];
    if (ui->radioButton_s1_amp->isChecked())
        tempS = para->S1[indexWL];
    if (ui->radioButton_s2_amp->isChecked())
        tempS = para->S2[indexWL];

     for (int i=0; i<static_cast<int>(para->nTheta); i++)
    {
        theta[i] = 180.0 * i /(para->nTheta-1);

        if (ui->radioButton_logYaxis->isChecked())
        {
             S[i] = log10(sqrt(tempS[i].real()*tempS[i].real() + tempS[i].imag()*tempS[i].imag()));
        }
        else
        {
            if (ui->radioButton_s1s2_abs->isChecked())
                S[i] = sqrt(tempS[i].real()*tempS[i].real() + tempS[i].imag()*tempS[i].imag());
            if (ui->radioButton_s1s2_real->isChecked())
                S[i] = tempS[i].real();
            if (ui->radioButton_s1s2_imag->isChecked())
                S[i] = tempS[i].imag();
        }
    }
    PlotS1S2(ui, para, theta, S);
}

//Plot S1S2 plots
void PlotData::PlotS1S2(Ui_MainWindow *ui, parameters *para, QVector<double> x, QVector<double> yS)
{
    //Clear previous graph
    ui->customPlot_s1s2->removeGraph(0);
    //Add new graph
    ui->customPlot_s1s2->addGraph();
    ui->customPlot_s1s2->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_s1s2->graph(0)->setData(x, yS);
    ui->customPlot_s1s2->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 3));
    ui->customPlot_s1s2->graph(0)->setPen( QPen( Qt::red, 1 ) );
    ui->customPlot_s1s2->graph(0)->rescaleAxes();
    if (ui->radioButton_linearYaxis->isChecked())
    {
        if (ui->radioButton_s1_amp->isChecked())
            ui->customPlot_s1s2->yAxis->setLabel("S1");
        if (ui->radioButton_s2_amp->isChecked())
            ui->customPlot_s1s2->yAxis->setLabel("S2");
    }
    if (ui->radioButton_logYaxis->isChecked())
    {
        if (ui->radioButton_s1_amp->isChecked())
            ui->customPlot_s1s2->yAxis->setLabel("Log (S1)");
        if (ui->radioButton_s2_amp->isChecked())
            ui->customPlot_s1s2->yAxis->setLabel("Log (S2)");
    }
    ui->customPlot_s1s2->replot();
    if (para->nWavel==1)
        ui->slider_s1s2_wavel->setDisabled(true);
    else
        ui->slider_s1s2_wavel->setDisabled(false);
}

//Assign values for Scattering cross section, g and Musp plots
void PlotData::AssignValuesAllOtherPlots(Ui_MainWindow *ui, parameters* para)
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
    double margin = (1.0 + ui->slider_concPercentChange->value() /200.0);

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
        if (ui->radioButton_monoDisperse->isChecked())
            ySizePara[i] = para->SizePara[i];

        if (ui->radioButton_logYaxis->isChecked())
        {
            yCsca[i] = log10(para->cSca[i] + tiny);    //scattering cross section
            yCext[i] = log10(para->cExt[i] + tiny);    //extinction cross section
            yCback[i] = log10(para->cBack[i] + tiny);  //backscattering cross section
            yG[i] = log10(para->g[i] +tiny);           //g - average cosine of phase function
            yMus[i] = log10(mus + tiny);              //scattering coefficient
            yMusp[i] = log10((mus * (1.0 - para->g[i])) + tiny);    //reduced scattering coefficient
            yF[i] = log10(para->forward[i] + tiny);
            yB[i] = log10(para->backward[i] + tiny);
            if (ui->radioButton_monoDisperse->isChecked())
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
    PlotMuspCurveForPowerLawFit(ui, para, x, yMusp);
    if (ui->radioButton_monoDisperse->isChecked())
        PlotSizeParameter(ui, x, ySizePara);
}

//Plot Scattering Cross Section plot
void PlotData::PlotScatteringCrossSection(Ui_MainWindow *ui, QVector<double> x, QVector<double> yCsca)
{
   //Clear previous graph
   ui->customPlot_csca->removeGraph(0);
   //Add new graph
   ui->customPlot_csca->addGraph();
   ui->customPlot_csca->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
   ui->customPlot_csca->graph(0)->setData(x, yCsca);
   ui->customPlot_csca->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
   ui->customPlot_csca->graph(0)->setPen( QPen( Qt::red, 1 ) );
   if (ui->radioButton_linearYaxis->isChecked())
       ui->customPlot_csca->yAxis->setLabel("Scattering Cross Section (μm²)");
   if (ui->radioButton_logYaxis->isChecked())
       ui->customPlot_csca->yAxis->setLabel("Log (Scattering Cross Section (μm²) )");
   ui->customPlot_csca->graph(0)->rescaleAxes();
   ui->customPlot_csca->replot();
}

//Plot Extinction Cross Section plot
void PlotData::PlotExtinctionCrossSection(Ui_MainWindow *ui, QVector<double> x, QVector<double> yCext)
{
   //Clear previous graph
   ui->customPlot_cext->removeGraph(0);
   //Add new graph
   ui->customPlot_cext->addGraph();
   ui->customPlot_cext->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
   ui->customPlot_cext->graph(0)->setData(x, yCext);
   ui->customPlot_cext->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
   ui->customPlot_cext->graph(0)->setPen( QPen( Qt::red, 1 ) );
   if (ui->radioButton_linearYaxis->isChecked())
       ui->customPlot_cext->yAxis->setLabel("Extinction Cross Section (μm²)");
   if (ui->radioButton_logYaxis->isChecked())
       ui->customPlot_cext->yAxis->setLabel("Log (Extinction Cross Section (μm²) )");
   ui->customPlot_cext->graph(0)->rescaleAxes();
   ui->customPlot_cext->replot();
}

//Plot Backscattering Cross Section plot
void PlotData::PlotBackscatteringCrossSection(Ui_MainWindow *ui, QVector<double> x, QVector<double> yCback)
{
   //Clear previous graph
   ui->customPlot_cback->removeGraph(0);
   //Add new graph
   ui->customPlot_cback->addGraph();
   ui->customPlot_cback->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
   ui->customPlot_cback->graph(0)->setData(x, yCback);
   ui->customPlot_cback->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
   ui->customPlot_cback->graph(0)->setPen( QPen( Qt::red, 1 ) );
   if (ui->radioButton_linearYaxis->isChecked())
       ui->customPlot_cback->yAxis->setLabel("Backscattering Cross Section (μm²)");
   if (ui->radioButton_logYaxis->isChecked())
       ui->customPlot_cback->yAxis->setLabel("Log (Backscattering Cross Section (μm²) )");
   ui->customPlot_cback->graph(0)->rescaleAxes();
   ui->customPlot_cback->replot();
}

//Plot Size Parameter
void PlotData::PlotSizeParameter(Ui_MainWindow *ui, QVector<double> x, QVector<double> ySizePara)
{
   //Clear previous graph
   ui->customPlot_sizePara->removeGraph(0);
   //Add new graph
   ui->customPlot_sizePara->addGraph();
   ui->customPlot_sizePara->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
   ui->customPlot_sizePara->graph(0)->setData(x, ySizePara);
   ui->customPlot_sizePara->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
   ui->customPlot_sizePara->graph(0)->setPen( QPen( Qt::red, 1 ) );
   if (ui->radioButton_linearYaxis->isChecked())
       ui->customPlot_sizePara->yAxis->setLabel("Size Parameter");
   if (ui->radioButton_logYaxis->isChecked())
       ui->customPlot_sizePara->yAxis->setLabel("Log (Size Parameter)");
   ui->customPlot_sizePara->graph(0)->rescaleAxes();
   ui->customPlot_sizePara->replot();
}

//Plot Mus plot
void PlotData::PlotMus(Ui_MainWindow *ui,QVector<double> x, QVector<double> yMus)
{
   //Clear previous graph
   ui->customPlot_mus->removeGraph(0);
   //Add new graph
   ui->customPlot_mus->addGraph();
   ui->customPlot_mus->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
   ui->customPlot_mus->graph(0)->setData(x, yMus);
   ui->customPlot_mus->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
   ui->customPlot_mus->graph(0)->setPen( QPen( Qt::red, 1 ) );
   if (ui->radioButton_linearYaxis->isChecked())
       ui->customPlot_mus->yAxis->setLabel("μs (mmˉˡ)");
   if (ui->radioButton_logYaxis->isChecked())
       ui->customPlot_mus->yAxis->setLabel("Log (μs (mmˉˡ) )");
   ui->customPlot_mus->graph(0)->rescaleAxes();
   ui->customPlot_mus->replot();
}

//Plot Musp plot
void PlotData::PlotMusp(Ui_MainWindow *ui, QVector<double> x, QVector<double> yMusp)
{
    ui->customPlot_musp->removeGraph(0);
    ui->customPlot_musp->addGraph();
    ui->customPlot_musp->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_musp->graph(0)->setData(x, yMusp);
    ui->customPlot_musp->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->customPlot_musp->graph(0)->setPen( QPen( Qt::red, 1 ) );
    if (ui->radioButton_linearYaxis->isChecked())
        ui->customPlot_musp->yAxis->setLabel("μs' (mmˉˡ)");
    if (ui->radioButton_logYaxis->isChecked())
        ui->customPlot_musp->yAxis->setLabel("Log (μs' (mmˉˡ) )");
    ui->customPlot_musp->graph(0)->rescaleAxes();
    ui->customPlot_musp->replot();
}

//Plot Musp for PowerLaw Fit
void PlotData::PlotMuspCurveForPowerLawFit(Ui_MainWindow *ui, parameters *para, QVector<double> x,
                                           QVector<double> yMusp)
{
    bool flag;

    //Clear previous graphs
    ui->customPlot_muspPowerLaw->removeGraph(1);
    ui->customPlot_muspPowerLaw->removeGraph(0);
    //Make another copy for muspfit
    ui->customPlot_muspPowerLaw->addGraph();
    ui->customPlot_muspPowerLaw->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_muspPowerLaw->graph(0)->setData(x, yMusp);
    ui->customPlot_muspPowerLaw->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->customPlot_muspPowerLaw->graph(0)->setPen( QPen( Qt::red, 1 ) );
    if (ui->radioButton_linearYaxis->isChecked())
        ui->customPlot_muspPowerLaw->yAxis->setLabel("μs' (mmˉˡ)");
    if (ui->radioButton_logYaxis->isChecked())
        ui->customPlot_muspPowerLaw->yAxis->setLabel("Log (μs' (mmˉˡ) )");
    ui->customPlot_muspPowerLaw->graph(0)->rescaleAxes();
    ui->customPlot_muspPowerLaw->replot();
    ui->customPlot_muspPowerLaw->addGraph();   //dummy

    if (para->nWavel==1)
        flag = true;
    else
        flag = false;

    ui->qwtslider_powerLaw_b->setDisabled(flag);
    ui->doubleSpinBox_powerLaw_b->setDisabled(flag);
    ui->pushButton_bestFit->setDisabled(flag);

    if ((!para->fittingComplex) || (para->nWavel==1))
    {
        ui->qwtslider_powerLaw_f->setDisabled(true);
        ui->doubleSpinBox_powerLaw_f->setDisabled(true);
    }
    else
    {
        ui->qwtslider_powerLaw_f->setDisabled(false);
        ui->doubleSpinBox_powerLaw_f->setDisabled(false);
    }
}

//Plot g (average phase function) plot
void PlotData::PlotG(Ui_MainWindow *ui, QVector<double> x, QVector<double> yG)
{
    //Clear previous graph
    ui->customPlot_g->removeGraph(0);
    //Add new graph
    ui->customPlot_g->addGraph();
    ui->customPlot_g->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_g->graph(0)->setData(x, yG);
    ui->customPlot_g->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->customPlot_g->graph(0)->setPen( QPen( Qt::red, 1 ) );
    ui->customPlot_g->graph(0)->rescaleAxes();
    if (ui->radioButton_linearYaxis->isChecked())
        ui->customPlot_g->yAxis->setLabel("g (Average Cosine of phase function)");
    if (ui->radioButton_logYaxis->isChecked())
        ui->customPlot_g->yAxis->setLabel("Log (g (Average Cosine of phase function))");
    ui->customPlot_g->graph(0)->rescaleAxes();
    ui->customPlot_g->replot();
}

//Plot Forward /Backward scattering percentage plot
void PlotData::PlotForwardBackward(Ui_MainWindow *ui, QVector<double> x, QVector<double> yF, QVector<double> yB, bool legendFlag)
{
    //Clear previous graphs
    ui->customPlot_forwardBackward->removeGraph(1);
    ui->customPlot_forwardBackward->removeGraph(0);
    ui->customPlot_forwardBackward->legend->clearItems();

    //Forward Plot
    ui->customPlot_forwardBackward->addGraph();
    ui->customPlot_forwardBackward->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_forwardBackward->graph(0)->setData(x, yF);
    ui->customPlot_forwardBackward->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->customPlot_forwardBackward->graph(0)->setPen( QPen( Qt::red, 1 ) );
    ui->customPlot_forwardBackward->graph(0)->setName("Forward Scat. %");
    ui->customPlot_forwardBackward->legend->setVisible(true);
    if (ui->radioButton_linearYaxis->isChecked())
    {
        ui->customPlot_forwardBackward->yAxis->setLabel("Forward & Backward Scattering %");
        ui->customPlot_forwardBackward->graph()->rescaleAxes();
        ui->customPlot_forwardBackward->yAxis->setRange(0,115);
    }
    if (ui->radioButton_logYaxis->isChecked())
    {
        ui->customPlot_forwardBackward->yAxis->setLabel("Log (Forward & Backward Scattering %)");
        ui->customPlot_forwardBackward->graph()->rescaleAxes();
        ui->customPlot_forwardBackward->yAxis->setRange(-3, 4);
    }
    //Backward plot
    ui->customPlot_forwardBackward->addGraph();
    ui->customPlot_forwardBackward->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_forwardBackward->graph(1)->setData(x, yB);
    ui->customPlot_forwardBackward->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->customPlot_forwardBackward->graph(1)->setPen( QPen( Qt::blue, 1 ) );
    ui->customPlot_forwardBackward->graph(1)->setName("Backward Scat. %");
    ui->customPlot_forwardBackward->legend->setVisible(true);
    if (legendFlag)
        ui->customPlot_forwardBackward->axisRect()->insetLayout()->setInsetAlignment(0,Qt::AlignCenter|Qt::AlignRight);
    else
        ui->customPlot_forwardBackward->axisRect()->insetLayout()->setInsetAlignment(0,Qt::AlignTop|Qt::AlignRight);
    ui->customPlot_forwardBackward->replot();
}

//Assign values for Musp power law plots
void PlotData::AssignValuesMuspPowerLawPlots(Ui_MainWindow *ui, parameters* para)
{
    QVector<double> x(static_cast<int>(para->nWavel));
    QVector<double> yMusp(static_cast<int>(para->nWavel));
    QVector<double> yFit(static_cast<int>(para->nWavel));
    double tiny = 1e-100;  //add a small number before log calculation to avoid NaN
    double tempError = 0.0;
    double error;
    double mus, fitA;
    double margin = (1.0 + ui->slider_concPercentChange->value() /200.0);

    for (int i=0; i<static_cast<int>(para->nWavel); i++)
    {
        mus = para->mus[i]*margin;  // set mus according to conc slider value

        x[i] = para->wavelArray[i];
        yMusp[i] = mus * (1.0 - para->g[i]);          //reduced scattering coefficient

        //Steve L Jacques,"Optical properties of biological tissues: a review" Phys. Med & Bio. 58(2013) R37-R61.
        //wavelength λ is normalized by a reference wavelength, 500 nm or 1000nm
        fitA = para->muspAtRefWavel *margin;
        if (para->fittingComplex)
            yFit[i] = fitA *(para->fRay*pow(x[i]/para->refWavel, -4.0) + (1-para->fRay)*pow(x[i]/para->refWavel, -para->bMie));
        else
            yFit[i] = fitA *pow(x[i]/para->refWavel, -para->bMie);   // A(lambda/lambdaRef)^-b

        error = yFit[i] - yMusp[i];
        tempError += error*error;

        if (ui->radioButton_logYaxis->isChecked())
        {
            yFit[i] = log10(yFit[i]+tiny);
            yMusp[i] = log10((mus * (1.0 - para->g[i])) + tiny);    //reduced scattering coefficient
        }
    }
    para->muspFittingError = tempError/para->nWavel;
    PlotMuspPowerLaw(ui, x,yMusp,yFit);
}

//Plot power law fit
void PlotData::PlotMuspPowerLaw(Ui_MainWindow *ui, QVector<double> x, QVector<double> yMusp, QVector<double> yFit)
{
    //Clear previous graphs
    ui->customPlot_muspPowerLaw->removeGraph(1);
    ui->customPlot_muspPowerLaw->removeGraph(0);
    ui->customPlot_muspPowerLaw->legend->clearItems();

    //Main plot
    ui->customPlot_muspPowerLaw->addGraph();
    ui->customPlot_muspPowerLaw->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_muspPowerLaw->graph(0)->setData(x, yMusp);
    ui->customPlot_muspPowerLaw->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->customPlot_muspPowerLaw->graph(0)->setPen( QPen( Qt::red, 1 ) );
    ui->customPlot_muspPowerLaw->graph(0)->setName("μs' Data");
    ui->customPlot_muspPowerLaw->legend->setVisible(true);
    if (ui->radioButton_linearYaxis->isChecked())
        ui->customPlot_muspPowerLaw->yAxis->setLabel("μs' (mmˉˡ)");
    if (ui->radioButton_logYaxis->isChecked())
        ui->customPlot_muspPowerLaw->yAxis->setLabel("Log (μs' (mmˉˡ))");

    //Power law fit plot
    ui->customPlot_muspPowerLaw->addGraph();
    ui->customPlot_muspPowerLaw->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_muspPowerLaw->graph(1)->setData(x,yFit);
    ui->customPlot_muspPowerLaw->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->customPlot_muspPowerLaw->graph(1)->setPen( QPen( Qt::blue, 1 ) );
    ui->customPlot_muspPowerLaw->graph(1)->setName("Best Fit");
    ui->customPlot_muspPowerLaw->legend->setVisible(true);
    ui->customPlot_muspPowerLaw->axisRect()->insetLayout()->setInsetAlignment(0,Qt::AlignTop|Qt::AlignRight);
    ui->customPlot_muspPowerLaw->graph()->rescaleAxes();
    ui->customPlot_muspPowerLaw->replot();
    ui->customPlot_muspPowerLaw->legend->setVisible(false);
}
