/**********************************************************************
** All QCustomPlot asignments and plotting actions are listed here
**********************************************************************/

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
    minX = 0.95 * (ui->lineEdit_Diameter->text().toDouble());
    maxX = 1.05 * (ui->lineEdit_Diameter->text().toDouble());

    ui->customPlot_Distribution->xAxis->setLabel("Diameter (μm)");
    ui->customPlot_Distribution->xAxis->setRange(minX, maxX);
    ui->customPlot_Distribution->yAxis->setLabel("Ns  (/vol. of 1mm³)");
    ui->customPlot_Distribution->clearPlottables();
    ui->customPlot_Distribution->replot();
}

//Initialize polar plot
void PlotData::InitializePhaseFunctionPolarPlot(Ui_MainWindow *ui, parameters *para )
{
//    SetSliderWL(ui);

    //Polar plot: phase function
    ui->qwtpolarplot_PhaseFunctionPolar->setAutoReplot( false );
    ui->qwtpolarplot_PhaseFunctionPolar->setPlotBackground( Qt::white );
    ui->qwtpolarplot_PhaseFunctionPolar->setScale( QwtPolar::Azimuth,0, 360, 360 / 12 );
    ui->qwtpolarplot_PhaseFunctionPolar->setScale( QwtPolar::Radius, 0, 1);
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
    mGrid->attach( ui->qwtpolarplot_PhaseFunctionPolar );
    para->polarCurve->detach();
    ui->qwtpolarplot_PhaseFunctionPolar->replot();
}

//Initialize distribution plot
void PlotData::InitializePhaseFunctionLinearPlot(Ui_MainWindow *ui)
{
    ui->customPlot_PhaseFunctionLinear->xAxis->setLabel("Angle (deg.)");
    ui->customPlot_PhaseFunctionLinear->xAxis->setRange(-180,180);
    ui->customPlot_PhaseFunctionLinear->xAxis->ticker()->setTickStepStrategy(QCPAxisTicker::TickStepStrategy::tssMeetTickCount);
    QVector<double> ticks;
    QVector<QString> labels;
    ticks << -180 << -135 << -90 << -45 << 0 << 45 << 90 << 135 << 180;
    labels << "-180" << "135" << "-90" << "-45" << "0" << "45" << "90" << "135" << "180";
    QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
    textTicker->addTicks(ticks, labels);
    ui->customPlot_PhaseFunctionLinear->xAxis->setTicker(textTicker);
    ui->customPlot_PhaseFunctionLinear->yAxis->setLabel("Magnitude");
    ui->customPlot_PhaseFunctionLinear->clearGraphs();
    ui->customPlot_PhaseFunctionLinear->replot();
    ui->customPlot_PhaseFunctionLinear->addGraph(); //dummy
    ui->customPlot_PhaseFunctionLinear->addGraph(); //dummy
    ui->customPlot_PhaseFunctionLinear->addGraph(); //dummy
    ui->customPlot_PhaseFunctionLinear->legend->setVisible(false);
}

//Initialize all other plots
void PlotData::InitializeAllOtherPlots(Ui_MainWindow *ui)
{
    double minX = ui->lineEdit_StartWL->text().toDouble();
    double maxX = ui->lineEdit_EndWL->text().toDouble();

    //Plot S1S2
    ui->customPlot_S1S2->xAxis->setLabel("Angle (deg.)");
    ui->customPlot_S1S2->xAxis->setRange(0,180);
    ui->customPlot_S1S2->xAxis->ticker()->setTickStepStrategy(QCPAxisTicker::TickStepStrategy::tssMeetTickCount);
    QVector<double> ticks;
    QVector<QString> labels;
    ticks << 0 << 45 << 90 << 135 << 180;
    labels << "0" << "45" << "90" << "135" << "180";
    QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
    textTicker->addTicks(ticks, labels);
    ui->customPlot_S1S2->xAxis->setTicker(textTicker);
    ui->customPlot_S1S2->yAxis->setLabel("S1");
    ui->customPlot_S1S2->clearGraphs();
    ui->customPlot_S1S2->replot();
    ui->customPlot_S1S2->addGraph();   //dummy

    //Plot Scattering Cross Section
    ui->customPlot_Csca->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_Csca->xAxis->setRange(minX, maxX);
    ui->customPlot_Csca->yAxis->setLabel("Scattering Cross Section (μm²)");
    ui->customPlot_Csca->clearGraphs();
    ui->customPlot_Csca->replot();
    ui->customPlot_Csca->addGraph();   //dummy

    //Plot Extinction Cross Section
    ui->customPlot_Cext->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_Cext->xAxis->setRange(minX, maxX);
    ui->customPlot_Cext->yAxis->setLabel("Extinction Cross Section (μm²)");
    ui->customPlot_Cext->clearGraphs();
    ui->customPlot_Cext->replot();
    ui->customPlot_Cext->addGraph();   //dummy

    //Plot Backscattering Cross Section
    ui->customPlot_Cback->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_Cback->xAxis->setRange(minX, maxX);
    ui->customPlot_Cback->yAxis->setLabel("Backscattering Cross Section (μm²)");
    ui->customPlot_Cback->clearGraphs();
    ui->customPlot_Cback->replot();
    ui->customPlot_Cback->addGraph();   //dummy

    //Plot Size Parameter
    ui->customPlot_SizePara->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_SizePara->xAxis->setRange(minX, maxX);
    ui->customPlot_SizePara->yAxis->setLabel("Size Parameter");
    ui->customPlot_SizePara->clearGraphs();
    ui->customPlot_SizePara->replot();
    ui->customPlot_SizePara->addGraph();   //dummy

    //Plot Mus
    ui->customPlot_Mus->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_Mus->xAxis->setRange(minX, maxX);
    ui->customPlot_Mus->yAxis->setLabel("μs (mmˉˡ)");
    ui->customPlot_Mus->clearGraphs();
    ui->customPlot_Mus->replot();
    ui->customPlot_Mus->addGraph();   //dummy

    //Plot G
    ui->customPlot_G->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_G->xAxis->setRange(minX, maxX);
    ui->customPlot_G->yAxis->setLabel("g (Average Cosine of phase function)");
    ui->customPlot_G->clearGraphs();
    ui->customPlot_G->replot();
    ui->customPlot_G->addGraph();   //dummy

    //Plot FB
    ui->customPlot_FB->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_FB->xAxis->setRange(minX, maxX);
    ui->customPlot_FB->yAxis->setLabel("Forward & Backward Scattering %");
    ui->customPlot_FB->clearGraphs();
    ui->customPlot_FB->replot();
    ui->customPlot_FB->addGraph(); //dummy
    ui->customPlot_FB->addGraph(); //dummy
    ui->customPlot_FB->legend->setVisible(false);

    //Plot Musp
    ui->customPlot_Musp->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_Musp->xAxis->setRange(minX, maxX);
    ui->customPlot_Musp->yAxis->setLabel("μs' (mmˉˡ)");
    ui->customPlot_Musp->clearGraphs();
    ui->customPlot_Musp->replot();
    ui->customPlot_Musp->addGraph();  //dummy

    //Plot Musp for fitting
    ui->customPlot_MuspPowerLaw->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_MuspPowerLaw->xAxis->setRange(minX, maxX);
    ui->customPlot_MuspPowerLaw->yAxis->setLabel("μs' (mmˉˡ)");
    ui->customPlot_MuspPowerLaw->clearGraphs();
    ui->customPlot_MuspPowerLaw->replot();
    ui->customPlot_MuspPowerLaw->addGraph();  //dummy
    ui->customPlot_MuspPowerLaw->addGraph();  //dummy
    ui->customPlot_MuspPowerLaw->legend->setVisible(false);
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
    QVector<double> xDist(para->nRadius);
    QVector<double> yDist(para->nRadius);
    double tiny = 1e-100;  //add a small number before log calculation to avoid NaN
    double margin = (1.0 + ui->slider_ConcPercentChange->value()/200.0);

    double totNumDen = 0.0;
    for (int i=0; i<para->nRadius; i++)
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

//plot distribution plot
void PlotData::PlotDistribution(Ui_MainWindow *ui, parameters *para, QVector<double> x, QVector<double> yDist)
{
    double minX, maxX, minY, maxY;

    minY = 0;
    maxY = yDist[0];
    for(int i = 0; i<para->nRadius; i++)
    {
        if (maxY < yDist[i])
               maxY = yDist[i];
    }
    maxY *= 1.1;

    ui->customPlot_Distribution->removeGraph(0);
    ui->customPlot_Distribution->clearPlottables();
    ui->customPlot_Distribution->replot();
    ui->customPlot_Distribution->addGraph();
    ui->customPlot_Distribution->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);

    ui->customPlot_Distribution->graph()->setData(x, yDist);
    ui->customPlot_Distribution->graph()->setPen( QPen( Qt::red, 1 ) );
    ui->customPlot_Distribution->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 3));

    QCPBars *bar = new QCPBars(ui->customPlot_Distribution->xAxis, ui->customPlot_Distribution->yAxis);
    bar->setData(x, yDist);
    ui->customPlot_Distribution->yAxis->setRange(minY, maxY);
    bar->setPen(Qt::NoPen);
    bar->setBrush(QColor(0,0,0,255));

    if (ui->radioButton_MonoDisperse->isChecked())
    {
        minX = 0.92 * 2.0 * para->minRadius;
        maxX = 1.08 * 2.0 * para->maxRadius;
        bar->setWidth(0.05*(maxX -minX));
        ui->customPlot_Distribution->xAxis->setRange(minX, maxX);        
        ui->customPlot_Distribution->xAxis->ticker()->setTickCount(3);
    }
    if (ui->radioButton_PolyDisperse->isChecked())
    {
        double shift = (para->maxRadius - para->minRadius)/(para->nRadius);
        if (shift ==0)      //Array of spherical scatterers with similar radius
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
            ui->customPlot_Distribution->xAxis->ticker()->setTickCount(5);
        }
        ui->customPlot_Distribution->xAxis->setRange(minX, maxX);
    }

    if (ui->radioButton_LinearYAxis->isChecked())
        ui->customPlot_Distribution->yAxis->setLabel("Concentration (Number Density: Ns)");
    if (ui->radioButton_LogYAxis->isChecked())
        ui->customPlot_Distribution->yAxis->setLabel("Log (Concentration (Number Density: Ns) )");
    if (ui->radioButton_LogXAxis->isChecked())
    {
        ui->customPlot_Distribution->xAxis->setScaleType(QCPAxis::stLogarithmic);
        ui->customPlot_Distribution->xAxis->setLabel("Diameter (μm) {Axis in Log Scale}");
    }
    if (ui->radioButton_LinearXAxis->isChecked())
    {
        ui->customPlot_Distribution->xAxis->setScaleType(QCPAxis::stLinear);
        ui->customPlot_Distribution->xAxis->setLabel("Diameter (μm)");
    }
    ui->customPlot_Distribution->replot();
}

//Assign values for phase function polar plot
void PlotData::AssignValuesPhaseFunctionPolarPlot(Ui_MainWindow *ui, parameters *para)
{
    QVector<double> phaseFunction(2*para->nTheta-1);
    QVector<double> theta(2*para->nTheta-1);
    int indexWL = ui->slider_WL_PFPolar->value();
    ui->label_CurrentWL_PFPolar->setText(QString::number(para->startWavel + indexWL*para->stepWavel));

    //Set minimum and maximum values
    para->maxPolarPtheta = 0;
    para->minPolarPtheta = 1e100;
    for (int i = 0; i < para->nWavel; i++)
    {
        for (int j = 0; j < para->nTheta; j++)
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
    for (int i=0; i<para->nTheta; i++)
    {
        theta[i] = 180.0 * i /(para->nTheta-1);
        theta[2*para->nTheta-2 -i] =  360.0 - (180.0 * i /(para->nTheta-1));
        if (ui->radioButton_PhaseAverage->isChecked())
        {
            phaseFunction[i] = para->phaseFunctionAve[indexWL][i];
            phaseFunction[2*para->nTheta-2 -i] = para->phaseFunctionAve[indexWL][i];
        }
        if (ui->radioButton_PhasePara->isChecked())
        {
            phaseFunction[i] = para->phaseFunctionPara[indexWL][i];
            phaseFunction[2*para->nTheta-2 -i] = para->phaseFunctionPara[indexWL][i];
        }
        if (ui->radioButton_PhasePerp->isChecked())
        {
            phaseFunction[i] = para->phaseFunctionPerp[indexWL][i];
            phaseFunction[2*para->nTheta-2 -i] = para->phaseFunctionPerp[indexWL][i];
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
    ui->qwtpolarplot_PhaseFunctionPolar->setScale(QwtPolar::Radius, 0, para->maxPolarPtheta);
    ui->qwtpolarplot_PhaseFunctionPolar->setScale(QwtPolar::Azimuth, 0, 360, 360/12);

    para->polarCurve->detach();
    ui->qwtpolarplot_PhaseFunctionPolar->replot();
    para->polarCurve->setStyle( QwtPolarCurve::Lines );
    para->polarCurve->setPen( QPen( Qt::red, 2 ) );
    para->polarCurve->setData(new Polar(radial,azimuth,2*para->nTheta-1,theta, phaseFunction));
    para->polarCurve->attach(ui->qwtpolarplot_PhaseFunctionPolar);
    if (ui->radioButton_PhaseLinear->isChecked())
        ui->qwtpolarplot_PhaseFunctionPolar->setScaleEngine( QwtPolar::Radius, new QwtLinearScaleEngine() );
    if (ui->radioButton_PhaseLog->isChecked())
    {
        ui->qwtpolarplot_PhaseFunctionPolar->setScale(QwtPolar::Radius, para->minPolarPtheta, para->maxPolarPtheta);
        ui->qwtpolarplot_PhaseFunctionPolar->setScaleEngine( QwtPolar::Radius, new QwtLogScaleEngine() );
    }
    ui->qwtpolarplot_PhaseFunctionPolar->replot();
    if (para->nWavel==1)
        ui->slider_WL_PFPolar->setDisabled(true);
    else
        ui->slider_WL_PFPolar->setDisabled(false);
}

//Assign values for S1/S2 plot
void PlotData::AssignValuesPhaseFunctionLinearPlot(Ui_MainWindow *ui, parameters *para)
{
    double tiny = 1e-100;  //add a small number before log calculation to avoid NaN

    QVector<double> phaseFunctionPara(2*para->nTheta-1);
    QVector<double> phaseFunctionPerp(2*para->nTheta-1);
    QVector<double> phaseFunctionAve(2*para->nTheta-1);
    QVector<double> theta(2*para->nTheta-1);
    int indexWL = ui->slider_WL_PFLinear->value();
    ui->label_CurrentWL_PFLinear->setText(QString::number(para->startWavel + indexWL*para->stepWavel));

    //Assign "rho(r)" and "theta" values in polar plot
    for (int i=0; i<para->nTheta; i++)
    {
        theta[i] = -180 + (180.0 * i /(para->nTheta-1));
        theta[2*para->nTheta-2 -i] =  180 - (180.0 * i /(para->nTheta-1));

        phaseFunctionPara[para->nTheta-1 -i] = para->phaseFunctionPara[indexWL][i];
        phaseFunctionPara[para->nTheta-1 +i] = para->phaseFunctionPara[indexWL][i];

        phaseFunctionPerp[para->nTheta-1 -i] = para->phaseFunctionPerp[indexWL][i];
        phaseFunctionPerp[para->nTheta-1 +i] = para->phaseFunctionPerp[indexWL][i];

        phaseFunctionAve[para->nTheta-1 -i] = para->phaseFunctionAve[indexWL][i];
        phaseFunctionAve[para->nTheta-1 +i] = para->phaseFunctionAve[indexWL][i];

        if (ui->radioButton_LogYAxis->isChecked())
        {
            phaseFunctionPara[para->nTheta-1 -i] = log10(para->phaseFunctionPara[indexWL][i]+tiny);
            phaseFunctionPara[para->nTheta-1 +i] = log10(para->phaseFunctionPara[indexWL][i]+tiny);

            phaseFunctionPerp[para->nTheta-1 -i] = log10(para->phaseFunctionPerp[indexWL][i]+tiny);
            phaseFunctionPerp[para->nTheta-1 +i] = log10(para->phaseFunctionPerp[indexWL][i]+tiny);

            phaseFunctionAve[para->nTheta-1 -i] = log10(para->phaseFunctionAve[indexWL][i]+tiny);
            phaseFunctionAve[para->nTheta-1 +i] = log10(para->phaseFunctionAve[indexWL][i]+tiny);
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
    ui->customPlot_PhaseFunctionLinear->removeGraph(2);
    ui->customPlot_PhaseFunctionLinear->removeGraph(1);
    ui->customPlot_PhaseFunctionLinear->removeGraph(0);
    ui->customPlot_PhaseFunctionLinear->legend->clearItems();

    //Average Plot
    ui->customPlot_PhaseFunctionLinear->addGraph();
    ui->customPlot_PhaseFunctionLinear->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_PhaseFunctionLinear->graph(0)->setData(x, yAve);
    ui->customPlot_PhaseFunctionLinear->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    ui->customPlot_PhaseFunctionLinear->graph(0)->setPen( QPen( Qt::red, 1 ) );
    ui->customPlot_PhaseFunctionLinear->graph(0)->setName("Ave.");

    //Para plot
    ui->customPlot_PhaseFunctionLinear->addGraph();
    ui->customPlot_PhaseFunctionLinear->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_PhaseFunctionLinear->graph(1)->setData(x, yPara);
    ui->customPlot_PhaseFunctionLinear->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    ui->customPlot_PhaseFunctionLinear->graph(1)->setPen( QPen( Qt::blue, 1 ) );
    ui->customPlot_PhaseFunctionLinear->graph(1)->setName("Para.");

    //Perp plot
    ui->customPlot_PhaseFunctionLinear->addGraph();
    ui->customPlot_PhaseFunctionLinear->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_PhaseFunctionLinear->graph(2)->setData(x, yPerp);
    ui->customPlot_PhaseFunctionLinear->graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 2));
    ui->customPlot_PhaseFunctionLinear->graph(2)->setPen( QPen( Qt::green, 1 ) );
    ui->customPlot_PhaseFunctionLinear->graph(2)->setName("Perp.");

    if (ui->radioButton_LinearYAxis->isChecked())    
        ui->customPlot_PhaseFunctionLinear->yAxis->setLabel("Magnitude");
    if (ui->radioButton_LogYAxis->isChecked())    
        ui->customPlot_PhaseFunctionLinear->yAxis->setLabel("Log (Magnitude)");

    ui->customPlot_PhaseFunctionLinear->yAxis->setRange(minY, maxY);
    ui->customPlot_PhaseFunctionLinear->replot();

    //Legend    
    ui->customPlot_PhaseFunctionLinear->legend->setVisible(true);
    ui->customPlot_PhaseFunctionLinear->axisRect()->insetLayout()->setInsetAlignment(0,Qt::AlignTop|Qt::AlignRight);
    ui->customPlot_PhaseFunctionLinear->replot();
}

//Assign values for S1/S2 plot
void PlotData::AssignValuesS1S2Plot(Ui_MainWindow *ui, parameters *para)
{
    QVector<double> S(para->nTheta);
    QVector<double> theta(para->nTheta);
    std::complex<double> *tempS;

    int indexWL = ui->slider_WL_S1S2->value();
    ui->label_CurrentWL_S1S2->setText(QString::number(para->startWavel + indexWL*para->stepWavel));

    tempS = new std::complex<double> [para->nTheta];
    if (ui->radioButton_S1->isChecked())
        tempS = para->S1[indexWL];
    if (ui->radioButton_S2->isChecked())
        tempS = para->S2[indexWL];

     for (int i=0; i<para->nTheta; i++)
    {
        theta[i] = 180.0 * i /(para->nTheta-1);

        if (ui->radioButton_LogYAxis->isChecked())
        {
             S[i] = log10(sqrt(tempS[i].real()*tempS[i].real() + tempS[i].imag()*tempS[i].imag()));
        }
        else
        {
            if (ui->radioButton_S1S2_Abs->isChecked())
                S[i] = sqrt(tempS[i].real()*tempS[i].real() + tempS[i].imag()*tempS[i].imag());
            if (ui->radioButton_S1S2_Real->isChecked())
                S[i] = tempS[i].real();
            if (ui->radioButton_S1S2_Imag->isChecked())
                S[i] = tempS[i].imag();
        }
    }
    PlotS1S2(ui, para, theta, S);
}

//Plot S1S2 plots
void PlotData::PlotS1S2(Ui_MainWindow *ui, parameters *para, QVector<double> x, QVector<double> yS)
{
    //Clear previous graph
    ui->customPlot_S1S2->removeGraph(0);
    //Add new graph
    ui->customPlot_S1S2->addGraph();
    ui->customPlot_S1S2->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_S1S2->graph(0)->setData(x, yS);
    ui->customPlot_S1S2->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 3));
    ui->customPlot_S1S2->graph(0)->setPen( QPen( Qt::red, 1 ) );
    ui->customPlot_S1S2->graph(0)->rescaleAxes();
    if (ui->radioButton_LinearYAxis->isChecked())
    {
        if (ui->radioButton_S1->isChecked())
            ui->customPlot_S1S2->yAxis->setLabel("S1");
        if (ui->radioButton_S2->isChecked())
            ui->customPlot_S1S2->yAxis->setLabel("S2");
    }
    if (ui->radioButton_LogYAxis->isChecked())
    {
        if (ui->radioButton_S1->isChecked())
            ui->customPlot_S1S2->yAxis->setLabel("Log (S1)");
        if (ui->radioButton_S2->isChecked())
            ui->customPlot_S1S2->yAxis->setLabel("Log (S2)");
    }
    ui->customPlot_S1S2->replot();
    if (para->nWavel==1)
        ui->slider_WL_S1S2->setDisabled(true);
    else
        ui->slider_WL_S1S2->setDisabled(false);
}

//Assign values for Scattering cross section, g and Musp plots
void PlotData::AssignValuesAllOtherPlots(Ui_MainWindow *ui, parameters* para)
{
    QVector<double> x(para->nWavel);
    QVector<double> yCsca(para->nWavel), yCext(para->nWavel), yCback(para->nWavel), ySizePara(para->nWavel);
    QVector<double> yG(para->nWavel),yMus(para->nWavel), yMusp(para->nWavel);
    QVector<double> yF(para->nWavel), yB(para->nWavel);

    double tiny = 1e-100;  //add a small number before log calculation to avoid NaN
    double mus;
    double margin = (1.0 + ui->slider_ConcPercentChange->value() /200.0);

    int fbLegendCheckLocation = 0.75*para->nWavel;
    bool fbLegnedFlag;
    double fbLimit =75;   //75%

    for (int i=0; i<para->nWavel; i++)
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
    PlotMuspCurveForPowerLawFit(ui, para, x, yMusp);
    if (ui->radioButton_MonoDisperse->isChecked())
        PlotSizeParameter(ui, x, ySizePara);
}

//Plot Scattering Cross Section plot
void PlotData::PlotScatteringCrossSection(Ui_MainWindow *ui, QVector<double> x, QVector<double> yCsca)
{
   //Clear previous graph
   ui->customPlot_Csca->removeGraph(0);
   //Add new graph
   ui->customPlot_Csca->addGraph();
   ui->customPlot_Csca->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
   ui->customPlot_Csca->graph(0)->setData(x, yCsca);
   ui->customPlot_Csca->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
   ui->customPlot_Csca->graph(0)->setPen( QPen( Qt::red, 1 ) );
   if (ui->radioButton_LinearYAxis->isChecked())
       ui->customPlot_Csca->yAxis->setLabel("Scattering Cross Section (μm²)");
   if (ui->radioButton_LogYAxis->isChecked())
       ui->customPlot_Csca->yAxis->setLabel("Log (Scattering Cross Section (μm²) )");
   ui->customPlot_Csca->graph(0)->rescaleAxes();
   ui->customPlot_Csca->replot();
}

//Plot Extinction Cross Section plot
void PlotData::PlotExtinctionCrossSection(Ui_MainWindow *ui, QVector<double> x, QVector<double> yCext)
{
   //Clear previous graph
   ui->customPlot_Cext->removeGraph(0);
   //Add new graph
   ui->customPlot_Cext->addGraph();
   ui->customPlot_Cext->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
   ui->customPlot_Cext->graph(0)->setData(x, yCext);
   ui->customPlot_Cext->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
   ui->customPlot_Cext->graph(0)->setPen( QPen( Qt::red, 1 ) );
   if (ui->radioButton_LinearYAxis->isChecked())
       ui->customPlot_Cext->yAxis->setLabel("Extinction Cross Section (μm²)");
   if (ui->radioButton_LogYAxis->isChecked())
       ui->customPlot_Cext->yAxis->setLabel("Log (Extinction Cross Section (μm²) )");
   ui->customPlot_Cext->graph(0)->rescaleAxes();
   ui->customPlot_Cext->replot();
}

//Plot Backscattering Cross Section plot
void PlotData::PlotBackscatteringCrossSection(Ui_MainWindow *ui, QVector<double> x, QVector<double> yCback)
{
   //Clear previous graph
   ui->customPlot_Cback->removeGraph(0);
   //Add new graph
   ui->customPlot_Cback->addGraph();
   ui->customPlot_Cback->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
   ui->customPlot_Cback->graph(0)->setData(x, yCback);
   ui->customPlot_Cback->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
   ui->customPlot_Cback->graph(0)->setPen( QPen( Qt::red, 1 ) );
   if (ui->radioButton_LinearYAxis->isChecked())
       ui->customPlot_Cback->yAxis->setLabel("Backscattering Cross Section (μm²)");
   if (ui->radioButton_LogYAxis->isChecked())
       ui->customPlot_Cback->yAxis->setLabel("Log (Backscattering Cross Section (μm²) )");
   ui->customPlot_Cback->graph(0)->rescaleAxes();
   ui->customPlot_Cback->replot();
}

//Plot Size Parameter
void PlotData::PlotSizeParameter(Ui_MainWindow *ui, QVector<double> x, QVector<double> ySizePara)
{
   //Clear previous graph
   ui->customPlot_SizePara->removeGraph(0);
   //Add new graph
   ui->customPlot_SizePara->addGraph();
   ui->customPlot_SizePara->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
   ui->customPlot_SizePara->graph(0)->setData(x, ySizePara);
   ui->customPlot_SizePara->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
   ui->customPlot_SizePara->graph(0)->setPen( QPen( Qt::red, 1 ) );
   if (ui->radioButton_LinearYAxis->isChecked())
       ui->customPlot_SizePara->yAxis->setLabel("Size Parameter");
   if (ui->radioButton_LogYAxis->isChecked())
       ui->customPlot_SizePara->yAxis->setLabel("Log (Size Parameter)");
   ui->customPlot_SizePara->graph(0)->rescaleAxes();
   ui->customPlot_SizePara->replot();
}

//Plot Mus plot
void PlotData::PlotMus(Ui_MainWindow *ui,QVector<double> x, QVector<double> yMus)
{
   //Clear previous graph
   ui->customPlot_Mus->removeGraph(0);
   //Add new graph
   ui->customPlot_Mus->addGraph();
   ui->customPlot_Mus->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
   ui->customPlot_Mus->graph(0)->setData(x, yMus);
   ui->customPlot_Mus->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
   ui->customPlot_Mus->graph(0)->setPen( QPen( Qt::red, 1 ) );
   if (ui->radioButton_LinearYAxis->isChecked())
       ui->customPlot_Mus->yAxis->setLabel("μs (mmˉˡ)");
   if (ui->radioButton_LogYAxis->isChecked())
       ui->customPlot_Mus->yAxis->setLabel("Log (μs (mmˉˡ) )");
   ui->customPlot_Mus->graph(0)->rescaleAxes();
   ui->customPlot_Mus->replot();
}

//Plot Musp plot
void PlotData::PlotMusp(Ui_MainWindow *ui, QVector<double> x, QVector<double> yMusp)
{
    ui->customPlot_Musp->removeGraph(0);
    ui->customPlot_Musp->addGraph();
    ui->customPlot_Musp->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_Musp->graph(0)->setData(x, yMusp);
    ui->customPlot_Musp->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->customPlot_Musp->graph(0)->setPen( QPen( Qt::red, 1 ) );
    if (ui->radioButton_LinearYAxis->isChecked())
        ui->customPlot_Musp->yAxis->setLabel("μs' (mmˉˡ)");
    if (ui->radioButton_LogYAxis->isChecked())
        ui->customPlot_Musp->yAxis->setLabel("Log (μs' (mmˉˡ) )");
    ui->customPlot_Musp->graph(0)->rescaleAxes();
    ui->customPlot_Musp->replot();
}

//Plot Musp for PowerLaw Fit
void PlotData::PlotMuspCurveForPowerLawFit(Ui_MainWindow *ui, parameters *para, QVector<double> x,
                                           QVector<double> yMusp)
{
    bool flag;

    //Clear previous graphs
    ui->customPlot_MuspPowerLaw->removeGraph(1);
    ui->customPlot_MuspPowerLaw->removeGraph(0);
    //Make another copy for muspfit
    ui->customPlot_MuspPowerLaw->addGraph();
    ui->customPlot_MuspPowerLaw->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_MuspPowerLaw->graph(0)->setData(x, yMusp);
    ui->customPlot_MuspPowerLaw->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->customPlot_MuspPowerLaw->graph(0)->setPen( QPen( Qt::red, 1 ) );
    if (ui->radioButton_LinearYAxis->isChecked())
        ui->customPlot_MuspPowerLaw->yAxis->setLabel("μs' (mmˉˡ)");
    if (ui->radioButton_LogYAxis->isChecked())
        ui->customPlot_MuspPowerLaw->yAxis->setLabel("Log (μs' (mmˉˡ) )");
    ui->customPlot_MuspPowerLaw->graph(0)->rescaleAxes();
    ui->customPlot_MuspPowerLaw->replot();
    ui->customPlot_MuspPowerLaw->addGraph();   //dummy

    if (para->nWavel==1)
        flag = true;
    else
        flag = false;

    ui->qwtslider_B->setDisabled(flag);
    ui->qwtslider_F->setDisabled(flag);
    ui->doubleSpinBox_B->setDisabled(flag);
    ui->doubleSpinBox_F->setDisabled(flag);
    ui->pushButton_BestFit->setDisabled(flag);
}

//Plot g (average phase function) plot
void PlotData::PlotG(Ui_MainWindow *ui, QVector<double> x, QVector<double> yG)
{
    //Clear previous graph
    ui->customPlot_G->removeGraph(0);
    //Add new graph
    ui->customPlot_G->addGraph();
    ui->customPlot_G->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_G->graph(0)->setData(x, yG);
    ui->customPlot_G->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->customPlot_G->graph(0)->setPen( QPen( Qt::red, 1 ) );
    ui->customPlot_G->graph(0)->rescaleAxes();
    if (ui->radioButton_LinearYAxis->isChecked())
        ui->customPlot_G->yAxis->setLabel("g (Average Cosine of phase function)");
    if (ui->radioButton_LogYAxis->isChecked())
        ui->customPlot_G->yAxis->setLabel("Log (g (Average Cosine of phase function))");
    ui->customPlot_G->graph(0)->rescaleAxes();
    ui->customPlot_G->replot();
}

//Plot Forward /Backward scattering percentage plot
void PlotData::PlotForwardBackward(Ui_MainWindow *ui, QVector<double> x, QVector<double> yF, QVector<double> yB, bool legendFlag)
{
    //Clear previous graphs
    ui->customPlot_FB->removeGraph(1);
    ui->customPlot_FB->removeGraph(0);
    ui->customPlot_FB->legend->clearItems();

    //Forward Plot
    ui->customPlot_FB->addGraph();
    ui->customPlot_FB->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_FB->graph(0)->setData(x, yF);
    ui->customPlot_FB->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->customPlot_FB->graph(0)->setPen( QPen( Qt::red, 1 ) );
    ui->customPlot_FB->graph(0)->setName("Forward Scat. %");
    ui->customPlot_FB->legend->setVisible(true);
    if (ui->radioButton_LinearYAxis->isChecked())
    {
        ui->customPlot_FB->yAxis->setLabel("Forward & Backward Scattering %");
        ui->customPlot_FB->graph()->rescaleAxes();
        ui->customPlot_FB->yAxis->setRange(0,115);
    }
    if (ui->radioButton_LogYAxis->isChecked())
    {
        ui->customPlot_FB->yAxis->setLabel("Log (Forward & Backward Scattering %)");
        ui->customPlot_FB->graph()->rescaleAxes();
        ui->customPlot_FB->yAxis->setRange(-3, 4);
    }
    //Backward plot
    ui->customPlot_FB->addGraph();
    ui->customPlot_FB->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_FB->graph(1)->setData(x, yB);
    ui->customPlot_FB->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->customPlot_FB->graph(1)->setPen( QPen( Qt::blue, 1 ) );
    ui->customPlot_FB->graph(1)->setName("Backward Scat. %");
    ui->customPlot_FB->legend->setVisible(true);
    if (legendFlag)
        ui->customPlot_FB->axisRect()->insetLayout()->setInsetAlignment(0,Qt::AlignCenter|Qt::AlignRight);
    else
        ui->customPlot_FB->axisRect()->insetLayout()->setInsetAlignment(0,Qt::AlignTop|Qt::AlignRight);
    ui->customPlot_FB->replot();
}

//Assign values for Musp power law plots
void PlotData::AssignValuesMuspPowerLawPlots(Ui_MainWindow *ui, parameters* para)
{
    QVector<double> x(para->nWavel), yMusp(para->nWavel);
    QVector<double> yFit(para->nWavel);
    double tiny = 1e-100;  //add a small number before log calculation to avoid NaN
    double tempError = 0.0;
    double error;
    double mus, fitA;
    double margin = (1.0 + ui->slider_ConcPercentChange->value() /200.0);

    for (int i=0; i<para->nWavel; i++)
    {
        mus = para->mus[i]*margin;  // set mus according to conc slider value

        x[i] = para->wavelArray[i];
        yMusp[i] = mus * (1.0 - para->g[i]);          //reduced scattering coefficient

        //Steve L Jacques,"Optical properties of biological tissues: a review" Phys. Med & Bio. 58(2013) R37-R61.
        //wavelength λ is normalized by a reference wavelength, 1000 nm
        fitA = para->fittedA *margin;
        yFit[i] = fitA *(para->fRay*pow(x[i]/1000.0, -4.0) + (1-para->fRay)*pow(x[i]/1000.0, -para->bMie));   // A(lambda/lambda1000)^-b
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

//Plot power law fit
void PlotData::PlotMuspPowerLaw(Ui_MainWindow *ui, QVector<double> x, QVector<double> yMusp, QVector<double> yFit)
{
    //Clear previous graphs
    ui->customPlot_MuspPowerLaw->removeGraph(1);
    ui->customPlot_MuspPowerLaw->removeGraph(0);
    ui->customPlot_MuspPowerLaw->legend->clearItems();

    //Main plot
    ui->customPlot_MuspPowerLaw->addGraph();
    ui->customPlot_MuspPowerLaw->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_MuspPowerLaw->graph(0)->setData(x, yMusp);
    ui->customPlot_MuspPowerLaw->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->customPlot_MuspPowerLaw->graph(0)->setPen( QPen( Qt::red, 1 ) );
    ui->customPlot_MuspPowerLaw->graph(0)->setName("μs' Data");
    ui->customPlot_MuspPowerLaw->legend->setVisible(true);
    if (ui->radioButton_LinearYAxis->isChecked())
        ui->customPlot_MuspPowerLaw->yAxis->setLabel("μs' (mmˉˡ)");
    if (ui->radioButton_LogYAxis->isChecked())
        ui->customPlot_MuspPowerLaw->yAxis->setLabel("Log (μs' (mmˉˡ))");

    //Power law fit plot
    ui->customPlot_MuspPowerLaw->addGraph();
    ui->customPlot_MuspPowerLaw->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->customPlot_MuspPowerLaw->graph(1)->setData(x,yFit);
    ui->customPlot_MuspPowerLaw->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->customPlot_MuspPowerLaw->graph(1)->setPen( QPen( Qt::blue, 1 ) );
    ui->customPlot_MuspPowerLaw->graph(1)->setName("Best Fit");
    ui->customPlot_MuspPowerLaw->legend->setVisible(true);
    ui->customPlot_MuspPowerLaw->axisRect()->insetLayout()->setInsetAlignment(0,Qt::AlignTop|Qt::AlignRight);
    ui->customPlot_MuspPowerLaw->graph()->rescaleAxes();
    ui->customPlot_MuspPowerLaw->replot();
    ui->customPlot_MuspPowerLaw->legend->setVisible(false);
}
