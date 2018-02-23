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
void PlotData::InitializePolarPlot(Ui_MainWindow *ui, parameters *para )
{
//    SetSliderWL(ui);

    //Polar plot: phase function
    ui->qwtpolarplot_PhaseFunction->setAutoReplot( false );
    ui->qwtpolarplot_PhaseFunction->setPlotBackground( Qt::white );
    ui->qwtpolarplot_PhaseFunction->setScale( QwtPolar::Azimuth,0, 360, 360 / 12 );
    ui->qwtpolarplot_PhaseFunction->setScale( QwtPolar::Radius, 0, 1);
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
    mGrid->attach( ui->qwtpolarplot_PhaseFunction );
    para->polarCurve->detach();
    ui->qwtpolarplot_PhaseFunction->replot();
}

//Initialize all other plots
void PlotData::InitializeAllOtherPlots(Ui_MainWindow *ui)
{
    double minX = ui->lineEdit_StartWL->text().toDouble();
    double maxX = ui->lineEdit_EndWL->text().toDouble();

    //Plot S1S2
    ui->customPlot_S1S2->xAxis->setLabel("Angle (deg.)");
    ui->customPlot_S1S2->xAxis->setRange(0, 180);
    ui->customPlot_S1S2->yAxis->setLabel("S1");
    ui->customPlot_S1S2->clearGraphs();
    ui->customPlot_S1S2->replot();
    ui->customPlot_S1S2->addGraph();   //dummy

    //Plot ScatCross Section
    ui->customPlot_ScatCross->xAxis->setLabel("Wavelength (nm)");
    ui->customPlot_ScatCross->xAxis->setRange(minX, maxX);
    ui->customPlot_ScatCross->yAxis->setLabel("Scattering Cross Section (μm²)");
    ui->customPlot_ScatCross->clearGraphs();
    ui->customPlot_ScatCross->replot();
    ui->customPlot_ScatCross->addGraph();   //dummy

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
    InitializeAllOtherPlots(ui);
    InitializePolarPlot(ui, para);
}

//Assign values for distribution plot
void PlotData::AssignValuesDistributionPlot(Ui_MainWindow *ui, parameters* para)
{
    QVector<double> xDist(para->nRadius), yDist(para->nRadius);
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
        ui->tabWidget_Ns->setTabText(0,"Concentration (Number Density: Ns)");

    }
    if (ui->radioButton_PolyDisperse->isChecked())
    {
        double shift = (para->maxRadius - para->minRadius)/(para->nRadius);
        minX = 2.0 * para->minRadius - shift;
        maxX = 2.0 * para->maxRadius + shift;
        bar->setWidth((para->maxRadius - para->minRadius)/para->nRadius);
        ui->customPlot_Distribution->xAxis->setRange(minX, maxX);
        ui->tabWidget_Ns->setTabText(0,"Scatterer Distribution and Concentration (Ns)");
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
void PlotData::AssignValuesPolarPlot(Ui_MainWindow *ui, parameters *para)
{
    QVector<double> phaseFunction(2*para->nTheta-1), theta(2*para->nTheta-1);
    int indexWL = ui->slider_PF_WL->value();
    ui->label_CurrentWL_PF->setText(QString::number(para->startWavel + indexWL*para->stepWavel));

    //Set minimum and maximum values
    para->maxPolarPtheta = 0;
    para->minPolarPtheta = 1e100;
    for (int i = 0; i < para->nWavel; i++)
    {
        for (int j = 0; j < para->nTheta; j++)
        {
            if (ui->radioButton_MonoDisperse->isChecked())
            {
                if (para->phaseFunctionAve[i][j]>para->maxPolarPtheta)
                    para->maxPolarPtheta = para->phaseFunctionAve[i][j];
                if (para->phaseFunctionPara[i][j]>para->maxPolarPtheta)
                    para->maxPolarPtheta = para->phaseFunctionPara[i][j];
                if (para->phaseFunctionPerp[i][j]>para->maxPolarPtheta)
                    para->maxPolarPtheta = para->phaseFunctionPerp[i][j];
                if (para->phaseFunctionAve[i][j]<para->minPolarPtheta)
                    para->minPolarPtheta = para->phaseFunctionAve[i][j];
                if (para->phaseFunctionPara[i][j]<para->minPolarPtheta)
                    para->minPolarPtheta = para->phaseFunctionPara[i][j];
                if (para->phaseFunctionPerp[i][j]<para->minPolarPtheta)
                    para->minPolarPtheta = para->phaseFunctionPerp[i][j];
            }

            if (ui->radioButton_PolyDisperse->isChecked())
            {
                if (para->phaseFunctionAve[i][j]>para->maxPolarPtheta)
                    para->maxPolarPtheta = para->phaseFunctionAve[i][j];
                if (para->phaseFunctionAve[i][j]<para->minPolarPtheta)
                    para->minPolarPtheta = para->phaseFunctionAve[i][j];
            }
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
    PlotPolar(ui, para, theta, phaseFunction);
}

//Plot phase function polar plot
void PlotData::PlotPolar(Ui_MainWindow *ui, parameters *para, QVector<double> theta, QVector<double> phaseFunction)
{
    //Polar plot: phase function
    QwtInterval radial( 0, para->maxPolarPtheta );
    QwtInterval azimuth( 0.0, 360.0 );
    ui->qwtpolarplot_PhaseFunction->setScale(QwtPolar::Radius, 0, para->maxPolarPtheta);
    ui->qwtpolarplot_PhaseFunction->setScale(QwtPolar::Azimuth, 0, 360, 360/12);

    para->polarCurve->detach();
    ui->qwtpolarplot_PhaseFunction->replot();
    para->polarCurve->setStyle( QwtPolarCurve::Lines );
    para->polarCurve->setPen( QPen( Qt::red, 2 ) );
    para->polarCurve->setData(new Polar(radial,azimuth,2*para->nTheta-1,theta, phaseFunction));
    para->polarCurve->attach(ui->qwtpolarplot_PhaseFunction);
    if (ui->radioButton_PhaseLinear->isChecked())
        ui->qwtpolarplot_PhaseFunction->setScaleEngine( QwtPolar::Radius, new QwtLinearScaleEngine() );
    if (ui->radioButton_PhaseLog->isChecked())
    {
        ui->qwtpolarplot_PhaseFunction->setScale(QwtPolar::Radius, para->minPolarPtheta, para->maxPolarPtheta);
        ui->qwtpolarplot_PhaseFunction->setScaleEngine( QwtPolar::Radius, new QwtLogScaleEngine() );
    }
    ui->qwtpolarplot_PhaseFunction->replot();
    if (para->nWavel==1)
        ui->slider_PF_WL->setDisabled(true);
    else
        ui->slider_PF_WL->setDisabled(false);
}

//Assign values for Scattering cross section, g and Musp plots
void PlotData::AssignValuesAllOtherPlots(Ui_MainWindow *ui, parameters* para)
{
    QVector<double> x(para->nWavel);
    QVector<double> yScat(para->nWavel), yG(para->nWavel),yMus(para->nWavel), yMusp(para->nWavel);
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
        yScat[i] = para->scatCross[i];    //scattering cross section
        yG[i] = para->g[i];               //g - average cosine of phase function
        yMus[i] = mus;          //scattering coefficient
        yMusp[i] = mus * (1.0 - para->g[i]);          //reduced scattering coefficient
        yF[i] = para->forward[i];
        yB[i] = para->backward[i];

        if (ui->radioButton_LogYAxis->isChecked())
        {
            yScat[i] = log10(para->scatCross[i] + tiny);    //scattering cross section
            yG[i] = log10(para->g[i] +tiny);                //g - average cosine of phase function
            yMus[i] = log10(mus + tiny);                    //scattering coefficient
            yMusp[i] = log10((mus * (1.0 - para->g[i])) + tiny);    //reduced scattering coefficient
            yF[i] = log10(para->forward[i] + tiny);
            yB[i] = log10(para->backward[i] + tiny);
            fbLimit = 2;   //100%
        }
    }    
    PlotG(ui, x, yG);
    if (yF[fbLegendCheckLocation]>fbLimit)
        fbLegnedFlag = true;
    else
        fbLegnedFlag = false;
    PlotForwardBackward(ui, x, yF, yB, fbLegnedFlag);

    PlotScatCross(ui, x, yScat);
    PlotMus(ui, x, yMus);
    PlotMusp(ui, x, yMusp);
    PlotMuspCurveForPowerLawFit(ui, para, x, yMusp);    
}

//Plot ScatCross plot
void PlotData::PlotScatCross(Ui_MainWindow *ui, QVector<double> x, QVector<double> yScat)
{
   //Clear previous graph
   ui->customPlot_ScatCross->removeGraph(0);
   //Add new graph
   ui->customPlot_ScatCross->addGraph();
   ui->customPlot_ScatCross->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
   ui->customPlot_ScatCross->graph(0)->setData(x, yScat);
   ui->customPlot_ScatCross->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
   ui->customPlot_ScatCross->graph(0)->setPen( QPen( Qt::red, 1 ) );
   if (ui->radioButton_LinearYAxis->isChecked())
       ui->customPlot_ScatCross->yAxis->setLabel("Scattering Cross Section (μm²)");
   if (ui->radioButton_LogYAxis->isChecked())
       ui->customPlot_ScatCross->yAxis->setLabel("Log (Scattering Cross Section (μm²) )");
   ui->customPlot_ScatCross->graph(0)->rescaleAxes();
   ui->customPlot_ScatCross->replot();
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

//Assign values for S1/S2 plot
void PlotData::AssignValuesS1S2Plot(Ui_MainWindow *ui, parameters *para)
{
    QVector<double> S(para->nTheta);
    QVector<double> theta(para->nTheta);
    std::complex<double> *tempS;

    int indexWL = ui->slider_S1S2_WL->value();
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
        ui->slider_S1S2_WL->setDisabled(true);
    else
        ui->slider_S1S2_WL->setDisabled(false);
}
