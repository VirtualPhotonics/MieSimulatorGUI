#ifndef PLOTDATA_H
#define PLOTDATA_H

#include "ui_mainwindow.h"
#include "parameters.h"
#include "lib/qwt/qwt_polar_grid.h"

class PlotData
{
public:
    PlotData(void);
    ~PlotData(void);

    void ClearPlots(Ui_MainWindow *ui, parameters *para);

    void InitializeDistributionPlot(Ui_MainWindow *ui);
    void InitializePhaseFunctionPolarPlot(Ui_MainWindow *ui , parameters *para);
    void InitializePhaseFunctionLinearPlot(Ui_MainWindow *ui);
    void InitializeAllOtherPlots(Ui_MainWindow *ui);
    void AssignValuesDistributionPlot(Ui_MainWindow *ui, parameters* para);
    void AssignValuesPhaseFunctionPolarPlot(Ui_MainWindow *ui, parameters *para);
    void AssignValuesPhaseFunctionLinearPlot(Ui_MainWindow *ui, parameters *para);
    void AssignValuesS1S2Plot(Ui_MainWindow *ui, parameters *para);
    void AssignValuesAllOtherPlots(Ui_MainWindow *ui, parameters* para);
    void AssignValuesMuspPowerLawPlots(Ui_MainWindow *ui, parameters* para);
    void PlotDistribution(Ui_MainWindow *ui, parameters *para, QVector<double> x, QVector<double> yDist);    
    void PlotPhaseFunctionLinear(Ui_MainWindow *ui, QVector<double> x, QVector<double> yPara,
                                           QVector<double> yPerp, QVector<double> yAve);
    void PlotS1S2(Ui_MainWindow *ui, parameters *para, QVector<double> x, QVector<double> yS);
    void PlotScatteringCrossSection(Ui_MainWindow *ui, QVector<double> x, QVector<double> yCsca);
    void PlotExtinctionCrossSection(Ui_MainWindow *ui, QVector<double> x, QVector<double> yCext);
    void PlotBackscatteringCrossSection(Ui_MainWindow *ui, QVector<double> x, QVector<double> yCback);
    void PlotSizeParameter(Ui_MainWindow *ui, QVector<double> x, QVector<double> ySizePara);
    void PlotMus(Ui_MainWindow *ui,QVector<double> x, QVector<double> yMus);
    void PlotMusp(Ui_MainWindow *ui, QVector<double> x, QVector<double> yMusp);
    void PlotMuspCurveForPowerLawFit(Ui_MainWindow *ui, parameters *para, QVector<double> x, QVector<double> yMusp);
    void PlotG(Ui_MainWindow *ui, QVector<double> x, QVector<double> yG);
    void PlotForwardBackward(Ui_MainWindow *ui, QVector<double> x, QVector<double> yF, QVector<double> yB, bool legendFlag);
    void PlotMuspPowerLaw(Ui_MainWindow *ui, QVector<double> x, QVector<double> yMusp, QVector<double> yFit);

private:
    void PlotPhaseFunctionPolar(Ui_MainWindow *ui, parameters *para, QVector<double> theta, QVector<double> phaseFunction,
                                double minPolarPtheta, double maxPolarPtheta);
    QwtPolarGrid *mGrid;
};

//Class for polar data
class Data: public QwtSeriesData<QwtPointPolar>
{
public:
    Data( const QwtInterval &radialInterval,
          const QwtInterval &azimuthInterval,
          size_t size,
          QVector<double> theta,
          QVector<double> phaseFunction):
              d_radialInterval( radialInterval ),
              d_azimuthInterval( azimuthInterval ),
              d_size( size )
    {
        d_radial = phaseFunction;
        d_azimuth = theta;
    }
    virtual size_t size() const
    {
        return d_size;
    }

protected:
    QwtInterval d_radialInterval;
    QwtInterval d_azimuthInterval;
    size_t d_size;
    QVector<double> d_radial;
    QVector<double> d_azimuth;
};

//Class for polar plot
class Polar: public Data
{
public:
    Polar( const QwtInterval &radialInterval,
            const QwtInterval &azimuthInterval,
           size_t size, QVector<double> &theta,
           QVector<double> &phaseFunction): Data(
                              radialInterval, azimuthInterval,
                              size, theta, phaseFunction)
    {
    }
    virtual QwtPointPolar sample( size_t i ) const
    {
        int j = static_cast<int>(i);
        return QwtPointPolar( d_azimuth[j],d_radial[j]);
    }

    virtual QRectF boundingRect() const
    {
        if ( d_boundingRect.width() < 0.0 )
            d_boundingRect = qwtBoundingRect( *this );
        return d_boundingRect;
    }
};

#endif // PLOTDATA_H
