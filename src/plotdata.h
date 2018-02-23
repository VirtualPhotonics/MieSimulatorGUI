#ifndef PLOTDATA_H
#define PLOTDATA_H

#include "ui_mainwindow.h"
#include "qcustomplot.h"
#include "parameters.h"
#include "calculate.h"
#include "qwt/qwt_polar_grid.h"
#include "qwt/qwt_polar_curve.h"
#include "qwt/qwt_scale_engine.h"

class PlotData
{
public:
    PlotData(void);
    ~PlotData(void);

    void ClearPlots(Ui_MainWindow *ui, parameters *para);

    void InitializeDistributionPlot(Ui_MainWindow *ui);
    void InitializePolarPlot(Ui_MainWindow *ui , parameters *para);
    void InitializeAllOtherPlots(Ui_MainWindow *ui);
    void AssignValuesDistributionPlot(Ui_MainWindow *ui, parameters* para);
    void AssignValuesPolarPlot(Ui_MainWindow *ui, parameters *para);
    void AssignValuesS1S2Plot(Ui_MainWindow *ui, parameters *para);
    void AssignValuesAllOtherPlots(Ui_MainWindow *ui, parameters* para);
    void AssignValuesMuspPowerLawPlots(Ui_MainWindow *ui, parameters* para);

    void PlotDistribution(Ui_MainWindow *ui, parameters *para, QVector<double> x, QVector<double> yDist);
    void PlotPolar(Ui_MainWindow *ui, parameters *para, QVector<double> theta, QVector<double> phaseFunction);
    void PlotScatCross(Ui_MainWindow *ui, QVector<double> x, QVector<double> yScat);
    void PlotMus(Ui_MainWindow *ui,QVector<double> x, QVector<double> yMus);
    void PlotMusp(Ui_MainWindow *ui, QVector<double> x, QVector<double> yMusp);
    void PlotMuspCurveForPowerLawFit(Ui_MainWindow *ui, parameters *para, QVector<double> x, QVector<double> yMusp);
    void PlotG(Ui_MainWindow *ui, QVector<double> x, QVector<double> yG);
    void PlotForwardBackward(Ui_MainWindow *ui, QVector<double> x, QVector<double> yF, QVector<double> yB, bool legendFlag);
    void PlotMuspPowerLaw(Ui_MainWindow *ui, QVector<double> x, QVector<double> yMusp, QVector<double> yFit);
    void PlotS1S2(Ui_MainWindow *ui, parameters *para, QVector<double> x, QVector<double> yS);

private:
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
        return QwtPointPolar( d_azimuth[i],d_radial[i]);
    }

    virtual QRectF boundingRect() const
    {
        if ( d_boundingRect.width() < 0.0 )
            d_boundingRect = qwtBoundingRect( *this );
        return d_boundingRect;
    }
};

#endif // PLOTDATA_H
