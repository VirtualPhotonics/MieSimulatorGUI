#ifndef PLOTDATA_H
#define PLOTDATA_H

#include "ui_mainwindow.h"
#include "parameters.h"

class PlotData
{
public:
    PlotData(void);
    ~PlotData(void);

    void ClearPlots(Ui_MainWindow *ui);
    void InitialSetupDistributionPlot(Ui_MainWindow *ui);
    void InitialSetupPhaseFunctionLinearPlot(Ui_MainWindow *ui);
    void InitialSetupPhaseFunctionPolarPlot(Ui_MainWindow *ui);
    void InitialSetupS1S2Plot(Ui_MainWindow *ui);
    void InitialSetupMuspPowerLawFit(Ui_MainWindow *ui);
    void InitialSetupOtherPlots(Ui_MainWindow *ui);
    void SetupPolarPlotForData(Ui_MainWindow *ui, Parameters *para );
    void AssignValuesDistributionPlot(Ui_MainWindow *ui, Parameters* para);
    void AssignValuesPhaseFunctionLinearPlot(Ui_MainWindow *ui, Parameters *para);
    void AssignValuesPhaseFunctionPolarPlot(Ui_MainWindow *ui, Parameters *para);
    void AssignValuesS1S2Plot(Ui_MainWindow *ui, Parameters *para);
    void AssignValuesMuspPowerLawPlots(Ui_MainWindow *ui, Parameters* para);
    void AssignValuesOtherPlots(Ui_MainWindow *ui, Parameters* para);
    double FindMinLogPolarPlot(Parameters *para);

private:
    void PlotDistribution(Ui_MainWindow *ui, Parameters *para, QVector<double> x, QVector<double> yDist);
    void PlotPhaseFunctionLinear(Ui_MainWindow *ui, QVector<double> x, QVector<double> yPara,
                                 QVector<double> yPerp, QVector<double> yAve);
    void PlotPhaseFunctionPolar(Ui_MainWindow *ui, Parameters *para, QVector<double> thetaNorth,
                                QVector<double> thetaSouth, QVector<double> phaseFunction);
    void PlotS1S2(Ui_MainWindow *ui, Parameters *para, QVector<double> x, QVector<double> yS);
    void PlotMuspCurveForPowerLaw(Ui_MainWindow *ui, QVector<double> x, QVector<double> yMusp);
    void PlotMuspPowerLaw(Ui_MainWindow *ui, QVector<double> x, QVector<double> yMusp, QVector<double> yFit);
    void PlotScatteringCrossSection(Ui_MainWindow *ui, QVector<double> x, QVector<double> yCsca);
    void PlotExtinctionCrossSection(Ui_MainWindow *ui, QVector<double> x, QVector<double> yCext);
    void PlotBackscatteringCrossSection(Ui_MainWindow *ui, QVector<double> x, QVector<double> yCback);
    void PlotSizeParameter(Ui_MainWindow *ui, QVector<double> x, QVector<double> ySizePara);
    void PlotMus(Ui_MainWindow *ui,QVector<double> x, QVector<double> yMus);
    void PlotG(Ui_MainWindow *ui, QVector<double> x, QVector<double> yG);
    void PlotMusp(Ui_MainWindow *ui, QVector<double> x, QVector<double> yMusp);
    void PlotForwardBackward(Ui_MainWindow *ui, QVector<double> x, QVector<double> yF, QVector<double> yB, bool legendFlag);

    void DrawPolarPlotGrid(QCustomPlot *customPlot, bool flagLinearLog);
    void CreateCircularGrid(QCustomPlot* customPlot, bool flagLinearLog);
    void CreateRadialGrid(QCustomPlot *customPlot, bool flagLinearLog);    
    void HideCartesianAxes(QCustomPlot* customPlot);
    void PlotSingleGraph(QCustomPlot* customPlot, const QVector<double>& x, const QVector<double>& y,
                         QColor color, const QString& name, int graphIndex, int sizeCircle);
    void RemoveGraphs(QCustomPlot *customPlot);
    void RemoveLegends(QCustomPlot *customPlot);
    void RemoveItems(QCustomPlot *customPlot);
    void RemovePlotables (QCustomPlot *customPlot);

    double mPolarMinRadius = 0;
    double mPolarMaxRadius = 1.0;
};

#endif // PLOTDATA_H
