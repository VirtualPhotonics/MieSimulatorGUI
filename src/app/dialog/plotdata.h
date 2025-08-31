#ifndef PLOTDATA_H
#define PLOTDATA_H

#include "lib/qcustomplot.h"
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
    void PlotPhaseFunctionPolar(Ui_MainWindow *ui, QVector<double> theta, QVector<double> yPara,
                                QVector<double> yPerp, QVector<double> yAve, int totalSize);
    void PlotS1S2(Ui_MainWindow *ui, Parameters *para, QVector<double> x, QVector<double> yS1, QVector<double> yS2);
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

    void PlotSingleGraph(QCustomPlot* customPlot, const QVector<double>& x, const QVector<double>& y,
                         QColor color, const QString& name, int graphIndex, int sizeCircle);
    void PlotSingleCurve(QCustomPlot* customPlot, const QVector<double>& xData, const QVector<double>& yData,
                         const QColor& color, const QString& name, int totalSize);
    void InitialSetupPlot(QCustomPlot *customPlot, const QString& xLabel, const QString& yLabel, double minX, double maxX);
    void DrawPolarPlotGrid(QCustomPlot *customPlot, bool flagLinearLog);
    void CreateCircularGrid(QCustomPlot* customPlot, bool flagLinearLog);
    void CreateRadialGrid(QCustomPlot *customPlot, bool flagLinearLog);
    void HideCartesianAxes(QCustomPlot* customPlot);
    void RemoveGraphs(QCustomPlot *customPlot);
    void RemoveLegends(QCustomPlot *customPlot);
    void RemoveItems(QCustomPlot *customPlot);
    void RemoveEllipseLineTextItems(QCustomPlot *customPlot);
    void RemovePlotables (QCustomPlot *customPlot);
    void RearrangePhaseFunctionData(Parameters *para, QVector<double> &theta, QVector<double> &phaseFuncPara,
                                    QVector<double> &phaseFuncPerp, QVector<double> &phaseFuncAve,
                                    int indexWL, bool flagLinearLog, bool flagThetaNegPosOrPos);

    //Variables to correctly draw the polar plot
    double mPolarMinRadius = 0;
    double mPolarMaxRadius = 1.0;
};

#endif // PLOTDATA_H
