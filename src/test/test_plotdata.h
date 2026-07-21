#ifndef TEST_PLOTDATA_H
#define TEST_PLOTDATA_H

#include <QObject>
#include "parameters.h"
#include "dialog/plotdata.h"

class TestPlotData : public QObject
{
    Q_OBJECT

public:
    TestPlotData() = default;;
    ~TestPlotData() = default;;

private slots:
    void init();
    void cleanup();

    void test_FindMinLogPolarPlot_singleWavelength();
    void test_FindMinLogPolarPlot_minimumInLaterWavelength();
    void test_RearrangePhaseFunctionData_negPosDegrees_linear();
    void test_RearrangePhaseFunctionData_zeroToTwoPiRadians_log();

private:
    Parameters *mPara;
    PlotData *mPlotData;
};

#endif // TEST_PLOTDATA_H