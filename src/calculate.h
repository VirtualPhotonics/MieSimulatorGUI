#ifndef CALCULATE_H
#define CALCULATE_H

#include "ui_mainwindow.h"
#include <complex>
#include <qmath.h>
#include "parameters.h"
#include "miesimulation.h"
#include "utilities.h"


class calculate
{
public:
    calculate(void);
    ~calculate(void);    

    void DoSimulation(Ui_MainWindow *ui, parameters *para);
    void ComputeMuspAtRefWavel(parameters *para);
    void CalculatePowerLawAutoFitSimple(parameters *para);
    void CalculatePowerLawAutoFitComplex(parameters *para);
    void DiameterRangeSetting(parameters *para, unsigned int distIndex);
    void SetSphereRadiusAndRefIndex(parameters *para, unsigned int index, bool flagVolOrConc);
    double CalculateForwardBackward(std::complex<double> *S1, std::complex<double> *S2, parameters *para, unsigned int start, unsigned int end);
    double CalculateG(std::complex<double> *S1, std::complex<double> *S2, parameters *para);

private:
    double _wavel;		            // wavel: wavelength of the light in the medium (in microns)
    double _mu;                     // mu: cos angle
    double _k;                      // k: wavevector
    std::complex<double> _cS2;		// complex parallel component: far field
    std::complex<double> _cS1;		// complex perpendicular component: far field
    double _qSca;                   // scattering efficiency
    double _qExt;                   // extinction efficiency
    double _qBack;                  // backscattering efficiency
};

#endif // CALCULATE_H
