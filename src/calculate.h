#ifndef CALCULATE_H
#define CALCULATE_H

#include "ui_mainwindow.h"
#include <complex>
#include "parameters.h"

class calculate
{
public:
    calculate(void);
    ~calculate(void);

    double wavel;                   // wavel: wavelength of the light in the medium (in microns)
    double mu;                      // mu: cos angle
    double k;                       // k: wavevector
    std::complex<double> cS2;		// complex parallel component: far field
    std::complex<double> cS1;		// complex perpendicular component: far field
    double qSca;                    // scattering efficiency
    double qExt;                    // extinction efficiency
    double qBack;                   // backscattering efficiency
    double minTheta = 0;
    double maxTheta = M_PI;

    void DoSimulation(Ui_MainWindow *ui, parameters *para);
    void ComputeMuspAtRefWavel(parameters *para);
    void CalculatePowerLawAutoFitSimple(parameters *para);
    void CalculatePowerLawAutoFitComplex(parameters *para);
    void DiameterRangeSetting(parameters *para, unsigned int distIndex);
    void SetSphereRadiusAndRefIndex(parameters *para, unsigned int index, bool flagVolOrConc);
    double CalculateForwardBackward(std::complex<double> *S1, std::complex<double> *S2, parameters *para, unsigned int start, unsigned int end);
    double CalculateG(std::complex<double> *S1, std::complex<double> *S2, parameters *para);
};

#endif // CALCULATE_H
