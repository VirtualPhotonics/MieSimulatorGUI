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

    double wavel;		// wavel: wavelength of the light in the medium (in microns)
    double mu;          // mu: cos angle
    double k;           // k: wavevector
    std::complex<double> cS2;		// complex parallel component: far field
    std::complex<double> cS1;		// complex perpendicular component: far field
    double qSca;                    // scattering efficiency
    double qExt;                    // extinction efficiency
    double qBack;                   // backscattering efficiency

    void DoSimulation(Ui_MainWindow *ui, parameters *para, double* mus1000, double *g1000);
    void CalculatePowerLawAutoFit(parameters *para);
    void DiameterRangeSetting(parameters *para, int distIndex);
    void SetSphereRadiusAndRefIndex(parameters *para, int index, bool flagVolOrConc);
    double CalculateForwardBackward(std::complex<double> *S1, std::complex<double> *S2, parameters *para, int start, int end);
    double CalculateG(std::complex<double> *S1, std::complex<double> *S2, parameters *para);
};

#endif // CALCULATE_H
