#ifndef CALCULATE_H
#define CALCULATE_H

#include <complex>
#include <QLabel>
#include "parameters.h"

class Calculate
{
public:
    Calculate();
    ~Calculate();

    double mWavel;                   // wavel: wavelength of the light in the medium (in microns)
    double mMu;                      // mu: cos angle
    double mK;                       // k: wavevector
    std::complex<double> mCS2;       // complex parallel component: far field
    std::complex<double> mCS1;       // complex perpendicular component: far field
    double mQSca;                    // scattering efficiency
    double mQExt;                    // extinction efficiency
    double mQBack;                   // backscattering efficiency
    double mMinTheta = 0;
    double mMaxTheta = M_PI;

    void DoSimulation(QLabel *progress, Parameters *para);
    void ComputeMuspAtRefWavel(Parameters *para);
    void CalculatePowerLawAutoFitSimple(Parameters *para);
    void CalculatePowerLawAutoFitComplex(Parameters *para);
    void DiameterRangeSetting(Parameters *para, unsigned int distIndex);
    void SetSphereRadiusAndRefIndex(Parameters *para, unsigned int index, bool flagVolOrConc);
    double CalculateForwardBackward(std::complex<double> *S1, std::complex<double> *S2,
                                    Parameters *para, unsigned int start, unsigned int end);
    double CalculateG(std::complex<double> *S1, std::complex<double> *S2, Parameters *para);
};

#endif // CALCULATE_H
