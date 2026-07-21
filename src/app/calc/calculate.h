#ifndef CALCULATE_H
#define CALCULATE_H

#include <complex>
#include <vector>
#include <QLabel>
#include "parameters.h"
#include "calc/miesimulation.h"

// Precomputed angle grid (cos/sin/Simpson weights) for a given nTheta.
struct ThetaGrid
{
    std::vector<double> cosTheta;        // cos(theta[t]), size nTheta
    std::vector<double> twoPiSinTheta;   // 2*pi*sin(theta[t]), size nTheta
    std::vector<double> weightFull;      // Simpson weight over the full [0, nTheta) range
    std::vector<double> weightForward;   // Simpson weight for the forward segment [0, forwardLen)
    std::vector<double> weightBackward;  // Simpson weight for the backward segment [0, backwardLen)
};

class Calculate
{
public:
    Calculate() = default;
    ~Calculate() = default;

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
    void BuildThetaGrid(Parameters *para, ThetaGrid &grid);
    double CalculateForwardBackward(std::complex<double> *S1, std::complex<double> *S2,
                                    Parameters *para, unsigned int start, unsigned int end);
    double CalculateForwardBackward(std::complex<double> *S1, std::complex<double> *S2,
                                    unsigned int start, unsigned int end,
                                    const ThetaGrid &grid, const std::vector<double> &weight);
    double CalculateG(std::complex<double> *S1, std::complex<double> *S2, Parameters *para);
    double CalculateG(std::complex<double> *S1, std::complex<double> *S2, Parameters *para, const ThetaGrid &grid);
    bool CheckIndependentScattering(Parameters *para, double &clearanceToWavelength,
                                    double &sizeParameter, double &volFraction,
                                    double &wavelength, double &clearance,
                                    QString &strRegime, bool flagVolFlag);

private:
    void ComputeMieForSphere(Parameters *para, unsigned int r, double wavelength,
                             const ThetaGrid &grid,
                             MieSimulation &mieSim, MieCoefficients &mieCoeff,
                             std::complex<double> *curS1, std::complex<double> *curS2,
                             double &xPara, double &piRadiusSquared);

    double PreparePowerLawFitData(Parameters *para, std::vector<double> &xRatio,
                                  std::vector<double> &y);
};

#endif // CALCULATE_H