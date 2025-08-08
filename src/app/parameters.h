#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <complex>
#include <qmath.h>
#include "lib/qwt/qwt_polar_curve.h"
#include <QRadioButton>
#include <QComboBox>

class parameters
{
public:
    parameters(void);
    ~parameters(void);

    bool CheckCommonParameters(QRadioButton *radioButton_MonoDisperse,
                               QRadioButton *radioButton_NumDen,
                               QRadioButton *radioButton_VolFrac);
    bool CheckDistributionParameters(QComboBox *comboBox_Distribution);


    double *radArray = nullptr;         // radArray: radius Array
    double *numDensityArray = nullptr;  // numDensityArray: Number Density array
    double meanRadius = 0.1;            // meanRadius: "mean radius" in Poly disperse and "radius" in mono disperse
    double minRadius = 0.0005;          // minRadius: minimum radius
    double maxRadius = 25;              // maxRadius: maximum radius
    double stdDev = 0.25;               // stdDev: Std. deviation
    unsigned int nRadius = 31;          // nRadius: number of radii

    double *scatRefRealArray = nullptr; // scatRefRealArray: refractive index of scatterer -Real  data array
    double *scatRefImagArray = nullptr; // scatRefImagArray: refractive index of scatterer -Imag  data array
    double scatRefReal = 1.377;         // scatRefReal: refractive index of scatterer - Real
    double scatRefImag = 0.0;           // scatRefImg: refractive index of scatterer - Imaginary ((n-ik): negative sign convention as Scott Prahl's Mie calculator)
    double medRef = 1.333;              // medRef: refractive index of medium
    double volFraction = 0.1;           // volFraction: Volume fraction of sphere volume
    double sphNumDensity = 1e8;         // sphNumDensity: Sphere concentration/volume (Number Density)

    double startWavel = 600;            // startWavel: starting wavelength
    double endWavel = 1000;             // endWavel: starting wavelength
    double stepWavel = 10;              // stepWavel: starting wavelength
    double *wavelArray = nullptr;       // wavel: wavelength array
    unsigned int nWavel = 41;           // nWavel: number of wavelength

    unsigned int nTheta = 361;          // nTheta: number of angles for phase function, S1/S2

    double **phaseFunctionAve = nullptr;    // phaseFunctionAve: Average phase function
    double **phaseFunctionPara = nullptr;   // phaseFunctionPara: parallel phase function
    double **phaseFunctionPerp = nullptr;   // phaseFunctionPerp: perpendicularphase function
    double *curPolarAng = nullptr;          // curPolarAng: current polar angles for phase function
    double *curPhaseFunc = nullptr;         // curPhaseFunc: current Phase function
    double *cSca = nullptr;                 // cSca: Scattering cross section
    double *cExt = nullptr;                 // cExt: Extinction cross section
    double *cBack = nullptr;                // cBack: Backscattering cross section
    double *SizePara = nullptr;             // SizePara: Size Parameter
    double *mus = nullptr;                  // mus: scattering coefficient
    double *g = nullptr;                    // g: average cosine of phase function
    std::complex<double> **S1 = nullptr;    // S1: amplitude matrix component S1
    std::complex<double> **S2 = nullptr;    // S2: amplitude matrix component S2
    double *forward = nullptr;              // Total intensity from forward hemisphere
    double *backward = nullptr;             // Total intensity from backward hemisphere

    double fRay = 1.0;                         //fRay: fitting parameter fRay
    double bMie = 4.0;                         //bMie: fitting parameter bMie
    double muspAtRefWavel[6] = {1,1,1,1,1,1};  //muspAtRefWavel:  musp value array at lambda = 500, 600, 700, 800, 900 and 1000nm
    unsigned int refWavelIdx = 0;              //refWavelIdx: current index of refWavel
    double muspFittingError = 0;               //muspfittingError: fitting error
    double refWavel = 1000;                    //refWavel: reference lambda: 500, 600, 700, 800, 900 and 1000nm
    bool fittingComplex = true;                //fittingComplex; False: Simple,  True: Ccomplex

    QwtPolarCurve *polarCurve = nullptr;       //polarCurve: polarCurve data
};

#endif // PARAMETERS_H
