#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <complex>
#include "qwt/qwt_polar_curve.h"

class parameters
{
public:
    parameters(void);
    ~parameters(void);

    double *radArray;           // radArray: radius Array
    double *numDensityArray;    // numDensityArray: Number Density array
    double meanRadius;          // meanRadius: "mean radius" in Poly disperse and "radius" in mono disperse
    double minRadius;           // minRadius: minimum radius
    double maxRadius;           // maxRadius: maximum radius
    double stdDev;              // stdDev: Std. deviation
    unsigned int nRadius;       // nRadius: number of radii

    double *scatRefRealArray;   // scatRefRealArray: refractive index of scatterer -Real  data array
    double *scatRefImagArray;   // scatRefImagArray: refractive index of scatterer -Imag  data array
    double scatRefReal;         // scatRefReal: refractive index of scatterer - Real
    double scatRefImag;         // scatRefImg: refractive index of scatterer - Imaginary ((n-ik): negative sign convention as Scott Prahl's Mie calculator)
    double medRef;              // medRef: refractive index of medium
    double volFraction;         // volFraction: Volume fraction of sphere volume
    double sphNumDensity;       // sphNumDensity: Sphere concentration/volume (Number Density)

    double startWavel;          // startWavel: starting wavelength
    double endWavel;            // endWavel: starting wavelength
    double stepWavel;           // stepWavel: starting wavelength
    double *wavelArray;         // wavel: wavelength array
    unsigned int nWavel;        // nWavel: number of wavelength

    double minTheta;            // minTheta: minimum angle
    double maxTheta;            // maxTheta: maximum angle
    double stepTheta;           // stepTheta: maximum angle
    unsigned int nTheta;        // nTheta: number of angles

    double **phaseFunctionAve;  // phaseFunctionAve: Average phase function
    double **phaseFunctionPara; // phaseFunctionPara: parallel phase function
    double **phaseFunctionPerp; // phaseFunctionPerp: perpendicularphase function
    double *curPolarAng;        // curPolarAng: current polar angles for phase function
    double *curPhaseFunc;       // curPhaseFunc: current Phase function
    double *cSca;               // cSca: Scattering cross section
    double *cExt;               // cExt: Extinction cross section
    double *cBack;              // cBack: Backscattering cross section
    double *SizePara;           // SizePara: Size Parameter
    double *mus;                // mus: scattering coefficient
    double *g;                  // g: average cosine of phase function
    std::complex<double> **S1;  // S1: amplitude matrix component S1
    std::complex<double> **S2;  // S2: amplitude matrix component S2
    double *forward;            // Total intensity from forward hemisphere
    double *backward;           // Total intensity from backward hemisphere

    double minPolarPtheta;      // minR: Minimum radius for polar plot
    double maxPolarPtheta;      // rMax: Maximum radius for polar plot

    double fRay;                //fRay: fitting parameter fRay
    double bMie;                //bMie: fitting parameter bMie
    double muspAtRefWavel[6];   //muspAtRefWavel:  musp value array at lambda = 500, 600, 700, 800, 900 and 1000nm
    unsigned int refWavelIdx;   //refWavelIdx: current index of refWavel
    double muspFittingError;    //muspfittingError: fitting error
    double refWavel;            //refWavel: reference lambda: 500, 600, 700, 800, 900 and 1000nm
    bool fittingComplex;        //fittingComplex; False: Simple,  True: Ccomplex

    QwtPolarCurve *polarCurve;  //polarCurve: polarCurve data
};

#endif // PARAMETERS_H
