#ifndef MIESIMULATION_H
#define MIESIMULATION_H

#include <complex>
#include <vector>

struct MieCoefficients
{
    std::vector<std::complex<double>> an;
    std::vector<std::complex<double>> bn;
    unsigned int nStop = 0;
    double qSca = 0.0;
    double qExt = 0.0;
    double qBack = 0.0;
};

class MieSimulation
{
public:
    MieSimulation() = default;
    ~MieSimulation() = default;

    // Angle-independent stage
    void ComputeCoefficientsForRealRefIndex(MieCoefficients &coeff, double xPara, double relRef);
    void ComputeCoefficientsForComplexRefIndex(MieCoefficients &coeff, double xPara, std::complex<double> cRelRef);

    // Angle-dependent stage
    void ComputeAngularS1S2(std::complex<double> *cS1, std::complex<double> *cS2,
                            const MieCoefficients &coeff, double mu);

    void FarFieldSolutionForRealRefIndex(std::complex<double> *cS1, std::complex<double> *cS2, double *qSca,
                                         double *qExt, double *qBack, double xPara, double relRef, double mu);
    void FarFieldSolutionForComplexRefIndex(std::complex<double> *cS1, std::complex<double> *cS2, double *qSca,
                                            double *qExt, double *qBack, double xPara,
                                            std::complex<double> cRelRef, double mu);
};

#endif // MIESIMULATION_H