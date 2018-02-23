#ifndef MIESIMULATION_H
#define MIESIMULATION_H

#include <complex>
#include "utilities.h"

class MieSimulation
{
public:
    MieSimulation(void);
    ~MieSimulation(void);

    void FarFieldSolutionForRealRefIndex(std::complex<double> *cS1, std::complex<double> *cS2, double *qSca, double xPara, double relRef, double mu);
    void FarFieldSolutionForComplexRefIndex(std::complex<double> *cS1, std::complex<double> *cS2, double *qSca, double xPara, std::complex<double> cRelRef, double mu);
};

#endif // MIESIMULATION_H
