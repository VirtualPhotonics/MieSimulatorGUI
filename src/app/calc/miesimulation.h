#ifndef MIESIMULATION_H
#define MIESIMULATION_H

#include <complex>

class MieSimulation
{
public:
    MieSimulation();
    ~MieSimulation();

    void FarFieldSolutionForRealRefIndex(std::complex<double> *cS1, std::complex<double> *cS2, double *qSca,
                                         double *qExt, double *qBack, double xPara, double relRef, double mu);
    void FarFieldSolutionForComplexRefIndex(std::complex<double> *cS1, std::complex<double> *cS2, double *qSca,
                                            double *qExt, double *qBack, double xPara,
                                            std::complex<double> cRelRef, double mu);
};

#endif // MIESIMULATION_H
