#ifndef UTILITIES_H
#define UTILITIES_H

#include <QTime>
#include <QCoreApplication>
#include <complex>

class Utilities
{
public:
    Utilities();
    ~Utilities();

    void Delay();
    double ComplexAbs(std::complex<double> c);
    double ComplexAbsSquared(std::complex<double> c);
    double SimpsonsWeight (unsigned int i, unsigned int n);
    double NiceStep(double range, int initialCircles);
    double FindMinMax(const QVector<double>& yPara, const QVector<double>& yPerp,
                      const QVector<double>& yAve, bool flagMinMax);
};

#endif // UTILITIES_H
