#ifndef UTILITIES_H
#define UTILITIES_H

#include <QTime>
#include <QCoreApplication>
#include <complex>

class Utilities
{
public:
    Utilities() = default;
    ~Utilities() = default;

    void Delay() const;
    double ComplexAbs(std::complex<double> c) const;
    double ComplexAbsSquared(std::complex<double> c) const;
    double SimpsonsWeight (unsigned int i, unsigned int n) const;
    double NiceStep(double range, int initialCircles) const;
    double FindMinMax(const QVector<double>& yPara, const QVector<double>& yPerp,
                      const QVector<double>& yAve, bool flagMinMax) const;
};

#endif // UTILITIES_H