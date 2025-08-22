//**********************************************************************
//** Some utility functions that are commonly used are listed here.
//**********************************************************************

#include "calc/utilities.h"

Utilities::Utilities()
{
}

Utilities::~Utilities()
{
}

//Delay function for display update
void Utilities::Delay()
{
    QTime dieTime= QTime::currentTime().addMSecs(1);
    while( QTime::currentTime() < dieTime )
    QCoreApplication::processEvents(QEventLoop::AllEvents, 50);
}

//Intensity calculation (Amplitude^2)
double Utilities::ComplexAbsSquared(std::complex<double> c)
{
    return (c.real()*c.real() + c.imag()*c.imag());
}

//Absolute value (amplitude) calculation
double Utilities::ComplexAbs(std::complex<double> c)
{
    return sqrt(ComplexAbsSquared(c));
}

//Simpson's 1/3 rule
double Utilities::SimpsonsWeight (unsigned int i, unsigned int n)
{
    if (i == 0 || i == n - 1)  return 1.0 / 3.0;
    if (i % 2 != 0)   return 4.0 / 3.0;
    return 2.0 / 3.0;
}

//Find nice step for polar plot ticks
double Utilities::NiceStep(double range, int initialCircles)
{
    double niceStep = 0.0;
    double roughStep = (initialCircles > 0) ? range / initialCircles : range;
    double exponent = floor(log10(roughStep));
    double fraction = roughStep / pow(10, exponent);

    if (fraction < 1.5)
    {
        niceStep = 1.0 * pow(10, exponent);
    }
    else
    {
        if (fraction < 3.0)
        {
            niceStep = 2.0 * pow(10, exponent);
        }
        else
        {
            if (fraction < 7.5)
            {
                niceStep = 5.0 * pow(10, exponent);
            }
            else
            {
                niceStep = 10.0 * pow(10, exponent);
            }
        }
    }
    return niceStep;
}

//Find min and max out of three variables
double Utilities::FindMinMax(const QVector<double>& yPara, const QVector<double>& yPerp,
                             const QVector<double>& yAve, bool flagMinMax)
{
    QVector<double> allY;
    allY << yPara << yPerp << yAve;

    if (allY.isEmpty())
        return std::nan("");

    if (flagMinMax)
        return *std::min_element(allY.constBegin(), allY.constEnd());
    else
        return *std::max_element(allY.constBegin(), allY.constEnd());
}
