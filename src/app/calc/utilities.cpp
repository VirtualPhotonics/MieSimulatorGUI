//**********************************************************************
//** Some utility functions that are commonly used are listed here.
//**********************************************************************

#include "calc/utilities.h"

utilities::utilities()
{
}

utilities::~utilities()
{
}

//Delay function for display update
void utilities::Delay()
{
    QTime dieTime= QTime::currentTime().addMSecs(1);
    while( QTime::currentTime() < dieTime )
    QCoreApplication::processEvents(QEventLoop::AllEvents, 50);
}

//Intensity calculation (Amplitude^2)
double utilities::ComplexAbsSquared(std::complex<double> c)
{
    return (c.real()*c.real() + c.imag()*c.imag());
}

//Absolute value (amplitude) calculation
double utilities::ComplexAbs(std::complex<double> c)
{
    return sqrt(ComplexAbsSquared(c));
}

//Simpson's 1/3 rule
double utilities::SimpsonsWeight (unsigned int i, unsigned int n)
{
    if (i == 0 || i == n - 1)  return 1.0 / 3.0;
    if (i % 2 != 0)   return 4.0 / 3.0;
    return 2.0 / 3.0;
}
