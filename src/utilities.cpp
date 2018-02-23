/**********************************************************************
** Some utility functions that are coommonly used, are listed here
**********************************************************************/

#include "utilities.h"

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
    QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
}

//Intensity calculation (Amplitude^2)
double utilities::ComplexAbsSquared(std::complex<double> a)
{
    return (a.real()*a.real() + a.imag()*a.imag());
}

//Absolute value (amplitude) calculation
double utilities::ComplexAbs(std::complex<double> a)
{
    return sqrt(ComplexAbsSquared(a));
}

//Simpson's 1/3 rule
double utilities::SimpsonsWeight (int i, int n)
{
    double weight = 1.0;
    if ((i!=0) && (i!=n-1))
    {
        if ((i+1)%2==0)
            weight = 4.0;
        else
            weight = 2.0;
    }
    return weight /3.0;
}
