/**********************************************************************
** Mie Calculations are listed in this file. This far field solution is
** based on Bohren and Huffman's (John Wiley 1983) BHMIE code
**********************************************************************/

#include "miesimulation.h"
#include <cmath>
#include <cstdlib>
#include <ctime>

#define max(a,b)    (((a) > (b)) ? (a) : (b))

MieSimulation::MieSimulation(void)
{
}


MieSimulation::~MieSimulation(void)
{
}

//This function provides S1,S2, qSca, qExt and qBack for given xParameter, relativeRefIndex and mu(cos(angle))
//cS1 - complex S1
//cS2 - complex S2
//qSca - scattering effieicncy
//qSca - extinction effieicncy
//qBack - backscattering effieicncy
//xPara - x Parameter (2*pi*rad*refmed/wavel)
//refRel - relative refractive index
//mu - cos(angle)
void MieSimulation::FarFieldSolutionForRealRefIndex(std::complex<double> *cS1, std::complex<double> *cS2,
                                                    double *qSca, double *qExt, double *qBack,
                                                    double xPara, double relRef, double mu)
{
    utilities util;

    //use conventional symbols
    double x = xPara;
    double m = relRef;
    double mx = m*x;
    double xstop = x+4.05*(pow(x,(1.0/3.0)))+2.0;
    int nstop = ceil(xstop);
    int ymod = ceil(fabs(mx));
    int nmx = max(xstop,ymod)+15;
    int arraySize = nstop+1;
    double x2 = x*x;

    double *Dn_mx = new double [nmx];
    Dn_mx[nmx-1]=0;
    for (int N = nmx-1; N>0; N--)
        Dn_mx[N-1] = (double(N)/mx)-(1.0/(Dn_mx[N]+double(N)/mx));

    // Legendre Polynomials
    double dPCost0 = 0.0;  //pi0
    double dPCost1 = 1.0;  //pi1

    // at the sphere boundary
    double j_x0 = cos(x);  // phi(-1)
    double y_x0 = -sin(x); // kai(-1)
    double j_x1 = sin(x);  // phi(0)
    double y_x1 = cos(x);  // kai(0)

    double *j_x = new double [arraySize];
    double *y_x = new double [arraySize];
    std::complex<double> *xi_x = new std::complex<double>  [arraySize];
    j_x[0] = j_x1;
    y_x[0] = y_x1;
    xi_x[0]=std::complex<double> (j_x1,-y_x1); // xi(1)

    //Initialize temp holders
    std::complex<double> tempS1 = std::complex<double>  (0.0, 0.0);
    std::complex<double> tempS2 = std::complex<double>  (0.0, 0.0);
    std::complex<double> tempQback = 0.0;
    double tempQsca = 0.0;
    double tempQext = 0.0;

    double *piCost = new double [arraySize];
    double *tauCost = new double [arraySize];
    int n=1;
    while ((n-1-nstop) < 0)
    {
        double fac0 = double(n);		// n
        double fac1 = fac0+1.0;         // n+1
        double fac2 = 2.0*fac1 -1.0;	// 2n+1
        double fac3 = fac2 -2.0;        // 2n-1
        double fac4 = fac2/(fac0*fac1); // (2n+1)/(n*(n+1))

        //Update Legrendre polynomials (Array indices = n)
        piCost[n-1] = dPCost1;
        tauCost[n-1] = (fac0*mu*dPCost1) - (fac1*dPCost0);

        dPCost1 = (fac2*mu*dPCost1/fac0) - (fac1*dPCost0/fac0);
        dPCost0 = piCost[n-1];

        //Update riccati Bessel functions  for x (Array indices = n+1)
        j_x[n] = fac3*j_x1/x-j_x0;		// phi recurrence
        y_x[n] = fac3*y_x1/x-y_x0;		// kai recurrence
        xi_x[n] = std::complex<double> (j_x[n],-y_x[n]);

        j_x0 = j_x1;
        j_x1 = j_x[n];
        y_x0 = y_x1;
        y_x1 = y_x[n];

        // Calculate an and bn  (According to Bohren and Huffman book)
        // Remark: GouGouesbet  uses size parameter as "ka" instead of "kx"
        double dervDn1 = (Dn_mx[n]/m) + (double(n)/x);
        double dervDn2 = (m*Dn_mx[n]) + (double(n)/x);

        std::complex<double> an = (dervDn1*j_x[n]-j_x[n-1])/ (dervDn1*xi_x[n]-xi_x[n-1]);
        std::complex<double> bn = (dervDn2*j_x[n]-j_x[n-1])/ (dervDn2*xi_x[n]-xi_x[n-1]);
        tempQback += fac2*pow(-1.0,(n-1.0))*(an-bn);
        tempQsca += fac2*(util.ComplexAbs(an)*util.ComplexAbs(an) + util.ComplexAbs(bn)*util.ComplexAbs(bn));
        tempQext += fac2*(an+bn).real();

        // Calculate cS1 and cS2
        tempS1 += fac4*(an*piCost[n-1]+bn*tauCost[n-1]);
        tempS2 += fac4*(an*tauCost[n-1]+bn*piCost[n-1]);

        n = n+1;
    }    
    *qBack = util.ComplexAbsSquared(tempQback)/x2;  //back scattering efficiency
    *qSca = 2.0 * tempQsca/x2;                      //scattering efficiency
    *qExt= 2.0 * tempQext/x2;                       //extinction efficiency
    *cS1 = tempS1;
    *cS2 = tempS2;

    delete[] Dn_mx;
    delete[] tauCost;
    delete[] piCost;
    delete[] y_x;
    delete[] j_x;
    delete[] xi_x;
}

//This function provides S1,S2, qSca, qExt and qBack for given xParameter, complex relativeRefIndex and mu(cos(angle))
//cS1 - complex S1
//cS2 - complex S2
//qSca - scattering effieicncy
//qSca - extinction effieicncy
//qBack - backscattering effieicncy
//xPara - x Parameter (2*pi*rad*refmed/wavel)
//cRefRel - complex relative refractive index
//mu - cos(angle)
void MieSimulation::FarFieldSolutionForComplexRefIndex(std::complex<double> *cS1, std::complex<double> *cS2,
                                                       double *qSca, double *qExt, double *qBack, double xPara, std::complex<double> cRelRef,
                                                       double mu)
{
    utilities util;

    //use conventional symbols
    double x = xPara;
    std::complex<double> m = cRelRef;
    std::complex<double> mx = m*x;
    double xstop = x+4.05*(pow(x,(1.0/3.0)))+2.0;
    int nstop = ceil(xstop);
    int ymod = ceil(util.ComplexAbs(mx));
    int nmx = max(xstop,ymod)+15;
    int arraySize = nstop+1;
    double x2 = x*x;

    std::complex<double> *Dn_mx = new std::complex<double> [nmx];
    Dn_mx[nmx-1] = 0;
    for (int N = nmx-1; N>0; N--)
        Dn_mx[N-1] = (double(N)/mx)-(1.0/(Dn_mx[N]+double(N)/mx));

    // Legendre Polynomials
    double dPCost0 = 0.0;  //pi0
    double dPCost1 = 1.0;  //pi1

    // at the sphere boundary
    double j_x0 = cos(x);  // phi(-1)
    double y_x0 = -sin(x); // kai(-1)
    double j_x1 = sin(x);  // phi(0)
    double y_x1 = cos(x);  // kai(0)

    double *j_x = new double [arraySize];
    double *y_x = new double [arraySize];
    std::complex<double> *xi_x = new std::complex<double>  [arraySize];
    j_x[0] = j_x1;
    y_x[0] = y_x1;
    xi_x[0]=std::complex<double> (j_x1,-y_x1); // xi(1)

    std::complex<double> tempS1 = std::complex<double>  (0.0, 0.0);
    std::complex<double> tempS2 = std::complex<double>  (0.0, 0.0);
    std::complex<double> tempQback = 0.0;
    double tempQsca = 0.0;
    double tempQext = 0.0;

    double *piCost = new double [arraySize];
    double *tauCost = new double [arraySize];
    int n=1;
    while ((n-1-nstop) < 0)
    {
        double fac0 = double(n);		// n
        double fac1 = fac0+1.0;		// n+1
        double fac2 = 2.0*fac1 -1.0;	// 2n+1
        double fac3 = fac2 -2.0;		// 2n-1
        double fac4 = fac2/(fac0*fac1);// (2n+1)/(n*(n+1))

        //Update Legrendre polynomials (Array indices = n)
        piCost[n-1]=dPCost1;
        tauCost[n-1]=(fac0*mu*dPCost1) - (fac1*dPCost0);

        dPCost1 = (fac2*mu*dPCost1/fac0) - (fac1*dPCost0/fac0);
        dPCost0 = piCost[n-1];

        //Update riccati Bessel functions  for x (Array indices = n+1)
        j_x[n] = fac3*j_x1/x-j_x0;		// phi recurrence
        y_x[n] = fac3*y_x1/x-y_x0;		// kai recurrence
        xi_x[n] = std::complex<double> (j_x[n],-y_x[n]);

        j_x0 = j_x1;
        j_x1 = j_x[n];
        y_x0 = y_x1;
        y_x1 = y_x[n];

        // Calculate an and bn  (According to Bohren and Huffman book)
        // Remark: GouGouesbet  uses size parameter as "ka" instead of "kx"
        std::complex<double> dervDn1 = (Dn_mx[n]/m) + (double(n)/x);
        std::complex<double> dervDn2 = (m*Dn_mx[n]) + (double(n)/x);

        std::complex<double> an = (dervDn1*j_x[n]-j_x[n-1])/ (dervDn1*xi_x[n]-xi_x[n-1]);
        std::complex<double> bn = (dervDn2*j_x[n]-j_x[n-1])/ (dervDn2*xi_x[n]-xi_x[n-1]);
        tempQback += fac2*pow(-1.0,(n-1.0))*(an-bn);
        tempQsca += fac2*(util.ComplexAbs(an)*util.ComplexAbs(an) + util.ComplexAbs(bn)*util.ComplexAbs(bn));
        tempQext += fac2*(an+bn).real();

        // Calculate cS1 and cS2
        tempS1 += fac4*(an*piCost[n-1]+bn*tauCost[n-1]);
        tempS2 += fac4*(an*tauCost[n-1]+bn*piCost[n-1]);

        n = n+1;
    }
    *qBack = util.ComplexAbsSquared(tempQback)/x2;  //back scattering efficiency
    *qSca= 2.0 * tempQsca/x2;                       //scattering efficiency
    *qExt= 2.0 * tempQext/x2;                       //extinction efficiency
    *cS1 = tempS1;
    *cS2 = tempS2;

    delete[] Dn_mx;
    delete[] tauCost;
    delete[] piCost;
    delete[] y_x;
    delete[] j_x;
    delete[] xi_x;
}
