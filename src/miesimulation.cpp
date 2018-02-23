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

//This function provides S1,S2 and qSca for given xParameter, relativeRefIndex and mu(cos(angle))
//cS1 - complex S1
//cS2 - complex S2
//qSca - scattering effieicncy
//xPara - x Parameter (2*pi*rad*refmed/wavel)
// refRel - relative refractive index
// mu - cos(angle)
void MieSimulation::FarFieldSolutionForRealRefIndex(std::complex<double> *cS1,std::complex<double> *cS2,
                                                    double *qSca,double xPara,double relRef,double mu)
{
    double xstop;
    double mx;
    double dPCost0, dPCost1;
    double fac0, fac1, fac2, fac3, fac4;
    double j_x0, j_x1, y_x0, y_x1;    
    double dervDn1, dervDn2;
    double sumQsca;
    std::complex<double> an, bn;
    std::complex<double> tempS1, tempS2;

    utilities util;

    //use conventional symbols
    double x = xPara;
    double m = relRef;

    mx = m*x;
    xstop = x+4.0*(pow(x,(1.0/3.0)))+2.0;
    int nstop = ceil(xstop);
    int ymod = ceil(fabs(mx));
    int nmx = max(xstop,ymod)+15;
    int arraySize = nstop+1;

	double *Dn_mx = new double [nmx];
	double *j_x = new double [arraySize];
	double *y_x = new double [arraySize];
	double *piCost = new double [arraySize];
	double *tauCost = new double [arraySize];
	std::complex<double> *xi_x = new std::complex<double>  [arraySize];

    Dn_mx[nmx-1]=0;
    for (int N = nmx-1; N>0; N--)
        Dn_mx[N-1] = (N/mx)-(1.0/(Dn_mx[N]+N/mx));

    // Legendre Polynomials
    dPCost0 = 0.0;  //pi0
    dPCost1 = 1.0;  //pi1

        // at the sphere boundary
    j_x0 = cos(x);  // phi(-1)
    y_x0 = -sin(x); // kai(-1)
    j_x1 = sin(x);  // phi(0)
    y_x1 = cos(x);  // kai(0)

    j_x[0] = j_x1;
    y_x[0] = y_x1;
    xi_x[0]=std::complex<double> (j_x1,-y_x1); // xi(1)

    tempS1 = std::complex<double>  (0.0, 0.0);
    tempS2 = std::complex<double>  (0.0, 0.0);
    sumQsca = 0.0;

    int n=1;
    while ((n-1-nstop) < 0)
    {
        fac0 = double(n);		// n
        fac1 = fac0+1.0;		// n+1
        fac2 = 2.0*fac1 -1.0;		// 2n+1
        fac3 = fac2 -2.0;               // 2n-1
        fac4 = fac2/(fac0*fac1);	// (2n+1)/(n*(n+1))

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
        dervDn1 = Dn_mx[n]/m +n/x;
        dervDn2 = m*Dn_mx[n] +n/x;

        an = (dervDn1*j_x[n]-j_x[n-1])/ (dervDn1*xi_x[n]-xi_x[n-1]);
        bn = (dervDn2*j_x[n]-j_x[n-1])/ (dervDn2*xi_x[n]-xi_x[n-1]);
        sumQsca += (fac2/(x*x))*(util.ComplexAbs(an)*util.ComplexAbs(an) + util.ComplexAbs(bn)*util.ComplexAbs(bn));

                // Calculate cS1 and cS2
        tempS1 += fac4*(an*piCost[n-1]+bn*tauCost[n-1]);
        tempS2 += fac4*(an*tauCost[n-1]+bn*piCost[n-1]);

        n = n+1;
    }

    *qSca= sumQsca;
    *cS1 = tempS1;
    *cS2 = tempS2;

    delete[] Dn_mx;
    delete[] tauCost;
    delete[] piCost;
    delete[] y_x;
    delete[] j_x;
    delete[] xi_x;
}

//This function provides S1,S2 and qSca for given xParameter, complex relativeRefIndex and mu(cos(angle))
//cS1 - complex S1
//cS2 - complex S2
//qSca - scattering effieicncy
//xPara - x Parameter (2*pi*rad*refmed/wavel)
// cRefRel - complex relative refractive index
// mu - cos(angle)
void MieSimulation::FarFieldSolutionForComplexRefIndex(std::complex<double> *cS1, std::complex<double> *cS2,
                                                       double *qSca, double xPara, std::complex<double> cRelRef,
                                                       double mu)
{
    int n;
    int nstop, ymod, nmx, arraySize;
    double xstop;
    std::complex<double> mx;
    double dPCost0, dPCost1;
    double fac0, fac1, fac2, fac3, fac4;
    double j_x0, j_x1, y_x0, y_x1;
    double *j_x;
    double *y_x;
    double *piCost;
    double *tauCost;
    std::complex<double> *Dn_mx;
    std::complex<double> dervDn1, dervDn2;
    double sumQSca;
    std::complex<double> *xi_x;
    std::complex<double> an, bn;
    std::complex<double> tempS1, tempS2;

    utilities util;

    //use conventional symbols
    double x = xPara;
    std::complex<double> m = cRelRef;

    mx = m*x;
    xstop = x+4.0*(pow(x,(1.0/3.0)))+2.0;
    nstop = ceil(xstop);
    ymod = ceil(util.ComplexAbs(mx));
    nmx = max(xstop,ymod)+15;
    arraySize = nstop+1;

    Dn_mx = new std::complex<double> [nmx];
    j_x = new double [arraySize];
    y_x = new double [arraySize];
    piCost = new double [arraySize];
    tauCost = new double [arraySize];
    xi_x = new std::complex<double>  [arraySize];

    Dn_mx[nmx-1] = 0;
    for (n = nmx-1; n>0; n--)
        Dn_mx[n-1] = (double(n)/mx)-(1.0/(Dn_mx[n]+(double(n)/mx)));

    // Legendre Polynomials
    dPCost0 = 0.0;  //pi0
    dPCost1 = 1.0;  //pi1

    // at the sphere boundary
    j_x0 = cos(x);  // phi(-1)
    y_x0 = -sin(x); // kai(-1)
    j_x1 = sin(x);  // phi(0)
    y_x1 = cos(x);  // kai(0)

    j_x[0] = j_x1;
    y_x[0] = y_x1;
    xi_x[0]=std::complex<double> (j_x1,-y_x1); // xi(1)

    tempS1 = std::complex<double>  (0.0, 0.0);
    tempS2 = std::complex<double>  (0.0, 0.0);
    sumQSca = 0.0;

    n=1;
    while ((n-1-nstop) < 0)
    {
        fac0 = double(n);		// n
        fac1 = fac0+1.0;		// n+1
        fac2 = 2.0*fac1 -1.0;		// 2n+1
        fac3 = fac2 -2.0;		// 2n-1
        fac4 = fac2/(fac0*fac1);	// (2n+1)/(n*(n+1))

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
        dervDn1 = (Dn_mx[n]/m) + (double(n)/x);
        dervDn2 = (m*Dn_mx[n]) + (double(n)/x);

        an = (dervDn1*j_x[n]-j_x[n-1])/ (dervDn1*xi_x[n]-xi_x[n-1]);
        bn = (dervDn2*j_x[n]-j_x[n-1])/ (dervDn2*xi_x[n]-xi_x[n-1]);
        sumQSca += fac2*(util.ComplexAbs(an)*util.ComplexAbs(an) + util.ComplexAbs(bn)*util.ComplexAbs(bn));

        // Calculate cS1 and cS2
        tempS1 += fac4*(an*piCost[n-1]+bn*tauCost[n-1]);
        tempS2 += fac4*(an*tauCost[n-1]+bn*piCost[n-1]);

        n = n+1;
    }

    *qSca= (2.0/(x*x))*sumQSca;
    *cS1 = tempS1;
    *cS2 = tempS2;

    delete[] Dn_mx;
    delete[] tauCost;
    delete[] piCost;
    delete[] y_x;
    delete[] j_x;
    delete[] xi_x;
}
