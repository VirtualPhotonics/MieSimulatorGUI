/**********************************************************************
** Calculations for all panels are listed in this file.
**********************************************************************/

#include "calc/calculate.h"
#include "calc/miesimulation.h"
#include "calc/utilities.h"
#include "parameters.h"

Calculate::Calculate()
{
}

Calculate::~Calculate()
{
}

//Calculate parameters using mie solution
void Calculate::DoSimulation(QLabel *progress, Parameters *para)
{
    MieSimulation mieSim;
    double curMus;
    double sumMus, sumMusG;
    double sumForward, sumBackward;
    double sumCsca, sumCext, sumCback;
    double tempPhase;
    double piRadiusSquared;
    double refRelRe = 1.0;
    double refRelIm = 0.0;
    double xPara = 0.0;
    double stepTheta = (mMaxTheta - mMinTheta)/static_cast<double>(para->nTheta - 1);

    std::complex<double> *curS1 = new std::complex<double> [para->nTheta];
    std::complex<double> *curS2 = new std::complex<double> [para->nTheta];
    std::complex<double> *sumS1 = new std::complex<double> [para->nTheta];
    std::complex<double> *sumS2 = new std::complex<double> [para->nTheta];
    double *sumPhaseFuncAve = new double [para->nTheta];
    double *sumPhaseFuncPara = new double [para->nTheta];
    double *sumPhaseFuncPerp = new double [para->nTheta];

    double sumNumDen = 0.0;
    for (unsigned int r = 0; r < para->nRadius; r++)
        sumNumDen += para->numDensityArray[r];

    //Calculate all Mie parameters
    Utilities util;
    for (unsigned int w = 0; w < para->nWavel; w++)
    {
        mWavel = para->wavelArray[w] / 1000.0;   //in microns
        progress->setText("<font color=\"red\">WL: <font>"+QString::number(1000 * mWavel)+"nm</font>");
        util.Delay();        
        sumMus = 0.0;
        sumMusG = 0.0;
        sumForward = 0.0;
        sumBackward = 0.0;
        sumCsca = 0.0;
        sumCext = 0.0;
        sumCback = 0.0;
        for (unsigned int t = 0; t < para->nTheta; t++)
        {
            sumS1[t] = std::complex<double>(0.0,0.0);
            sumS2[t] = std::complex<double>(0.0,0.0);
            sumPhaseFuncAve[t] = 0.0;
            sumPhaseFuncPara[t] = 0.0;
            sumPhaseFuncPerp[t] = 0.0;
        }
        for (unsigned int r = 0; r < para->nRadius; r++)
        {
            mK = 2 * M_PI * para->medRefArray[r] / mWavel;
            xPara = mK * para->radArray[r];
            piRadiusSquared = M_PI * para->radArray[r] * para->radArray[r];

            refRelRe = para->scatRefRealArray[r] / para->medRefArray[r];
            refRelIm = para->scatRefImagArray[r] / para->medRefArray[r];

            if (refRelIm == 0.0)  //FarFieldSolutionForRealRefIndex is ~2x faster than FarFieldSolutionForComplexRefIndex
            {
                for (unsigned int t = 0; t < para->nTheta; t++)
                {
                    mMu = cos(mMinTheta + t * stepTheta);
                    mieSim.FarFieldSolutionForRealRefIndex(&mCS1, &mCS2, &mQSca, &mQExt, &mQBack, xPara, refRelRe, mMu);
                    curS1[t] = mCS1;
                    curS2[t] = mCS2;
                }
            }
            else
            {
                for (unsigned int t = 0; t < para->nTheta; t++)
                {
                    mMu = cos(mMinTheta + t * stepTheta);
                    mieSim.FarFieldSolutionForComplexRefIndex(&mCS1, &mCS2, &mQSca, &mQExt, &mQBack, xPara,
                                                           std::complex<double>(refRelRe,-refRelIm), mMu);  //multiply by -1 to use "n-ik" convention
                    curS1[t] = mCS1;
                    curS2[t] = mCS2;
                }
            }
            //Ref: Schmitt and Kumar, Applied Optics 37(13) 1998
            //Mus calculation
            curMus = piRadiusSquared * mQSca * para->numDensityArray[r];
            sumMus += curMus;          //Σμs
            //G calculation
            sumMusG += CalculateG(curS1, curS2, para) * curMus;
            //Forward Scattering
            sumForward += CalculateForwardBackward(curS1, curS2, para, 0, ((para->nTheta-1)/2)+1)* curMus;  //0-90 deg
            //Backward Scattering
            sumBackward += CalculateForwardBackward(curS1, curS2, para, ((para->nTheta-1)/2), para->nTheta)* curMus; //90-180 deg

            //Coefficients
            sumCsca += piRadiusSquared * mQSca * para->numDensityArray[r];
            sumCext += piRadiusSquared * mQExt * para->numDensityArray[r];
            sumCback += piRadiusSquared * mQBack * para->numDensityArray[r];

            //S1 and S2
            double factor = 1.0 / (M_PI * xPara * xPara * mQSca);   // 1/(pi *X^2 * qSca);
            for (unsigned int t = 0; t < para->nTheta; t++)
            {
                sumS1[t] += curS1[t] * curMus;
                sumS2[t] += curS2[t] * curMus;
                //Phase function - Parallel
                tempPhase = util.ComplexAbsSquared(curS2[t])*factor;    //|S2|^2/(pi X^2 Qsca)
                sumPhaseFuncPara[t] +=  tempPhase*curMus;
                //Phase function - Perpendicular
                tempPhase = util.ComplexAbsSquared(curS1[t])*factor;    //|S1|^2/(pi X^2 Qsca)
                sumPhaseFuncPerp[t] +=  tempPhase*curMus;
                //Phase function - Average
                tempPhase = 0.5*(util.ComplexAbsSquared(curS1[t])+util.ComplexAbsSquared(curS2[t]))*factor;  //(|S1|^2+|S2|^2)/(2 pi X^2 Qsca)
                sumPhaseFuncAve[t] +=  tempPhase*curMus;
            }
        }
        //Normalize
        para->cSca[w]  = sumCsca / sumNumDen ;
        para->cExt[w]  = sumCext / sumNumDen;
        para->cBack[w] = sumCback / sumNumDen;
        para->sizePara[w] = xPara;
        para->mus[w] = sumMus * 1e-6;   //1e-6--> 1micron2 to 1mm2
        para->g[w] = sumMusG /sumMus;
        para->forward[w] = sumForward*100.0/(sumForward+sumBackward);    //Not necessary to divide by sumMus in ratio calculation
        para->backward[w] = sumBackward*100.0/(sumForward+sumBackward);  //Not necessary to divide by sumMus in ratio calculation
        for (unsigned int t = 0; t < para->nTheta; t++)
        {
            para->S1[w][t] = sumS1[t] /sumMus;
            para->S2[w][t] = sumS2[t] /sumMus;
            para->phaseFunctionAve[w][t] = sumPhaseFuncAve[t] /sumMus;
            para->phaseFunctionPara[w][t] = sumPhaseFuncPara[t] /sumMus;
            para->phaseFunctionPerp[w][t] = sumPhaseFuncPerp[t] /sumMus;
        }
    }
    delete[] sumPhaseFuncPerp;
    delete[] sumPhaseFuncPara;
    delete[] sumPhaseFuncAve;
    delete[] sumS2;
    delete[] sumS1;
    delete[] curS2;
    delete[] curS1;
}

//Compute Musp at reference wavelength
void Calculate::ComputeMuspAtRefWavel(Parameters *para)
{
    MieSimulation mieSim;
    double curMus;
    double sumMus, sumMusG;
    double piRadiusSquared;
    double refRelRe = 1.0;
    double refRelIm = 0.0;
    double xPara = 0.0;
    double stepTheta = (mMaxTheta - mMinTheta)/static_cast<double>(para->nTheta - 1);

    std::complex<double> *curS1 = new std::complex<double> [para->nTheta];
    std::complex<double> *curS2 = new std::complex<double> [para->nTheta];
    std::complex<double> *sumS1 = new std::complex<double> [para->nTheta];
    std::complex<double> *sumS2 = new std::complex<double> [para->nTheta];

    for (unsigned int i = 0; i< 6; i++)
    {
        double lambda = 0.5 + 0.1*i;            // in um

        sumMus = 0.0;
        sumMusG = 0.0;
        for (unsigned int r = 0; r < para->nRadius; r++)
        {
            //Calculate mus and g at reference wavelength for musp fitting plot
            mK = 2 * M_PI * para->medRefArray[r] /lambda;    //k at reference wavelength

            xPara = mK * para->radArray[r];
            piRadiusSquared = M_PI * para->radArray[r] * para->radArray[r];

            refRelRe = para->scatRefRealArray[r] / para->medRefArray[r];
            refRelIm = para->scatRefImagArray[r] / para->medRefArray[r];

            if (refRelIm == 0.0)  //FarFieldSolutionForRealRefIndex is ~2x faster than FarFieldSolutionForComplexRefIndex
            {
                for (unsigned int t = 0; t < para->nTheta; t++)
                {
                    mMu = cos(mMinTheta + t * stepTheta);
                    mieSim.FarFieldSolutionForRealRefIndex(&mCS1, &mCS2, &mQSca, &mQExt, &mQBack, xPara, refRelRe, mMu);
                    curS1[t] = mCS1;
                    curS2[t] = mCS2;
                }
            }
            else
            {
                for (unsigned int t = 0; t < para->nTheta; t++)
                {
                    mMu = cos(mMinTheta + t * stepTheta);
                    mieSim.FarFieldSolutionForComplexRefIndex(&mCS1, &mCS2, &mQSca, &mQExt, &mQBack, xPara,
                                                           std::complex<double>(refRelRe,-refRelIm), mMu);  //multiply by -1 to use "n-ik" convention
                    curS1[t] = mCS1;
                    curS2[t] = mCS2;
                }
            }
            //Mus calculation
            curMus = piRadiusSquared * mQSca * para->numDensityArray[r] * 1e-6;  //1e-6--> 1micron2 to 1mm2
            sumMus += curMus;          //Σμs
            //G calculation
            sumMusG += CalculateG(curS1, curS2, para)*curMus;
        }
        para->muspAtRefWavel[i]= sumMus*(1-(sumMusG /sumMus));
    }
    delete[] sumS2;
    delete[] sumS1;
    delete[] curS2;
    delete[] curS1;
}

//Calculate bestfit and plot for Simple Algorithm
void Calculate::CalculatePowerLawAutoFitSimple(Parameters *para)
{
    double bMie, yFit;
    double error, sumError;
    double minError = 1e100;
    double curB = 0.0;
    double x, y;

   for (int j=0; j<=400; j++)
   {
       bMie = j*0.01;   //Range: [0, 4]
       sumError = 0.0;
       for (unsigned int k=0; k<para->nWavel; k++)
       {
           //Steve L Jacques,"Optical properties of biological tissues: a review" Phys. Med & Bio. 58(2013) R37-R61.
           //wavelength λ is normalized by a reference wavelength
           x = para->wavelArray[k];
           yFit = ( para->muspAtRefWavel[para->refWavelIdx]*pow(x/para->refWavel, -bMie));

           y = para->mus[k] * (1.0 - para->g[k]);
           error = yFit - y;
           sumError += error*error;
       }
       if (sumError < minError)
       {
           minError = sumError;
           curB = bMie;
       }
   }
    para->bMie = curB;
}

//Calculate bestfit and plot for Complex algorithm
void Calculate::CalculatePowerLawAutoFitComplex(Parameters *para)
{
    double bMie, fRay, yFit;
    double error, sumError;
    double minError = 1e100;
    double curB = 0.0, curF = 0.0;
    double x, y;

    for (int i=0; i<=100; i++)
    {
       fRay = i*0.01;       //Range: [0, 1]
       for (int j=0; j<=400; j++)
       {
           bMie = j*0.01;   //Range: [0, 4]
           sumError = 0.0;
           for (unsigned int k=0; k<para->nWavel; k++)
           {
               //Steve L Jacques,"Optical properties of biological tissues: a review" Phys. Med & Bio. 58(2013) R37-R61.
               //wavelength λ is normalized by a reference wavelength, 1000 nm
               x = para->wavelArray[k];
               yFit = ( para->muspAtRefWavel[para->refWavelIdx] *(fRay*pow(x/para->refWavel, -4.0) + (1.0-fRay)*pow(x/para->refWavel, -bMie)));

               y = para->mus[k] * (1.0 - para->g[k]);
               error = yFit - y;
               sumError += error*error;
           }
           if (sumError < minError)
           {
               minError = sumError;
               curB = bMie;
               curF = fRay;
           }
       }
    }
    //When b=4, fRay can take any value between 0 and 1. This condition sets fRay to unity.
    if (curB == 4.0)
    {
       curF = 1.0;
       curB = 0.0;
    }
    para->fRay = curF;
    para->bMie = curB;
}

//Calculate forward and backward scattering percentage
double Calculate::CalculateForwardBackward(std::complex<double> *S1,
                                           std::complex<double> *S2,
                                           Parameters *para,
                                           unsigned int start,
                                           unsigned int end)
{
    double S;
    double theta, twoPiSinTheta;
    double sum = 0.0;
    double stepTheta = (mMaxTheta - mMinTheta)/static_cast<double>(para->nTheta - 1);

    Utilities util;
    for (unsigned int i = start; i <end; i++)
    {
        S = (util.ComplexAbsSquared(S1[i]) + util.ComplexAbsSquared(S2[i]))/2.0;
        theta = mMinTheta + i * stepTheta;
        twoPiSinTheta = 2.0 * M_PI * sin(theta);
        sum += S * twoPiSinTheta * util.SimpsonsWeight (i-start, end-start);
    }
    return sum;
}

//Calculate average cosine of phase function
double Calculate::CalculateG(std::complex<double> *S1, std::complex<double> *S2, Parameters *para)
{
    double sVal;
    double theta, twoPiSinTheta;
    double num = 0.0;
    double den = 0.0;
    double stepTheta = (mMaxTheta - mMinTheta)/static_cast<double>(para->nTheta - 1);

    Utilities util;
    for (unsigned int i = 0; i < para->nTheta; i++)
    {
        sVal = (util.ComplexAbsSquared(S1[i]) + util.ComplexAbsSquared(S2[i]))/2.0;
        theta = mMinTheta + i * stepTheta;
        twoPiSinTheta = 2.0 * M_PI * sin(theta);
        num += sVal * cos(theta) * twoPiSinTheta * util.SimpsonsWeight (i, para->nTheta);
        den += sVal * twoPiSinTheta * util.SimpsonsWeight (i, para->nTheta);
    }
    return num/den;
}

//Set sphere parameters
void Calculate::SetSphereRadiusAndRefIndex(Parameters *para, unsigned int index, bool flagVolOrConc)
{
    const double volumeConst = 4.0 * M_PI / 3.0;
    const double sqrtTwoPi = sqrt(2.0 * M_PI);
    const double sig = para->stdDev;
    const double twoSigSq = 2.0 * sig * sig;

    if (para->nRadius == 1)  //Mono Disperse
    {
        para->radArray[0] = para->meanRadius;
        if (flagVolOrConc)  //If volume fraction is selected, update number density
        {
            double singleSphereVolume = volumeConst * pow(para->meanRadius, 3);
            para->sphNumDensity = std::round(1e9 * para->volFraction / singleSphereVolume);
        }        

        para->numDensityArray[0] = para->sphNumDensity;
        para->scatRefRealArray[0] = para->scatRefReal;
        para->scatRefImagArray[0] = para->scatRefImag;
        para->medRefArray[0] = para->medRef;
    }
    else                    //Poly Disperse
    {        
        double totalSphereVolume = 0.0;
        double totalFuncSum = 0.0;
        std::vector<double> funcArray(para->nRadius, 0.0);
        double stepR = (para->maxRadius - para->minRadius) / (para->nRadius - 1);

        for (unsigned int i = 0; i < para->nRadius; i++)
        {
            double r = para->minRadius + i * stepR;
            para->radArray[i] = r;

            if (index == 0) // Log-Normal
            {
                double lnR = std::log(r);
                double lnMean = std::log(para->meanRadius);
                double diff = lnR - lnMean;
                funcArray[i] = (std::exp(-(diff * diff) / twoSigSq)) / (r * sig * sqrtTwoPi);
            }
            else if (index == 1) // Gaussian
            {
                double diff = r - para->meanRadius;
                funcArray[i] = (std::exp(-(diff * diff) / twoSigSq)) / (sig * sqrtTwoPi);
            }
            totalSphereVolume += funcArray[i] * volumeConst * pow(r, 3);
            totalFuncSum += funcArray[i];
        }

        // Compute factor to compute number density
        double factor = 1e-12;
        if (flagVolOrConc && totalSphereVolume > 0)
        {
            factor = 1e9 * para->volFraction /totalSphereVolume;
        }
        else if (totalFuncSum > 0)
        {
            factor = para->sphNumDensity/totalFuncSum;
        }

        // Apply factor to compute data
        for (unsigned int i=0; i<para->nRadius; i++)
        {
            para->numDensityArray[i] = std::round(funcArray[i]*factor);
            para->scatRefRealArray[i] = para->scatRefReal;
            para->scatRefImagArray[i] = para->scatRefImag;
            para->medRefArray[i] = para->medRef;
        }        
    }
}

//Selction of discrete sphere sizes for polydisperse distribution
//This process is used to obtain the best distribution for assigned mean diameter
void Calculate::DiameterRangeSetting(Parameters *para, unsigned int index)
{
    // index: 0 = Log-normal, 1 = Gaussian, 3 = Monodisperse
    if (index == 3)
    {
        para->minRadius = para->meanRadius;
        para->maxRadius = para->meanRadius;
        return;
    }

    double mu, sigma;
    sigma = para->stdDev;

    // Define how much of the distribution tail to include.
    // A Z-score of 3.29 covers 99.9% of a Normal distribution.
    // Scale it based on nRadius if fewer bins require tighter bounds.
    double zScore = 3.0 + (para->nRadius / 100.0);

    switch (index)
    {
        case 0: // Log-normal Distribution
        {
            // If para->meanRadius is the Arithmetic Mean (E[X]),
            // convert it to the underlying Normal mu.
            // E[X] = exp(mu + sigma^2 / 2)
            mu = std::log(para->meanRadius) - (sigma * sigma / 2.0);

            para->minRadius = std::exp(mu - zScore * sigma);
            para->maxRadius = std::exp(mu + zScore * sigma);

            // Safety floor
            if (para->minRadius < 1e-10) para->minRadius = 1e-10;
            break;
        }

        case 1: // Gaussian (Normal) Distribution
        {
            mu = para->meanRadius;

            para->minRadius = mu - zScore * sigma;
            para->maxRadius = mu + zScore * sigma;

            // Physical constraint: Radius cannot be negative
            if (para->minRadius < 1e-10) para->minRadius = 1e-10;
            break;
        }
    }
}

// Determines the scattering regime based on Tien et. al, A.R. Heat Trandfer 1(1987) & Galy et al. JQSRT 246(2020)
bool Calculate::CheckIndependentScattering(Parameters *para, double &clearanceToWavelength, double &sizeParameter,
                                           double &volFraction, double &criticalWavength, QString &strRegime)
{
    double const volumeConstant = (4.0/3.0) * M_PI ;
    double effectiveRadius = 0.0;
    double interParticleDistance;

    // Calculate wavelengths in microns
    double criticalWavelength = para->endWavel / (para->medRef * 1000.0);

    if (para->nRadius == 1)       //monodisperse
    {
        effectiveRadius = para->meanRadius;
        double singleSphVolume = volumeConstant * pow(effectiveRadius, 3);
        volFraction = singleSphVolume * para->numDensityArray[0] / 1e9;
        interParticleDistance = 1e3 / pow(para->numDensityArray[0], 1.0 / 3.0);
    }
    else                         //polydisperse
    {
        double totalVolume = 0.0;
        double totalNumDensity = 0.0;
        for (unsigned int i = 0; i < para->nRadius; i++)
        {
            double singleSphVolume = volumeConstant * pow(para->radArray[i], 3);
            totalVolume += singleSphVolume * para->numDensityArray[i];
            totalNumDensity += para->numDensityArray[i];
        }
        volFraction = totalVolume / 1e9;
        interParticleDistance = 1e3 / pow(totalNumDensity, 1.0 / 3.0);

        double averageVolume = totalVolume / totalNumDensity;
        effectiveRadius = pow(3.0 * averageVolume / (4.0 * M_PI), 1.0/3.0);
    }

    clearanceToWavelength = (interParticleDistance - 2 * effectiveRadius) / criticalWavelength;
    sizeParameter = 2.0 * M_PI * effectiveRadius / criticalWavelength;

    // Determine the threshold for clearance based on the size parameter (Galy 2020)
    double requiredClearance = (sizeParameter <= 2.0) ? 2.0 : 5.0;

    bool isDependent = false;

    if (volFraction > 0.1) // High concentration regime
    {
        strRegime = "High Concentration Regime";
        isDependent = true;
    }
    else if (volFraction > 0.006) // Transitional regime
    {
        strRegime = "Transitional Regime";
        isDependent = (clearanceToWavelength <= requiredClearance);
    }
    else // Low concentration regime
    {
        strRegime = "Low Concentration Regime";
        isDependent = (sizeParameter > 0.388) ? (clearanceToWavelength <= requiredClearance) : false;
    }
    return isDependent;
}
