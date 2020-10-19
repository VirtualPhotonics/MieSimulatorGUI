#include "calculate.h"


calculate::calculate()
{
}

calculate::~calculate()
{
}

//Calculate parameters using mie solution
void calculate::DoSimulation(Ui_MainWindow *ui, parameters *para)
{
    MieSimulation sim;
    double curMus;
    double sumMus, sumMusG;
    double sumForward, sumBackward;
    double sumCsca, sumCext, sumCback;
    double tempPhase;
    double piRadiusSquared;
    double refRelRe = 1.0;
    double refRelIm = 0.0;
    double xPara = 0.0;

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
    utilities util;
    for (unsigned int w = 0; w < para->nWavel; w++)
    {
        wavel = para->wavelArray[w]/1000;   //in microns
        ui->label_Progress->setText("<font color=\"red\">WL: <font>"+QString::number(1000*wavel)+"nm</font>");
        util.Delay();
        k = 2 * M_PI * para->medRef /wavel;
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
            xPara = k * para->radArray[r];            
            piRadiusSquared = M_PI * para->radArray[r] * para->radArray[r];

            refRelRe = para->scatRefRealArray[r] / para->medRef;
            refRelIm = para->scatRefImagArray[r] / para->medRef;

            if (refRelIm == 0.0)  //FarFieldSolutionForRealRefIndex is ~2x faster than FarFieldSolutionForComplexRefIndex
            {
                for (unsigned int t = 0; t < para->nTheta; t++)
                {
                    mu = cos(para->minTheta + t * para->stepTheta);
                    sim.FarFieldSolutionForRealRefIndex(&cS1, &cS2, &qSca, &qExt, &qBack, xPara, refRelRe, mu);
                    curS1[t] = cS1;
                    curS2[t] = cS2;
                }
            }
            else
            {
                for (unsigned int t = 0; t < para->nTheta; t++)
                {
                    mu = cos(para->minTheta + t * para->stepTheta);
                    sim.FarFieldSolutionForComplexRefIndex(&cS1, &cS2, &qSca, &qExt, &qBack, xPara,
                                                           std::complex<double>(refRelRe,-refRelIm), mu);  //multiply by -1 to use "n-ik" convention
                    curS1[t] = cS1;
                    curS2[t] = cS2;
                }
            }
            //Ref: Schmitt and Kumar, Applied Optics 37(13) 1998
            //Mus calculation
            curMus = piRadiusSquared * qSca * para->numDensityArray[r];
            sumMus += curMus;          //Σμs
            //G calculation
            sumMusG += CalculateG(curS1, curS2, para) * curMus;
            //Forward Scattering
            sumForward += CalculateForwardBackward(curS1, curS2, para, 0, (para->nTheta-1)/2)* curMus;
            //Backward Scattering
            sumBackward += CalculateForwardBackward(curS1, curS2, para, ((para->nTheta-1)/2), para->nTheta)* curMus;
            //Coefficients
            sumCsca += piRadiusSquared * qSca * para->numDensityArray[r];
            sumCext += piRadiusSquared * qExt * para->numDensityArray[r];
            sumCback += piRadiusSquared * qBack * para->numDensityArray[r];

            //S1 and S2
            double factor = 1.0 / (M_PI * xPara * xPara * qSca);   // 1/(pi *X^2 * qSca);
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
        para->SizePara[w] = xPara;
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

void calculate::ComputeMuspAtRefWavel(parameters *para)
{
    MieSimulation sim;
    double curMus;
    double sumMus, sumMusG;
    double piRadiusSquared;
    double refRelRe = 1.0;
    double refRelIm = 0.0;
    double xPara = 0.0;

    std::complex<double> *curS1 = new std::complex<double> [para->nTheta];
    std::complex<double> *curS2 = new std::complex<double> [para->nTheta];
    std::complex<double> *sumS1 = new std::complex<double> [para->nTheta];
    std::complex<double> *sumS2 = new std::complex<double> [para->nTheta];

    for (unsigned int i = 0; i< 6; i++)
    {
        double lambda = 0.5 + 0.1*i;            // in um
        //Calculate mus and g at reference wavelength for musp fitting plot
        k = 2 * M_PI * para->medRef /lambda;    //k at reference wavelength

        sumMus = 0.0;
        sumMusG = 0.0;
        for (unsigned int r = 0; r < para->nRadius; r++)
        {
            xPara = k * para->radArray[r];
            piRadiusSquared = M_PI * para->radArray[r] * para->radArray[r];

            refRelRe = para->scatRefRealArray[r] / para->medRef;
            refRelIm = para->scatRefImagArray[r] / para->medRef;

            if (refRelIm == 0.0)  //FarFieldSolutionForRealRefIndex is ~2x faster than FarFieldSolutionForComplexRefIndex
            {
                for (unsigned int t = 0; t < para->nTheta; t++)
                {
                    mu = cos(para->minTheta + t * para->stepTheta);
                    sim.FarFieldSolutionForRealRefIndex(&cS1, &cS2, &qSca, &qExt, &qBack, xPara, refRelRe, mu);
                    curS1[t] = cS1;
                    curS2[t] = cS2;
                }
            }
            else
            {
                for (unsigned int t = 0; t < para->nTheta; t++)
                {
                    mu = cos(para->minTheta + t * para->stepTheta);
                    sim.FarFieldSolutionForComplexRefIndex(&cS1, &cS2, &qSca, &qExt, &qBack, xPara,
                                                           std::complex<double>(refRelRe,-refRelIm), mu);  //multiply by -1 to use "n-ik" convention
                    curS1[t] = cS1;
                    curS2[t] = cS2;
                }
            }
            //Mus calculation
            curMus = piRadiusSquared * qSca * para->numDensityArray[r] * 1e-6;  //1e-6--> 1micron2 to 1mm2
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
void calculate::CalculatePowerLawAutoFitSimple(parameters *para)
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
void calculate::CalculatePowerLawAutoFitComplex(parameters *para)
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

//Calculate forward
double calculate::CalculateForwardBackward(std::complex<double> *S1,
                                           std::complex<double> *S2,
                                           parameters *para,
                                           unsigned int start,
                                           unsigned int end)
{
    double S;
    double theta, twoPiSinTheta;
    double sum = 0.0;

    utilities util;
    for (unsigned int i = start; i <end; i++)
    {
        S = (util.ComplexAbsSquared(S1[i]) + util.ComplexAbsSquared(S2[i]))/2.0;
        theta = para->minTheta + i * para->stepTheta;
        twoPiSinTheta = 2.0 * M_PI * sin(theta);
        sum += S * twoPiSinTheta * util.SimpsonsWeight (i-start, end-start);
    }
    return sum;
}

//Calculate average
double calculate::CalculateG(std::complex<double> *S1, std::complex<double> *S2, parameters *para)
{
    double S;
    double theta, twoPiSinTheta;
    double num = 0.0;
    double den = 0.0;

    utilities util;
    for (unsigned int i = 0; i < para->nTheta; i++)
    {
        S = (util.ComplexAbsSquared(S1[i]) + util.ComplexAbsSquared(S2[i]))/2.0;
        theta = para->minTheta + i * para->stepTheta;
        twoPiSinTheta = 2.0 * M_PI * sin(theta);
        num += S * cos(theta) * twoPiSinTheta * util.SimpsonsWeight (i, para->nTheta);
        den += S * twoPiSinTheta * util.SimpsonsWeight (i, para->nTheta);
    }
    return num/den;
}

//Set sphere parameters
void calculate::SetSphereRadiusAndRefIndex(parameters *para, unsigned int index, bool flagVolOrConc)
{
    double temp, factor;
    double totalSphereVolume = 0.0;
    double totalFuncSum = 0.0;
    double stepR;
    double *funcArray;

    if (para->nRadius == 1)  //Mono Disperse
    {
        para->radArray[0] = para->meanRadius;
        para->numDensityArray[0] = para->sphNumDensity;
        para->scatRefRealArray[0] = para->scatRefReal;
        para->scatRefImagArray[0] = para->scatRefImag;
    }
    else                    //Poly Disperse
    {
        funcArray = new double[para->nRadius];
        for (unsigned int i = 0; i < para->nRadius; i++)
			funcArray[i] = 1.0;
        stepR = (para->maxRadius - para->minRadius)/(para->nRadius -1);

        switch (index)
        {
        case 0:     //Apply log normal distribution
            for (unsigned int i=0; i<para->nRadius; i++)
            {
                para->radArray[i] = para->minRadius + i*stepR;
                temp = log(para->radArray[i])-log(para->meanRadius);
                funcArray[i] = (exp(-(temp*temp)/(2.0 * para->stdDev * para->stdDev)))/
                                (para->radArray[i] * para->stdDev * sqrt (2.0 * M_PI)) ;
                totalSphereVolume += funcArray[i] * 4.0 * M_PI * para->radArray[i] *
                                 para->radArray[i] * para->radArray[i] /3.0;
                totalFuncSum += funcArray[i];
            }
            break;
        case 1:     //Apply Gaussian distribution
            for (unsigned int i=0; i<para->nRadius; i++)
            {
                para->radArray[i] = para->minRadius + i*stepR;
                temp = para->radArray[i]-para->meanRadius;
                funcArray[i]  = (exp(-(temp*temp/(2.0 * para->stdDev * para->stdDev))))/
                                   (para->stdDev * sqrt (2.0 * M_PI)) ;
                totalSphereVolume += funcArray[i] * 4.0 * M_PI * para->radArray[i] *
                                   para->radArray[i] * para->radArray[i] /3.0;
                totalFuncSum += funcArray[i];
            }
            break;
        }

        if (flagVolOrConc)
            factor = 1e9 * para->volFraction /totalSphereVolume;   //1mm3 x volume Fraction / Total volume of spheres
        else
            factor = para->sphNumDensity/totalFuncSum;

        for (unsigned int i=0; i<para->nRadius; i++)
        {
            para->numDensityArray[i] = funcArray[i]*factor;
            para->scatRefRealArray[i] = para->scatRefReal;
            para->scatRefImagArray[i] = para->scatRefImag;
        }
        delete []funcArray;
    }
}

//Selction of discrete sphere sizes for poly disperse distribution
//This process is used to obtain the best distribution for assigned mean diameter
void calculate::DiameterRangeSetting(parameters *para, unsigned int index)
{
    double curR, dR;
    double cutoffPercent = 0.0;
    double mu, sigma = 0.0;
    double modeR;
    double curY, maxY, minY;
    int i;

    if (index == 3) //Mono sphere Distribution
    {
        para->minRadius = para->meanRadius;
        para->maxRadius = para->meanRadius;
    }
    else
    {
        //Calculate a percentage to cutoff points for Log normal and Gaussian distribution
        if (para->nRadius <52)
            cutoffPercent = -0.0005*para->nRadius + 0.031;   //from 2:52 change from  3% to 0.5%
        else
            cutoffPercent = -0.00009*para->nRadius + 0.00968;   //from 52:102 change from  0.5% to 0.05%
        sigma = para->stdDev;
    }

    switch (index)
    {
    case 0:     //Apply log normal distribution.  Source: http://en.wikipedia.org/wiki/Log-normal_distribution

        //find minR and maxR for a stand curve of mu=0 and multiply it by mean radius
        mu = 0;
        modeR = exp(mu - sigma*sigma);
        maxY = (exp(-((log(modeR)-mu)*(log(modeR)-mu)/(2.0 * sigma * sigma))))/(modeR * sigma * sqrt (2.0 * M_PI));
        minY = maxY*cutoffPercent;

        //backward search
        i = 1;
        dR = sigma*modeR/1e3;
        do
        {
            curR = modeR - i*dR;
            curY = (exp(-((log(curR)-mu)*(log(curR)-mu)/(2.0 * sigma * sigma))))/(curR * sigma * sqrt (2.0 * M_PI));
            i++;
        }
        while (curY>minY);
        para->minRadius = curR*para->meanRadius;
        if (para->minRadius <= 0.0)
            para->minRadius = 1e-10;

        //forward search
        i = 1;
        dR = sigma*modeR/1e3;
        do
        {
            curR = modeR + i*dR;
            curY = (exp(-((log(curR)-mu)*(log(curR)-mu)/(2.0 * sigma * sigma))))/(curR * sigma * sqrt (2.0 * M_PI));
            i++;
        }
        while (curY>minY);
        para->maxRadius = curR*para->meanRadius;
        break;

    case 1:     //Apply Gaussian distribution

        mu = 0;
        maxY = 1.0/(sigma * sqrt (2.0 * M_PI));
        minY = maxY*cutoffPercent;
        dR = sigma/1e3;

        //forward search to find right end
        i = 1;
        do
        {
            curR = mu + i*dR;
            curY = (exp(-((curR-mu)*(curR-mu)/(2.0 * sigma * sigma))))/(sigma * sqrt (2.0 * M_PI));
            i++;
        }
        while (curY>minY);
        para->maxRadius = para->meanRadius + curR;
        //Left end setting
        if (curR < para->meanRadius)
            para->minRadius = para->meanRadius - curR;
        else
            para->minRadius = 1e-8;
        break;    
    }
}
