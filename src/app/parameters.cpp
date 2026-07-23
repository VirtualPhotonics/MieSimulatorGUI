//**********************************************************************
//** "parameters" variable is used to store all variables in a single
//** variable. This allows easy portability.
//**********************************************************************

#include "parameters.h"
#include <QComboBox>
#include <QRadioButton>
#include <QMessageBox>
#include <vector>

//Check the validity of common Parameters
bool Parameters::CheckCommonParameters(QRadioButton *radioButton_MonoDisperse,
                                       QRadioButton *radioButton_NumDen,
                                       QRadioButton *radioButton_VolFrac,
                                       QComboBox *comboBox_Distribution) const
{
    bool monoDisperseSelection = false;
    bool numDenSelection = false;
    bool volFracSelection = false;

    if (radioButton_MonoDisperse->isChecked())
    {
        monoDisperseSelection = true;
    }
    if (radioButton_NumDen->isChecked())
    {
        numDenSelection = true;
    }
    if (radioButton_VolFrac->isChecked())
    {
        volFracSelection = true;
    }

    //comboBox item order (Log Normal=0, Gaussian=1, Custom=2) matches distType enum values
    int comboBoxIndex = comboBox_Distribution->currentIndex();

    ParameterValidationResult check = CheckValidityCommonParameters(monoDisperseSelection,
                                                                    numDenSelection,
                                                                    volFracSelection,
                                                                    comboBoxIndex);
    if (!check.isValid)
    {
        QMessageBox msgBoxError;
        msgBoxError.setWindowTitle("Error");
        msgBoxError.setIcon(QMessageBox::Critical);
        msgBoxError.setText(check.errorMessage);
        msgBoxError.exec();
        return false;
    }
    else
    {
        return true;
    }
}

ParameterValidationResult Parameters::CheckValidityCommonParameters(
    bool monoDisperseSelection,
    bool numDenSelection,
    bool volFracSelection,
    int comboBoxIndex) const
{
    ParameterValidationResult result;
    result.isValid = true;

    if ((scatRefReal <= 0.0) || (medRef <= 0.0))
    {
        result.isValid = false;
        result.errorMessage = "Refractive index (Real) cannot be zero.";
        return result;
    }
    if ((scatRefReal/medRef == 1.0))
    {
        result.isValid = false;
        result.errorMessage = "Relative refractive index cannot be 1.0.";
        return result;
    }
    double m = scatRefReal / medRef;
    if ((m < 0.05) || (m > 5.0))
    {
        result.isValid = false;
        result.errorMessage = "Unrealistic relative refractive index! Check sphere and medium refractive index values.";
        return result;
    }
    if (scatRefImag >= 5.0)
    {
        result.isValid = false;
        result.errorMessage = "Imaginary refractive index must be negative or less than 5.0.";
        return result;
    }
    //Avoid zeros
    if ((startWavel <= 0.0) || (endWavel <= 0.0))
    {
        result.isValid = false;
        result.errorMessage = "The starting or ending wavelength cannot be zero.";
        return result;
    }
    //Avoid zeros
    if (stepWavel <= 0.0)
    {
        result.isValid = false;
        result.errorMessage = "Wavelength step cannot be zero.";
        return result;
    }
    if (startWavel < 50.0)
    {
        result.isValid = false;
        result.errorMessage = "Current minimum wavlength is 50nm.";
        return result;
    }
    if (endWavel > 3000.0)
    {
        result.isValid = false;
        result.errorMessage = "Current maximum wavlength is 3000nm.";
        return result;
    }
    if (endWavel - startWavel < 0.0)
    {
        result.isValid = false;
        result.errorMessage = "The starting wavelength is greater than the ending wavelength.";
        return result;
    }
    if (numDenSelection)
    {
        if (sphNumDensity < 1.0)
        {
            result.isValid = false;
            result.errorMessage = "Sphere concentration cannot be less than 1.";
            return result;
        }
    }
    if (volFracSelection)
    {
        if (volFraction <= 0.0)
        {
            result.isValid = false;
            result.errorMessage = "Volume Fraction cannot be zero.";
            return result;
        }

        if (volFraction >= 0.74048)  //Maximum packing factor = PI/(3*sqrt(2))
        {
            result.isValid = false;
            result.errorMessage = "Volume Fraction must not exceed the maximum packing factor (~0.74).";
            return result;
        }
    }
    if ((meanRadius < 0.00005) ||(meanRadius >150))
    {
        result.isValid = false;
        result.errorMessage = "Diameter is out of range! Enter a value between 0.0001μm and 300μm.";
        return result;
    }
    if (monoDisperseSelection)
    {
        if (numDenSelection)
        {
            double singleSphVolume = (4.0/3.0) * M_PI * pow(meanRadius, 3);
            if (sphNumDensity*singleSphVolume/1e9 >= M_PI/(3*sqrt(2)))  //Maximum packing factor = PI/(3*sqrt(2))
            {
                result.isValid = false;
                result.errorMessage = "'Concentration x Sphere Volume' exceeds the maximum packing factor! Reduce Concentration (Conc).";
                return result;
            }
        }
        if (volFracSelection)
        {
            double singleSphVolume = (4.0/3.0) * M_PI * pow(meanRadius, 3);
            double sphNumDensity = std::round(1e9 * volFraction / singleSphVolume);
            if (sphNumDensity < 1)
            {
                result.isValid = false;
                result.errorMessage = "Unrealistc Volume Fraction! Number density is less than 1";
                return result;
            }
        }
    }
    else //Poly Disperse
    {
        if (numDenSelection && comboBoxIndex != Custom)
        {
            double totalNumDensity = static_cast<double>(nRadius);
            if (sphNumDensity < totalNumDensity)
            {
                result.isValid = false;
                result.errorMessage = "Concentration (Ns) is too low! Increase concentration. ";
                return result;
            }
        }
        if (volFracSelection && comboBoxIndex != Custom)
        {
            double sigma = stdDev;
            double zScore = 3.0 + (nRadius / 100.0);
            double localMinRadius, localMaxRadius;

            if (comboBoxIndex == LogNormal)
            {
                double mu = std::log(meanRadius) - (sigma * sigma / 2.0);
                localMinRadius = std::exp(mu - zScore * sigma);
                localMaxRadius = std::exp(mu + zScore * sigma);
            }
            else //Gaussian
            {
                localMinRadius = meanRadius - zScore * sigma;
                localMaxRadius = meanRadius + zScore * sigma;
            }
            if (localMinRadius < 1e-10) localMinRadius = 1e-10;

            double stepR = (localMaxRadius - localMinRadius) / (nRadius - 1);
            double totalNumDensity = static_cast<double>(nRadius);
            double totalVolume = 0.0;
            for (unsigned int i = 0; i < nRadius; i++)
            {
                double r = localMinRadius + i * stepR;
                double singleSphVolume = (4.0/3.0) * M_PI * pow(r, 3);
                totalVolume += singleSphVolume * 1.0;  //numDensity assumed = 1 per bin
            }
            double computedMinVolFraction = totalVolume / 1e9;
            Q_UNUSED(totalNumDensity);

            if (volFraction < computedMinVolFraction)
            {
                result.isValid = false;
                result.errorMessage = "Unrealistic Volume Fraction! Volume Fraction is too low";
                return result;
            }
        }
    }
    return result;
}

//Check the validity of Distribution parameters
bool Parameters::CheckDistributionParameters(QComboBox *comboBox_Distribution) const
{
    int comboBoxIndex =0;

    if (comboBox_Distribution->currentIndex() == LogNormal)
        comboBoxIndex = 0;
    if (comboBox_Distribution->currentIndex() == Gaussian)
        comboBoxIndex = 1;
    if (comboBox_Distribution->currentIndex() == Custom)
        comboBoxIndex = 2;

    ParameterValidationResult check = CheckValidityDistributionParameters(comboBoxIndex);

    if (!check.isValid)
    {
        QMessageBox msgBoxError;
        msgBoxError.setWindowTitle("Error");
        msgBoxError.setIcon(QMessageBox::Critical);
        msgBoxError.setText(check.errorMessage);
        msgBoxError.exec();
        return false;
    }
    else
    {
        return true;
    }
}

ParameterValidationResult Parameters::CheckValidityDistributionParameters(int comboBoxIndex) const
{
    ParameterValidationResult result;
    result.isValid = true;

    if (stdDev == 0.0)
    {
        result.isValid = false;
        result.errorMessage = "Standard Deviation is zero! Use 'Mono Disperse'.";
        return result;
    }
    if(comboBoxIndex == 0)
    {
        if (stdDev > 3.0)
        {
            result.isValid = false;
            result.errorMessage = "Large standard deviation provides an abnormal Log Normal distribution! Current limit for Log Normal is 3.0μm.";
            return result;
        }
        if (stdDev < 1e-5)
        {
            result.isValid = false;
            result.errorMessage = "The standard deviation is too small! Use 'Mono Disperse'.";
            return result;
        }
    }
    if(comboBoxIndex == 1)
    {
        if (stdDev > 50.0)
        {
            result.isValid = false;
            result.errorMessage = "The standard deviation is too large! Current limit for Gaussian is 50.0μm.";
            return result;
        }
        if (stdDev < 1e-8)
        {
            result.isValid = false;
            result.errorMessage = "The standard deviation is too small! Use 'Mono Disperse'.";
            return result;
        }
    }
    if (nRadius == 1)
    {
        result.isValid = false;
        result.errorMessage = "Discrete sphere size is 1! Use 'Mono Disperse'.";
        return result;
    }
    if ((nRadius < 2.0) ||(nRadius >201.0))
    {
        result.isValid = false;
        result.errorMessage = "Number of sphere sizes is out of range! Enter a value between 2 and 201.";
        return result;
    }
    if ((meanRadius < 0.0005) ||(meanRadius >25))
    {
        result.isValid = false;
        result.errorMessage = "Diameter is out of range! Enter a value between 0.001μm and 50μm.";
        return result;
    }
    if (stdDev/meanRadius < 1.999e-5)
    {
        result.isValid = false;
        result.errorMessage = "Standard deviation to mean diameter ratio is smaller than 1e-5! Use 'Mono Disperse'.";
        return result;
    }
    return result;
}


bool Parameters::CheckNumDensityDistribution(QComboBox *comboBox_Distribution,
                                             QRadioButton *radioButton_NumDen,
                                             QRadioButton *radioButton_VolFrac) const
{
    int comboBoxIndex = 0;
    if (comboBox_Distribution->currentIndex() == LogNormal)
        comboBoxIndex = LogNormal;
    if (comboBox_Distribution->currentIndex() == Gaussian)
        comboBoxIndex = Gaussian;
    if (comboBox_Distribution->currentIndex() == Custom)
        comboBoxIndex = Custom;

    bool flagVolOrConc = radioButton_VolFrac->isChecked();
    Q_UNUSED(radioButton_NumDen);

    ParameterValidationResult check = CheckValidityNumDensityDistribution(comboBoxIndex, flagVolOrConc);

    if (!check.isValid)
    {
        QMessageBox msgBoxError;
        msgBoxError.setWindowTitle("Error");
        msgBoxError.setIcon(QMessageBox::Critical);
        msgBoxError.setText(check.errorMessage);
        msgBoxError.exec();
        return false;
    }
    else
    {
        return true;
    }
}

ParameterValidationResult Parameters::CheckValidityNumDensityDistribution(int comboBoxIndex, bool flagVolOrConc) const
{
    ParameterValidationResult result;
    result.isValid = true;

    if (comboBoxIndex == Custom)
        return result;

    //Mirrors Calculate::DiameterRangeSetting()
    double sigma = stdDev;
    double zScore = 3.0 + (nRadius / 100.0);
    double localMinRadius, localMaxRadius;

    if (comboBoxIndex == LogNormal)
    {
        double mu = std::log(meanRadius) - (sigma * sigma / 2.0);
        localMinRadius = std::exp(mu - zScore * sigma);
        localMaxRadius = std::exp(mu + zScore * sigma);
    }
    else //Gaussian
    {
        localMinRadius = meanRadius - zScore * sigma;
        localMaxRadius = meanRadius + zScore * sigma;
    }
    if (localMinRadius < 1e-10) localMinRadius = 1e-10;

    //Mirrors Calculate::SetSphereRadiusAndRefIndex()
    const double volumeConst = 4.0 * M_PI / 3.0;
    const double sqrtTwoPi = std::sqrt(2.0 * M_PI);
    const double twoSigSq = 2.0 * sigma * sigma;
    double stepR = (localMaxRadius - localMinRadius) / (nRadius - 1);

    double totalSphereVolume = 0.0;
    double totalFuncSum = 0.0;
    std::vector<double> funcArray(nRadius, 0.0);

    for (unsigned int i = 0; i < nRadius; i++)
    {
        double r = localMinRadius + i * stepR;

        if (comboBoxIndex == LogNormal)
        {
            double diff = std::log(r) - std::log(meanRadius);
            funcArray[i] = std::exp(-(diff * diff) / twoSigSq) / (r * sigma * sqrtTwoPi);
        }
        else //Gaussian
        {
            double diff = r - meanRadius;
            funcArray[i] = std::exp(-(diff * diff) / twoSigSq) / (sigma * sqrtTwoPi);
        }
        totalSphereVolume += funcArray[i] * volumeConst * pow(r, 3);
        totalFuncSum += funcArray[i];
    }

    double factor = 1e-12;
    if (flagVolOrConc && totalSphereVolume > 0)
    {
        factor = 1e9 * volFraction / totalSphereVolume;
    }
    else if (totalFuncSum > 0)
    {
        factor = sphNumDensity / totalFuncSum;
    }

    double totalRoundedNumDensity = 0.0;
    for (unsigned int i = 0; i < nRadius; i++)
        totalRoundedNumDensity += std::round(funcArray[i] * factor);

    if (totalRoundedNumDensity < 1.0)
    {
        result.isValid = false;
        result.errorMessage = "Volume Fraction is too low! Increase";
        return result;
    }
    return result;
}

//Check the validity of Distribution parameters
bool Parameters::CheckPackingVolume() const
{
    double totalVolume = 0.0;
    for (unsigned int i = 0; i< nRadius; i++)
    {
        double singleSphVolume = (4.0/3.0) * M_PI * pow(radArray[i], 3);
        totalVolume += singleSphVolume * numDensityArray[i];
    }

    ParameterValidationResult check = CheckValidityPackingVolume(totalVolume/1e9);

    if (!check.isValid)
    {
        QMessageBox msgBoxError;
        msgBoxError.setWindowTitle("Error");
        msgBoxError.setIcon(QMessageBox::Critical);
        msgBoxError.setText(check.errorMessage);
        msgBoxError.exec();
        return false;
    }
    else
    {
        return true;
    }
}

// Check packing totalVolume limit for polydisperse distribution
ParameterValidationResult Parameters::CheckValidityPackingVolume(double totalVolume) const
{
    ParameterValidationResult result;
    result.isValid = true;

    if (totalVolume >= M_PI/(3*sqrt(2)))  //Maximum packing factor = PI/(3*sqrt(2))
    {
        result.isValid = false;
        result.errorMessage = "Total sphere volume in 1mm³ exceeds the maximum packing factor. Reduce Concentration.";
        return result;
    }
    else
    {
        return result;
    }
}