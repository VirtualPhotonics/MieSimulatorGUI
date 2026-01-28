//**********************************************************************
//** "parameters" variable is used to store all variables in a single
//** variable. This allows easy portability.
//**********************************************************************

#include "parameters.h"
#include <QComboBox>
#include <QRadioButton>
#include <QMessageBox>

Parameters::Parameters()
{
}

Parameters::~Parameters()
{
}

//Check the validity of common Parameters
bool Parameters::CheckCommonParameters(QRadioButton *radioButton_MonoDisperse,
                                       QRadioButton *radioButton_NumDen,
                                       QRadioButton *radioButton_VolFrac)
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

    ParameterValidationResult check = CheckValidityCommonParameters(monoDisperseSelection,
                                                                    numDenSelection,
                                                                    volFracSelection);
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
                                                    bool volFracSelection)
{
    ParameterValidationResult result;
    result.isValid = true;

    if ((scatRefReal <= 0.0) || (medRef <= 0.0))
    {
        result.isValid = false;
        result.errorMessage = "Refractive index cannot be zero.";
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
        if (sphNumDensity <= 0.0)
        {
            result.isValid = false;
            result.errorMessage = "Sphere concentration cannot be zero.";
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
            result.errorMessage = "Volume Fraction must not exceed the maximum packing factor.";
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
    }
    return result;
}

//Check the validity of Distribution parameters
bool Parameters::CheckDistributionParameters(QComboBox *comboBox_Distribution)
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

ParameterValidationResult Parameters::CheckValidityDistributionParameters(int comboBoxIndex)
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
    if ((nRadius < 2.0) ||(nRadius >101.0))
    {
        result.isValid = false;
        result.errorMessage = "Number of sphere sizes is out of range! Enter a value between 2 and 101.";
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

//Check the validity of Distribution parameters
bool Parameters::CheckPackingVolume()
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
ParameterValidationResult Parameters::CheckValidityPackingVolume(double totalVolume)
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
