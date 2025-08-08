//**********************************************************************
//** "parameters" variable is used to store all variables in a single
//** variable. This allows easy portability.
//**********************************************************************

#include "parameters.h"
#include "qcombobox.h"
#include "qradiobutton.h"
#include <QMessageBox>

parameters::parameters(void)
{
}

parameters::~parameters(void)
{
}

//Check the validity of common parameters
bool parameters::CheckCommonParameters(QRadioButton *radioButton_MonoDisperse,
                                       QRadioButton *radioButton_NumDen,
                                       QRadioButton *radioButton_VolFrac)
{
    QMessageBox msgBox, msgBoxWarn;
    msgBox.setWindowTitle("Error");
    msgBox.setIcon(QMessageBox::Critical);
    msgBoxWarn.setWindowTitle("Warning");
    msgBoxWarn.setIcon(QMessageBox::Warning);

    if ((scatRefReal <= 0.0) || (medRef <= 0.0))
    {
        msgBox.setText("Refractive index cannot be zero");
        msgBox.exec();
        return 1;
    }
    if ((scatRefReal/medRef == 1.0))
    {
        msgBox.setText("Relative refractive index cannot be 1.0");
        msgBox.exec();
        return 1;
    }
    double m = scatRefReal / medRef;
    if ((m < 0.05) || (m > 5.0))
    {
        msgBoxWarn.setText("Unrealistic relative refractive index.");
        msgBoxWarn.setInformativeText("Check sphere and medium refractive index values");
        msgBoxWarn.exec();
        return 1;
    }
    if (scatRefImag >= 5.0)
    {
        msgBox.setText("Imaginary refractive index must be negative or less than 5.0.");
        msgBox.exec();
        return 1;
    }
    //Avoid zeros
    if ((startWavel <= 0.0) || (endWavel <= 0.0))
    {
        msgBox.setText("The starting or ending wavelength cannot be zero");
        msgBox.exec();
        return 1;
    }
    //Avoid zeros
    if (stepWavel <= 0.0)
    {
        msgBox.setText("Wavelength step cannot be zero");
        msgBox.exec();
        return 1;
    }
    if (startWavel < 50.0)
    {
        msgBox.setText("Minimum wavlength is 50nm");
        msgBox.exec();
        return 1;
    }
    if (endWavel > 3000.0)
    {
        msgBox.setText("Maximum wavlength is 3000nm");
        msgBox.exec();
        return 1;
    }
    if (endWavel - startWavel < 0.0)
    {
        msgBox.setText("The starting wavelength is greater than the ending wavelength");
        msgBox.exec();
        return 1;
    }
    if (radioButton_NumDen->isChecked())
    {
        if (sphNumDensity <= 0.0)
        {
            msgBox.setText("Sphere concentration cannot be zero");
            msgBox.exec();
            return 1;
        }
    }
    if (radioButton_VolFrac->isChecked())
    {
        if (volFraction <= 0.0)
        {
            msgBox.setText("Volume Fraction cannot be zero");
            msgBox.exec();
            return 1;
        }

        if (volFraction >= 1.0)
        {
            msgBox.setText("Volume Fraction cannot exceed 1.0");
            msgBox.exec();
            return 1;
        }
    }
    if ((meanRadius < 0.00005) ||(meanRadius >150))
    {
        msgBox.setText("Diameter is out of range!");
        msgBox.setInformativeText("Enter a value between 0.0001μm and 300μm");
        msgBox.exec();
        return 1;
    }
    if (radioButton_MonoDisperse->isChecked())
    {
        if (radioButton_NumDen->isChecked())
        {
            double volume = 4.0 * M_PI *meanRadius * meanRadius * meanRadius / 3.0;
            if (sphNumDensity*volume >= 1e9)
            {
                msgBoxWarn.setText("Concentration x Sphere Volume exceeds 1mm³.");
                msgBoxWarn.setInformativeText("Reduce Concentration.");
                msgBoxWarn.exec();
                return 1;
            }
        }
    }
    return 0;
}

bool parameters::CheckDistributionParameters(QComboBox *comboBox_Distribution)
{
    QMessageBox msgBox;
    msgBox.setWindowTitle("Error");

    if (stdDev == 0.0)
    {
        msgBox.setText("Standard Deviation is zero.");
        msgBox.setInformativeText("Use 'Mono Disperse'.");
        msgBox.exec();
        return 1;
    }
    if(comboBox_Distribution->currentIndex() == 0)
    {
        if (stdDev > 3.0)
        {
            msgBox.setIcon(QMessageBox::Critical);
            msgBox.setText("Large standard deviation provides an abnormal Log Normal distribution.");
            msgBox.setInformativeText("The limit was set to 3.0μm.");
            msgBox.exec();
            return 1;
        }
        if (stdDev < 1e-5)
        {
            msgBox.setIcon(QMessageBox::Critical);
            msgBox.setText("The standard deviation is too small.");
            msgBox.setInformativeText("Use 'Mono Disperse'.");
            msgBox.exec();
            return 1;
        }
    }

    if(comboBox_Distribution->currentIndex() == 1)
    {
        if (stdDev > 50.0)
        {
            msgBox.setIcon(QMessageBox::Critical);
            msgBox.setText("The standard deviation is too large.");
            msgBox.setInformativeText("The limit was set to 50.0μm.");
            msgBox.exec();
            return 1;
        }
        if (stdDev < 1e-8)
        {
            msgBox.setIcon(QMessageBox::Critical);
            msgBox.setText("The standard deviation is too small.");
            msgBox.setInformativeText("Use 'Mono Disperse'.");
            msgBox.exec();
            return 1;
        }
    }
    if (nRadius == 1)
    {
        msgBox.setText("Discrete sphere size is 1. ");
        msgBox.setInformativeText("Use 'Mono Disperse'.");
        msgBox.exec();
        return 1;
    }
    if ((nRadius < 2.0) ||(nRadius >101.0))
    {
        msgBox.setText("Number of sphere sizes is out of range.");
        msgBox.setInformativeText("Enter a value between 2 and 101.");
        msgBox.exec();
        return 1;
    }
    if ((meanRadius < 0.0005) ||(meanRadius >25))
    {
        msgBox.setText("Diameter is out of range.");
        msgBox.setInformativeText("Enter a value between 0.001μm and 50μm.");
        msgBox.exec();
        return 1;
    }
    if (stdDev/meanRadius < 1.999e-5)
    {
        msgBox.setIcon(QMessageBox::Critical);
        msgBox.setText("Standard deviation to mean diameter ratio is smaller than 1e-5.");
        msgBox.setInformativeText("Use 'Mono Disperse'");
        msgBox.exec();
        return 1;
    }
    return 0;
}
