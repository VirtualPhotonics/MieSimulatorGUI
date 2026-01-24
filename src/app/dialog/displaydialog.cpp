/**********************************************************************
** This Displaydialog will allow us to capture current data and display
** it in a textBrowser.
**********************************************************************/

#include "displaydialog.h"
#include "ui_displaydialog.h"

DisplayDialog::DisplayDialog(QWidget *parent) :
	QDialog(parent),
	ui(new Ui::DisplayDialog)
{
	ui->setupUi(this);
}

DisplayDialog::~DisplayDialog()
{
	delete ui;
}

void DisplayDialog::on_pushButton_Close_clicked()
{
	this->close();
}

void DisplayDialog::DisplayData(QRadioButton *radioButton_MonoDisperse,
                                QRadioButton *radioButton_PolyDisperse,
                                QRadioButton *radioButton_NumDen,
                                QRadioButton *radioButton_VolFrac,
                                QComboBox *comboBox_Distribution,
                                QSlider *slider_ConcPercentChange,
                                QSlider *slider_WL_PFPolar,
                                QCheckBox *checkBox_PhasePolarAve,
                                QCheckBox *checkBox_PhasePolarPara,
                                QCheckBox *checkBox_PhasePolarPerp,
                                Parameters *para)
{
    double margin = (1.0 + slider_ConcPercentChange->value() / 200.0);

	//Simulation parameters
    ui->textBrowser_Display->setText("Simulation parameters:\n");

    if (radioButton_MonoDisperse->isChecked())
	{
        ui->textBrowser_Display->append("Distribution: Mono Disperse");
		ui->textBrowser_Display->append("Diameter of spheres: " + QString::number(2.0*para->meanRadius) + " um");
        if (radioButton_NumDen->isChecked())
			ui->textBrowser_Display->append("Concentration (Spheres in a volume of 1mm^3): "
				+ QString::number(para->sphNumDensity * margin));
        if (radioButton_VolFrac->isChecked())
			ui->textBrowser_Display->append("Volume Fraction (Total sphere volume / 1mm^3): "
				+ QString::number(para->volFraction * margin));
	}

    if (radioButton_PolyDisperse->isChecked())
	{
        int currentIndex = comboBox_Distribution->currentIndex();
        if (currentIndex == para->LogNormal) //Log normal distribution
			ui->textBrowser_Display->append("Distribution: Poly Disperse - Log Normal");
        if (currentIndex == para->Gaussian) //Gaussian distribution
			ui->textBrowser_Display->append("Distribution: Poly Disperse - Gaussian");
        if (currentIndex == para->Custom) //Custom distribution
			ui->textBrowser_Display->append("Distribution: Poly Disperse - Custom");
		ui->textBrowser_Display->append("Number of discrete sphere sizes: "
            + QString::number(para->nRadius));
        if (currentIndex != para->Custom)
		{
			ui->textBrowser_Display->append("Mean diameter of spheres: "
				+ QString::number(2.0*para->meanRadius) + " um");
			ui->textBrowser_Display->append("Std. deviation: " + QString::number(para->stdDev) + " um\n");
            if (radioButton_NumDen->isChecked())
				ui->textBrowser_Display->append("Total Concentration (Spheres in a volume of 1mm^3): " 
					+ QString::number(para->sphNumDensity * margin));
            if (radioButton_VolFrac->isChecked())
				ui->textBrowser_Display->append("Volume Fraction (Total sphere volume / 1mm^3): " 
					+ QString::number(para->volFraction * margin));
		}
	}
    ui->textBrowser_Display->append("Wavelength Range: " + QString::number(para->startWavel) + "nm to "
        + QString::number(para->endWavel) + "nm in "+ QString::number(para->stepWavel) + "nm steps\n");

	ui->textBrowser_Display->append("\nSphere Data:");
    ui->textBrowser_Display->append("Dia.(um)\tNum. Den.(in a vol. of 1mm^3)\tRef."
                                    " index of sphere (real | imag)\t Ref. index of medium");
    for (unsigned int i = 0; i<para->nRadius; i++)
		ui->textBrowser_Display->append(QString::number(2.0 * para->radArray[i]) + "\t"
			+ QString::number(para->numDensityArray[i] * margin) + "\t"
            + QString::number(para->scatRefRealArray[i]) + "\t"
            + QString::number(para->scatRefImagArray[i]) + "\t"
            + QString::number(para->medRefArray[i]));

     ui->textBrowser_Display->append("\n\nOutput:");
    if (radioButton_MonoDisperse->isChecked())
    {
        ui->textBrowser_Display->append("\nSize Parameter (2*pi*R/lambda):");
        ui->textBrowser_Display->append("WL(nm)\t2*pi*R/lambda");
        for (unsigned int i = 0; i<para->nWavel; i++)
            ui->textBrowser_Display->append(QString::number(para->wavelArray[i]) + "\t" +
                                            QString::number(para->sizePara[i]) );
    }

    ui->textBrowser_Display->append("\nScattering Coefficient (us), g (Average Cosine of Phase Function) and Reduced Scattering Coefficient (us'):");
    ui->textBrowser_Display->append("WL(nm)\tus(mm^-1)\tg\tus'(mm^-1)");
    for (unsigned int i = 0; i<para->nWavel; i++)
        ui->textBrowser_Display->append(QString::number(para->wavelArray[i]) + "\t" +
                                        QString::number(para->mus[i] * margin) + "\t" +
                                        QString::number(para->g[i]) + "\t" +
                                        QString::number(para->mus[i] * (1 - para->g[i])*margin));

    ui->textBrowser_Display->append("\nScattering (Csca), Extinction (Cext) and Backscattering (Cback) Cross Sections:");
    ui->textBrowser_Display->append("WL(nm)\tCsca(um^2)\tCext(um^2)\tCback(um^2)");
    for (unsigned int i = 0; i<para->nWavel; i++)
        ui->textBrowser_Display->append(QString::number(para->wavelArray[i]) + "\t" +
                                        QString::number(para->cSca[i]) + "\t" +
                                        QString::number(para->cExt[i]) + "\t" +
                                        QString::number(para->cBack[i]));

    ui->textBrowser_Display->append("\nForward and Backward Scattering percentages:");
	ui->textBrowser_Display->append("WL(nm)\tForward %\tBackward %");
    for (unsigned int i = 0; i<para->nWavel; i++)
        ui->textBrowser_Display->append(QString::number(para->wavelArray[i]) + "\t" +
                                        QString::number(para->forward[i]) + "\t" +
                                        QString::number(para->backward[i]));

    int indexWL = slider_WL_PFPolar->value();
	double currentWL = para->startWavel + indexWL*para->stepWavel;

    ui->textBrowser_Display->append("\nPhase Function @ Wavelength of " + QString::number(currentWL) + " nm:");
    if (checkBox_PhasePolarAve->isChecked())
    {
        ui->textBrowser_Display->append("\nAngle(deg)\tPhase Function (Ave)");
        for (unsigned int i = 0; i<para->nTheta; i++)
            ui->textBrowser_Display->append(QString::number(180.0 * i / (para->nTheta - 1)) + "\t" +
                                            QString::number(para->phaseFunctionAve[indexWL][i]));
    }
    if (checkBox_PhasePolarPara->isChecked())
    {
        ui->textBrowser_Display->append("\nAngle(deg)\tPhase Function (Para)");
        for (unsigned int i = 0; i<para->nTheta; i++)
            ui->textBrowser_Display->append(QString::number(180.0 * i / (para->nTheta - 1)) + "\t" +
                                            QString::number(para->phaseFunctionPara[indexWL][i]));
    }
    if (checkBox_PhasePolarPerp->isChecked())
    {
        ui->textBrowser_Display->append("\nAngle(deg)\tPhase Function (Perp)");
        for (unsigned int i = 0; i<para->nTheta; i++)
            ui->textBrowser_Display->append(QString::number(180.0 * i / (para->nTheta - 1)) + "\t" +
                                            QString::number(para->phaseFunctionPerp[indexWL][i]));
    }
}
