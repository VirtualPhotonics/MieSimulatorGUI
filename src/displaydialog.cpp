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

void DisplayDialog::DisplayData(Ui_MainWindow *mainUi, parameters *para)
{
	double margin = (1.0 + mainUi->slider_ConcPercentChange->value() / 200.0);

	//Simulation parameters
    ui->textBrowser_Display->setText("test");
	ui->textBrowser_Display->setText("Simulation parameters:\n");

	if (mainUi->radioButton_MonoDisperse->isChecked())
	{
		ui->textBrowser_Display->append("Distribution: Mono Disperse");
		ui->textBrowser_Display->append("Diameter of spheres: " + QString::number(2.0*para->meanRadius) + " um");
		if (mainUi->radioButton_Conc_mm3->isChecked())
			ui->textBrowser_Display->append("Concentration (Spheres in a volume of 1mm^3): "
				+ QString::number(para->sphNumDensity * margin));
		if (mainUi->radioButton_VolFrac->isChecked())
			ui->textBrowser_Display->append("Volume Fraction (Total sphere volume / 1mm^3): "
				+ QString::number(para->volFraction * margin));
	}

	if (mainUi->radioButton_PolyDisperse->isChecked())
	{
		int currentIndex = mainUi->comboBox_Distribution->currentIndex();
		if (currentIndex == 0) //Log normal distribution
			ui->textBrowser_Display->append("Distribution: Poly Disperse - Log Normal");
		if (currentIndex == 1) //Gaussian distribution
			ui->textBrowser_Display->append("Distribution: Poly Disperse - Gaussian");
		if (currentIndex == 2) //Custom distribution
			ui->textBrowser_Display->append("Distribution: Poly Disperse - Custom");
		ui->textBrowser_Display->append("Number of discrete sphere sizes: "
			+ QString::number(para->nRadius));
		if (currentIndex != 2)
		{
			ui->textBrowser_Display->append("Mean diameter of spheres: "
				+ QString::number(2.0*para->meanRadius) + " um");
			ui->textBrowser_Display->append("Std. deviation: " + QString::number(para->stdDev) + " um\n");
			if (mainUi->radioButton_Conc_mm3->isChecked())
				ui->textBrowser_Display->append("Total Concentration (Spheres in a volume of 1mm^3): " 
					+ QString::number(para->sphNumDensity * margin));
			if (mainUi->radioButton_VolFrac->isChecked())
				ui->textBrowser_Display->append("Volume Fraction (Total sphere volume / 1mm^3): " 
					+ QString::number(para->volFraction * margin));
		}
	}
	ui->textBrowser_Display->append("Wavelength Range: " + QString::number(para->startWavel) + "nm to " 
		+ QString::number(para->endWavel) + "nm in "+ QString::number(para->stepWavel) + "nm steps\n");
	ui->textBrowser_Display->append("Refractive index of the medium: " + QString::number(para->medRef) );

	ui->textBrowser_Display->append("\nSphere Data:");
	ui->textBrowser_Display->append("Dia.(um)\tNum. Den.(in a vol. of 1mm^3)\tRef. Index (real | imag)");
	for (int i = 0; i<para->nRadius; i++)
		ui->textBrowser_Display->append(QString::number(2.0 * para->radArray[i]) + "\t"
			+ QString::number(para->numDensityArray[i] * margin) + "\t"
			+ QString::number(para->scatRefRealArray[i]) + "\t" + QString::number(para->scatRefImagArray[i]) );

     ui->textBrowser_Display->append("\n\nOutput:");
    if (mainUi->radioButton_MonoDisperse->isChecked())
    {
        ui->textBrowser_Display->append("\nSize Parameter:");
        ui->textBrowser_Display->append("WL(nm)\tSize Parameter");
        for (int i = 0; i<para->nWavel; i++)
            ui->textBrowser_Display->append(QString::number(para->wavelArray[i]) + "\t"
                + QString::number(para->SizePara[i]) );        
    }

    ui->textBrowser_Display->append("\nScattering Coefficient (us):");
	ui->textBrowser_Display->append("WL(nm)\tScattering Coefficient (mm^-1)");
	for (int i = 0; i<para->nWavel; i++)
		ui->textBrowser_Display->append(QString::number(para->wavelArray[i]) + "\t" 
			+ QString::number(para->mus[i] * margin) );
	

	ui->textBrowser_Display->append("\ng (Average Cosine of Phase Function):");
	ui->textBrowser_Display->append("WL(nm)\tg (Average Cosine of Phase Function)");
	for (int i = 0; i<para->nWavel; i++)
		ui->textBrowser_Display->append(QString::number(para->wavelArray[i]) + "\t" 
			+ QString::number(para->g[i]) );


	ui->textBrowser_Display->append("\nReduced Scattering Coefficient (us'):");
	ui->textBrowser_Display->append("WL(nm)\tReduced Scattering Coefficient (mm^-1)");
	for (int i = 0; i<para->nWavel; i++)
		ui->textBrowser_Display->append(QString::number(para->wavelArray[i]) + "\t" 
			+ QString::number(para->mus[i] * (1 - para->g[i])*margin) );


	ui->textBrowser_Display->append("\nScattering Cross Section (Csca):");
	ui->textBrowser_Display->append("WL(nm)\tScattering Cross Section (um^2)");
	for (int i = 0; i<para->nWavel; i++)
		ui->textBrowser_Display->append(QString::number(para->wavelArray[i]) + "\t" 
            + QString::number(para->cSca[i]) );

    ui->textBrowser_Display->append("\nExtinction Cross Section (Cext):");
    ui->textBrowser_Display->append("WL(nm)\tExtinction Cross Section (um^2)");
    for (int i = 0; i<para->nWavel; i++)
        ui->textBrowser_Display->append(QString::number(para->wavelArray[i]) + "\t"
            + QString::number(para->cExt[i]) );

    ui->textBrowser_Display->append("\nBackscattering Cross Section (Cback):");
    ui->textBrowser_Display->append("WL(nm)\tBackcattering Cross Section (um^2)");
    for (int i = 0; i<para->nWavel; i++)
        ui->textBrowser_Display->append(QString::number(para->wavelArray[i]) + "\t"
            + QString::number(para->cBack[i]) );    

	ui->textBrowser_Display->append("\nForward and Backward Scattering %");
	ui->textBrowser_Display->append("WL(nm)\tForward %\tBackward %");
	for (int i = 0; i<para->nWavel; i++)
		ui->textBrowser_Display->append(QString::number(para->wavelArray[i]) + "\t" 
			+ QString::number(para->forward[i]) + "\t" + QString::number(para->backward[i]));

	ui->textBrowser_Display->append("\nPhase Function:");
	QVector<double> phaseFunction(para->nTheta), theta(para->nTheta);
    int indexWL = mainUi->slider_WL_PFPolar->value();
	double currentWL = para->startWavel + indexWL*para->stepWavel;
	for (int i = 0; i<para->nTheta; i++)
	{
		theta[i] = 180.0 * i / (para->nTheta - 1);
		if (mainUi->radioButton_PhaseAverage->isChecked())
		{
			phaseFunction[i] = para->phaseFunctionAve[indexWL][i];
		}
		if (mainUi->radioButton_PhasePara->isChecked())
		{
			phaseFunction[i] = para->phaseFunctionPara[indexWL][i];
		}
		if (mainUi->radioButton_PhasePerp->isChecked())
		{
			phaseFunction[i] = para->phaseFunctionPerp[indexWL][i];
		}
	}
	if (mainUi->radioButton_PhaseAverage->isChecked())
		ui->textBrowser_Display->append("Angle(deg)\t Phase Function (Ave) @ Wavelength of " + QString::number(currentWL) + " nm");	
    if (mainUi->radioButton_PhasePara->isChecked())
        ui->textBrowser_Display->append("Angle(deg)\t Phase Function (Para) @ Wavelength of " + QString::number(currentWL) + " nm");
    if (mainUi->radioButton_PhasePerp->isChecked())
        ui->textBrowser_Display->append("Angle(deg)\t Phase Function (Perp) @ Wavelength of " + QString::number(currentWL) + " nm");
	for (int i = 0; i<para->nTheta; i++)
		ui->textBrowser_Display->append(QString::number(theta[i]) + "\t" + QString::number(phaseFunction[i]));
}
