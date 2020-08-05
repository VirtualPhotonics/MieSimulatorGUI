//***********************************************************************
//** This Displaydialog will allow us to capture current data and display
//** it in a textBrowser.
//***********************************************************************

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

void DisplayDialog::on_pushButton_close_clicked()
{
	this->close();
}

void DisplayDialog::DisplayData(Ui_MainWindow *mainUi, parameters *para)
{
    double margin = (1.0 + mainUi->slider_concPercentChange->value() / 200.0);

	//Simulation parameters
    ui->textBrowser_display->setText("test");
    ui->textBrowser_display->setText("Simulation parameters:\n");

    if (mainUi->radioButton_monoDisperse->isChecked())
	{
        ui->textBrowser_display->append("Distribution: Mono Disperse");
        ui->textBrowser_display->append("Diameter of spheres: " + QString::number(2.0*para->meanRadius) + " um");
        if (mainUi->radioButton_conc_mm3->isChecked())
            ui->textBrowser_display->append("Concentration (Spheres in a volume of 1mm^3): "
				+ QString::number(para->sphNumDensity * margin));
        if (mainUi->radioButton_volFrac->isChecked())
            ui->textBrowser_display->append("Volume Fraction (Total sphere volume / 1mm^3): "
				+ QString::number(para->volFraction * margin));
	}

    if (mainUi->radioButton_polyDisperse->isChecked())
	{
        int currentIndex = mainUi->comboBox_distribution->currentIndex();
		if (currentIndex == 0) //Log normal distribution
            ui->textBrowser_display->append("Distribution: Poly Disperse - Log Normal");
		if (currentIndex == 1) //Gaussian distribution
            ui->textBrowser_display->append("Distribution: Poly Disperse - Gaussian");
		if (currentIndex == 2) //Custom distribution
            ui->textBrowser_display->append("Distribution: Poly Disperse - Custom");
        ui->textBrowser_display->append("Number of discrete sphere sizes: "
			+ QString::number(para->nRadius));
		if (currentIndex != 2)
		{
            ui->textBrowser_display->append("Mean diameter of spheres: "
				+ QString::number(2.0*para->meanRadius) + " um");
            ui->textBrowser_display->append("Std. deviation: " + QString::number(para->stdDev) + " um\n");
            if (mainUi->radioButton_conc_mm3->isChecked())
                ui->textBrowser_display->append("Total Concentration (Spheres in a volume of 1mm^3): "
					+ QString::number(para->sphNumDensity * margin));
            if (mainUi->radioButton_volFrac->isChecked())
                ui->textBrowser_display->append("Volume Fraction (Total sphere volume / 1mm^3): "
					+ QString::number(para->volFraction * margin));
		}
	}
    ui->textBrowser_display->append("Wavelength Range: " + QString::number(para->startWavel) + "nm to "
		+ QString::number(para->endWavel) + "nm in "+ QString::number(para->stepWavel) + "nm steps\n");
    ui->textBrowser_display->append("Refractive index of the medium: " + QString::number(para->medRef) );

    ui->textBrowser_display->append("\nSphere Data:");
    ui->textBrowser_display->append("Dia.(um)\tNum. Den.(in a vol. of 1mm^3)\tRef. Index (real | imag)");
    for (unsigned int i = 0; i<para->nRadius; i++)
        ui->textBrowser_display->append(QString::number(2.0 * para->radArray[i]) + "\t"
			+ QString::number(para->numDensityArray[i] * margin) + "\t"
			+ QString::number(para->scatRefRealArray[i]) + "\t" + QString::number(para->scatRefImagArray[i]) );

     ui->textBrowser_display->append("\n\nOutput:");
    if (mainUi->radioButton_monoDisperse->isChecked())
    {
        ui->textBrowser_display->append("\nSize Parameter:");
        ui->textBrowser_display->append("WL(nm)\t2*pi*R/lambda");
        for (unsigned int i = 0; i<para->nWavel; i++)
            ui->textBrowser_display->append(QString::number(para->wavelArray[i]) + "\t" +
                                            QString::number(para->SizePara[i]) );
    }

    ui->textBrowser_display->append("\nScattering Coefficient, g (Average Cosine of Phase Function) and Reduced Scattering Coefficient:");
    ui->textBrowser_display->append("WL(nm)\tus(mm^-1)\tg\tus'(mm^-1)");
    for (unsigned int i = 0; i<para->nWavel; i++)
        ui->textBrowser_display->append(QString::number(para->wavelArray[i]) + "\t" +
                                        QString::number(para->mus[i] * margin) + "\t" +
                                        QString::number(para->g[i]) + "\t" +
                                        QString::number(para->mus[i] * (1 - para->g[i])*margin));

    ui->textBrowser_display->append("\nScattering (Csca), Extinction (Cext) and Backscattering (Cback) Cross Sections:");
    ui->textBrowser_display->append("WL(nm)\tCsca(um^2)\tCext(um^2)\tCback(um^2)");
    for (unsigned int i = 0; i<para->nWavel; i++)
        ui->textBrowser_display->append(QString::number(para->wavelArray[i]) + "\t" +
                                        QString::number(para->cSca[i]) + "\t" +
                                        QString::number(para->cExt[i]) + "\t" +
                                        QString::number(para->cBack[i]));

    ui->textBrowser_display->append("\nForward and Backward Scattering percentages:");
    ui->textBrowser_display->append("WL(nm)\tForward %\tBackward %");
    for (unsigned int i = 0; i<para->nWavel; i++)
        ui->textBrowser_display->append(QString::number(para->wavelArray[i]) + "\t" +
                                        QString::number(para->forward[i]) + "\t" +
                                        QString::number(para->backward[i]));

    int indexWL = mainUi->slider_pFunctionPolar_wavel->value();
	double currentWL = para->startWavel + indexWL*para->stepWavel;

    ui->textBrowser_display->append("\nPhase Function @ Wavelength of " + QString::number(currentWL) + " nm:");
    if (mainUi->radioButton_pFunction_ave->isChecked())
        ui->textBrowser_display->append("Angle(deg)\tPhase Function (Ave)");
    if (mainUi->radioButton_pFunction_para->isChecked())
        ui->textBrowser_display->append("Angle(deg)\tPhase Function (Para)");
    if (mainUi->radioButton_pFunction_perp->isChecked())
        ui->textBrowser_display->append("Angle(deg)\tPhase Function (Perp)");
    for (unsigned int i = 0; i<para->nTheta; i++)
    {
        if (mainUi->radioButton_pFunction_ave->isChecked())
            ui->textBrowser_display->append(QString::number(180.0 * i / (para->nTheta - 1)) + "\t" +
                                            QString::number(para->phaseFunctionAve[indexWL][i]));
        if (mainUi->radioButton_pFunction_para->isChecked())
            ui->textBrowser_display->append(QString::number(180.0 * i / (para->nTheta - 1)) + "\t" +
                                            QString::number(para->phaseFunctionPara[indexWL][i]));
        if (mainUi->radioButton_pFunction_perp->isChecked())
            ui->textBrowser_display->append(QString::number(180.0 * i / (para->nTheta - 1)) + "\t" +
                                            QString::number(para->phaseFunctionPerp[indexWL][i]));
    }
}
