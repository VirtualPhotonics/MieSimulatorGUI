/**********************************************************************
** This Optionsdialog will allow users to select and save current data.
**********************************************************************/

#include "optionsdialog.h"

OptionsDialog::OptionsDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::OptionsDialog)
{    
    flagScatPara = false;
    flagPhaseFunction = false;
    flagS1 = false;
    flagS2 = false;
    ui->setupUi(this);
}

OptionsDialog::~OptionsDialog()
{
    delete ui;
}

void OptionsDialog::on_pushButton_ScatPara_clicked()
{
    flagScatPara = true;
    this->close();
}

void OptionsDialog::on_pushButton_PhaseFunction_clicked()
{
    flagPhaseFunction = true;
    this->close();
}

void OptionsDialog::on_pushButton_S1_clicked()
{
    flagS1 = true;
    this->close();
}

void OptionsDialog::on_pushButton_S2_clicked()
{
    flagS2 = true;
    this->close();
}

void OptionsDialog::on_pushButton_Cancel_clicked()
{
    this->close();
}

void OptionsDialog::SaveData(Ui_MainWindow *ui, parameters *para)
{
    setModal(true);
    exec();

    double margin = (1.0 + ui->slider_ConcPercentChange->value() /200.0);
    if (flagScatPara)
    {
        QString fileName = QFileDialog::getSaveFileName(this,
                 tr("Save Data"), "Mie_ScatteringParameters",
                 tr("Text File (*.txt);;All Files (*)"));

        if (fileName.isEmpty())
            return;
        else
        {
            SaveScatPara(ui, para,fileName,margin);
            RememberLastDirectory(fileName);
        }
    }
    if (flagPhaseFunction)
    {
        QString fileName = QFileDialog::getSaveFileName(this,
                 tr("Save Data"), "Mie_PhaseFunctionData",
                 tr("Text File (*.txt);;All Files (*)"));

        if (fileName.isEmpty())
            return;
        else
        {
            SavePhaseFunction(ui, para,fileName);
            RememberLastDirectory(fileName);
        }
    }
    if (flagS1)
    {
        QString fileName = QFileDialog::getSaveFileName(this,
                 tr("Save Data"), "Mie_S1Data",
                 tr("Text File (*.txt);;All Files (*)"));

        if (fileName.isEmpty())
            return;
        else
        {
            SaveS1(para,fileName);
            RememberLastDirectory(fileName);
        }
    }
    if (flagS2)
    {
        QString fileName = QFileDialog::getSaveFileName(this,
                 tr("Save Data"), "Mie_S2Data",
                 tr("Text File (*.txt);;All Files (*)"));

        if (fileName.isEmpty())
            return;
        else
        {
            SaveS2(para,fileName);
            RememberLastDirectory(fileName);
        }
    }

}

void OptionsDialog::SaveScatPara(Ui_MainWindow *ui, parameters *para, QString fileName, double margin)
{    

    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        QMessageBox::information(this, tr("Unable to open file"),file.errorString());
        return;
    }

    /***************Writing Data***********************/
    QTextStream out(&file);

    //Simulation parameters
    out << "Simulation parameters:\n\n";

    if(ui->radioButton_MonoDisperse->isChecked())
    {
        out << "Distribution: Mono Disperse \n";
        out << "Diameter of spheres: " << 2.0*para->meanRadius <<" um\n";
        if (ui->radioButton_NumDen->isChecked())
            out << "Concentration (Spheres in a volume of 1mm^3): " << para->sphNumDensity * margin<<"\n";
        if (ui->radioButton_VolFrac->isChecked())
            out << "Volume Fraction (Total sphere volume / 1mm^3): " << para->volFraction * margin<<"\n";
    }

    if(ui->radioButton_PolyDisperse->isChecked())
    {
        int currentIndex = ui->comboBox_Distribution->currentIndex();
        if (currentIndex == 0) //Log normal distribution
            out << "Distribution: Poly Disperse - Log Normal \n";
        if (currentIndex == 1) //Gaussian distribution
            out << "Distribution: Poly Disperse - Gaussian \n";
        if (currentIndex == 2) //Custom distribution
            out << "Distribution: Poly Disperse - Custom \n";
        out << "Number of discrete sphere sizes: " << para->nRadius <<"\n";
        if (currentIndex != 2)
        {
            out << "Mean diameter of spheres: " << 2.0*para->meanRadius <<" um\n";
            out << "Std. deviation: " << para->stdDev <<" um\n";
            if (ui->radioButton_NumDen->isChecked())
                out << "Total Concentration (Spheres in a volume of 1mm^3): " << para->sphNumDensity * margin <<"\n";
            if (ui->radioButton_VolFrac->isChecked())
                out << "Volume Fraction (Total sphere volume / 1mm^3): " << para->volFraction * margin <<"\n";
        }
    }

    out << "Wavelength Range: " << para->startWavel << "nm to "  << para->endWavel << "nm in "
        << para->stepWavel <<"nm steps\n";
    out << "\nRefractive index of the medium: " << para->medRef <<"\n\n";

    out << "Sphere Data:\n";
    out << "Dia.(um)\tNum. Den.(in a vol. of 1mm^3)\tRef. Index (real | imag)\n";
    for (unsigned int i=0; i<para->nRadius; i++)
        out << 2.0 * para->radArray[i] <<"\t" << para->numDensityArray[i]*margin<<"\t"
            <<para->scatRefRealArray[i]<<"\t" << para->scatRefImagArray[i] <<"\n";

    out << "\n\nOutput:\n\n";
    if (ui->radioButton_MonoDisperse->isChecked())
    {
        out << "Size Parameter (2*pi*R/lambda):\n";
        out << "WL(nm)\t2*pi*R/lambda\n";
        for (unsigned int i=0; i<para->nWavel; i++)
            out << para->wavelArray[i] <<"\t" << para->SizePara[i] <<"\n";
    }

    out << "\nScattering Coefficient (us), g and Reduced Scattering Coefficient (us'):\n";
    out << "WL(nm)\tus(mm^-1)\tg\tus'(mm^-1)\n";
    for (unsigned int i=0; i<para->nWavel; i++)
        out << para->wavelArray[i] <<"\t" << para->mus[i]* margin << "\t" << para->g[i] <<"\t" << para->mus[i]*(1-para->g[i])* margin <<"\n";
    out << "\n";

    out << "Scattering (Csca), Extinction (Cext) and Backscattering (Cback) Cross Sections:\n";
    out << "WL(nm)\tCsca(um^2)\tCext(um^2)\tCback(um^2)\n";
    for (unsigned int i=0; i<para->nWavel; i++)
        out << para->wavelArray[i] <<"\t" << para->cSca[i] <<"\t" << para->cExt[i] <<"\t" << para->cBack[i] <<"\n";
    out << "\n";

    out << "Forward and Backward Scattering percentages:\n";
    out << "WL(nm)\tForward %\tBackward %\n";
    for (unsigned int i=0; i<para->nWavel; i++)
        out << para->wavelArray[i] <<"\t" << para->forward[i] <<"\t" << para->backward[i] <<"\n";
    out << "\n";
}

void OptionsDialog::SavePhaseFunction(Ui_MainWindow *ui, parameters *para, QString fileName)
{

    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        QMessageBox::information(this, tr("Unable to open file"),file.errorString());
        return;
    }

    /***************Writing Data***********************/
    QTextStream out(&file);

    out << "Phase Function ";
    if (ui->radioButton_PhaseAverage->isChecked())
        out << "(Ave): ";
    if (ui->radioButton_PhasePara->isChecked())
         out << "(Para): ";
    if (ui->radioButton_PhasePerp->isChecked())
         out << "(Perp): ";
     out << "\n\tWL(nm)-->\nAngle(deg)";
     for (unsigned int i=0; i<para->nWavel; i++)
        out << para->wavelArray[i] <<"\t";
     out << "\n";

    for (int j=0; j<static_cast<int>(para->nTheta); j++)
    {
        out << 180.0 * j /(para->nTheta-1) <<"\t";
        for (unsigned int i=0; i<para->nWavel; i++)
        {
            if (ui->radioButton_PhaseAverage->isChecked())
                out << para->phaseFunctionAve[i][j] <<"\t";
            if (ui->radioButton_PhasePara->isChecked())
                out << para->phaseFunctionPara[i][j];
            if (ui->radioButton_PhasePerp->isChecked())
                out << para->phaseFunctionPerp[i][j];
        }
        out << "\n";
    }
}

void OptionsDialog::SaveS1(parameters *para, QString fileName)
{

    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        QMessageBox::information(this, tr("Unable to open file"),file.errorString());
        return;
    }

    /***************Writing Data***********************/
    QTextStream out(&file);

    out << "S1 (Real, Imaginary):";
    out << "\n\tWL(nm)-->\nAngle(deg)";
    for (unsigned int i=0; i<para->nWavel; i++)
       out << para->wavelArray[i] <<"\t\t";
    out << "\n";

    for (int j=0; j<static_cast<int>(para->nTheta); j++)
    {
        out << 180.0 * j /(para->nTheta-1) <<"\t";
        for (unsigned int i=0; i<para->nWavel; i++)
            out << para->S1[i][j].real() <<"\t" << para->S1[i][j].imag() <<"\t";
         out << "\n";
    }
}

void OptionsDialog::SaveS2(parameters *para, QString fileName)
{

    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        QMessageBox::information(this, tr("Unable to open file"),file.errorString());
        return;
    }

    /***************Writing Data***********************/
    QTextStream out(&file);

    out << "S2 (Real, Imaginary):";
    out << "\n\tWL(nm)-->\nAngle(deg)";
    for (unsigned int i=0; i<para->nWavel; i++)
       out << para->wavelArray[i] <<"\t\t";
    out << "\n";

    for (int j=0; j<static_cast<int>(para->nTheta); j++)
    {
        out << 180.0 * j /(para->nTheta-1) <<"\t";
        for (unsigned int i=0; i<para->nWavel; i++)
            out << para->S2[i][j].real() <<"\t" << para->S2[i][j].imag() <<"\t";
         out << "\n";
    }
}

void OptionsDialog::RememberLastDirectory(QString fileName)
{
    int pos = fileName.lastIndexOf('/');
    QDir::setCurrent(fileName.left(pos));
}
