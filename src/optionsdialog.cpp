/**********************************************************************
** This Optionsdialog will allow users to select and save current data.
**********************************************************************/

#include "optionsdialog.h"

OptionsDialog::OptionsDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::OptionsDialog)
{    
    ui->setupUi(this);    
    cSca = false;
    cExt = false;
    cBack = false;
    scatCoeff = false;
    redScatCoeff = false;
    fitPara = false;
    g = false;
    phaseFunc = false;
    S1 = false;
    S2 = false;
    forward = false;
    backward = false;
    ApplyAllEnableDisable();
}

OptionsDialog::~OptionsDialog()
{
    delete ui;
}

void OptionsDialog::on_pushButton_Apply_clicked()
{
    this->close();
}

void OptionsDialog::on_pushButton_Cancel_clicked()
{
    cSca = false;
    cExt = false;
    cBack = false;
    scatCoeff = false;
    redScatCoeff = false;
    fitPara = false;
    g = false;
    phaseFunc = false;
    S1 = false;
    S2 = false;
    forward = false;
    backward = false;
    this->close();
}

void OptionsDialog::CheckSelectAll()
{
    if (cSca && scatCoeff && redScatCoeff && fitPara
            && g && phaseFunc && S1 && S2 && forward
            && backward && cExt && cBack)
        ui->checkBox_SelectAll->setChecked(true);
    else
        ui->checkBox_SelectAll->setChecked(false);
}

void OptionsDialog::ApplyAllEnableDisable()
{
    if (cSca || scatCoeff || redScatCoeff || fitPara
            || g || phaseFunc || S1 || S2 || forward
            || backward ||cExt || cBack)
        ui->pushButton_Apply->setEnabled(true);
    else
        ui->pushButton_Apply->setEnabled(false);
}

void OptionsDialog::SaveData(Ui_MainWindow *ui, parameters *para)
{
    bool flagSave;

    if (para->bMie == 0.0)
        DisableFitPara();
    else
       EnableFitPara();

    setModal(true);
    exec();

    if (cSca || scatCoeff || redScatCoeff || fitPara
            || g || phaseFunc || S1 || S2 || forward
            || backward ||cExt || cBack)
        flagSave = true;
    else
        flagSave = false;

    if (flagSave)
    {
        QString fileName = QFileDialog::getSaveFileName(this,
                 tr("Save Data"), "MieData",
                 tr("Text File (*.txt);;All Files (*)"));

        if (fileName.isEmpty())
            return;
        else
        {
            QFile file(fileName);
            if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
            {
                QMessageBox::information(this, tr("Unable to open file"),file.errorString());
                return;
            }

            /***************Writing Data***********************/
            double margin = (1.0 + ui->slider_ConcPercentChange->value() /200.0);
            QTextStream out(&file);

            //Simulation parameters
            out << "Simulation parameters:\n\n";

            if(ui->radioButton_MonoDisperse->isChecked())
            {
                out << "Distribution: Mono Disperse \n";
                out << "Diameter of spheres: " << 2.0*para->meanRadius <<" um\n";
                if (ui->radioButton_Conc_mm3->isChecked())
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
                    if (ui->radioButton_Conc_mm3->isChecked())
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
            for (int i=0; i<para->nRadius; i++)
                out << 2.0 * para->radArray[i] <<"\t" << para->numDensityArray[i]*margin<<"\t"
                    <<para->scatRefRealArray[i]<<"\t" << para->scatRefImagArray[i] <<"\n";            

            if (scatCoeff)
            {
                out << "Scattering Coefficient (us):\n";
                out << "WL(nm)\tScattering Coefficient (mm^-1)\n";
                for (int i=0; i<para->nWavel; i++)
                    out << para->wavelArray[i] <<"\t" << para->mus[i]* margin <<"\n";
                out << "\n";
            }

            if (g)
            {
                out << "g (Average Cosine of Phase Function):\n";
                out << "WL(nm)\tg (Average Cosine of Phase Function)\n";
                for (int i=0; i<para->nWavel; i++)
                    out << para->wavelArray[i] <<"\t" << para->g[i] <<"\n";
                out << "\n";
            }

            if (redScatCoeff)
            {
                out << "Reduced Scattering Coefficient (us'):\n";
                out << "WL(nm)\tReduced Scattering Coefficient (mm^-1)\n";
                for (int i=0; i<para->nWavel; i++)
                    out << para->wavelArray[i] <<"\t" << para->mus[i]*(1-para->g[i])* margin <<"\n";
                out << "\n";
            }

            if (cSca)
            {
                out << "Scattering Cross Section (Csca):\n";
                out << "WL(nm)\tScattering Cross Section (um^2)\n";
                for (int i=0; i<para->nWavel; i++)
                    out << para->wavelArray[i] <<"\t" << para->cSca[i] <<"\n";
                out << "\n";
            }

            if (cExt)
            {
                out << "Extinction Cross Section (Cext):\n";
                out << "WL(nm)\tExtinction Cross Section (um^2)\n";
                for (int i=0; i<para->nWavel; i++)
                    out << para->wavelArray[i] <<"\t" << para->cExt[i] <<"\n";
                out << "\n";
            }

            if (cBack)
            {
                out << "Backscattering Cross Section (Cback):\n";
                out << "WL(nm)\tBackscattering Cross Section (um^2)\n";
                for (int i=0; i<para->nWavel; i++)
                    out << para->wavelArray[i] <<"\t" << para->cBack[i] <<"\n";
                out << "\n";
            }

            if (fitPara)
            {
                out << "us' Power Law Fitting Parmeters\n";
                out << " fRay: " << para->fRay <<"\n";
                out << " bMie: " << para->bMie <<"\n";
                out << " A: " << para->fittedA <<"\n";
                out << " Mean Square Error: " << para->muspFittingError <<"\n";
                out << "\n";
            }

            if (forward)
            {
                out << "Forward Scattering:\n";
                out << "WL(nm)\tForward Scattering %\n";
                for (int i=0; i<para->nWavel; i++)
                    out << para->wavelArray[i] <<"\t" << para->forward[i] <<"\n";
                out << "\n";
            }

            if (backward)
            {
                out << "Backward Scattering:\n";
                out << "WL(nm)\tBackward Scattering %\n";
                for (int i=0; i<para->nWavel; i++)
                    out << para->wavelArray[i] <<"\t" << para->backward[i] <<"\n";
                out << "\n";
            }            

            if (phaseFunc)
            {
                out << "Phase Function:\n";
                QVector<double> phaseFunction(2*para->nTheta-1), theta(2*para->nTheta-1);
                int indexWL = ui->slider_WL_PFPolar->value();
                double currentWL = para->startWavel + indexWL*para->stepWavel;
                for (int i=0; i<para->nTheta; i++)
                {
                    theta[i] = 180.0 * i /(para->nTheta-1);
                    theta[2*para->nTheta-2 -i] =  360.0 - (180.0 * i /(para->nTheta-1));
                    if (ui->radioButton_PhaseAverage->isChecked())
                    {
                        phaseFunction[i] = para->phaseFunctionAve[indexWL][i];
                        phaseFunction[2*para->nTheta-2 -i] = para->phaseFunctionAve[indexWL][i];
                    }
                    if (ui->radioButton_PhasePara->isChecked())
                    {
                        phaseFunction[i] = para->phaseFunctionPara[indexWL][i];
                        phaseFunction[2*para->nTheta-2 -i] = para->phaseFunctionPara[indexWL][i];
                    }
                    if (ui->radioButton_PhasePerp->isChecked())
                    {
                        phaseFunction[i] = para->phaseFunctionPerp[indexWL][i];
                        phaseFunction[2*para->nTheta-2 -i] = para->phaseFunctionPerp[indexWL][i];
                    }
                }
                out << "Angle(deg)\t Phase Function";
                if (ui->radioButton_PhaseAverage->isChecked())
                    out << " (Ave) ";
                if (ui->radioButton_PhasePara->isChecked())
                     out << " (Para) ";
                if (ui->radioButton_PhasePerp->isChecked())
                     out << " (Perp) ";
                out << "@ Wavelength of " << currentWL << " nm\n";
                for (int i=0; i<2*para->nTheta-1; i++)
                    out << theta[i] <<"\t" << phaseFunction[i] <<"\n";
                out << "\n";
            }

            if (S1)
            {                
                out << "S1:\n";
                QVector<double>  theta(para->nTheta);
                int indexWL = ui->slider_WL_S1S2->value();
                double currentWL = para->startWavel + indexWL*para->stepWavel;
                out << "Angle(deg)\t S1(real | imag)" << " @ Wavelength of " << currentWL << " nm\n";
                for (int i=0; i<para->nTheta; i++)
                {
                    theta[i] = 180.0 * i /(para->nTheta-1);
                    out << theta[i] <<"\t" << para->S1[indexWL][i].real() <<"\t"
                        << para->S1[indexWL][i].imag()<<"\n";
                }
                out << "\n";
            }

            if (S2)
            {            
                QVector<double>  theta(para->nTheta);
                out << "S2:\n";
                int indexWL = ui->slider_WL_S1S2->value();
                double currentWL = para->startWavel + indexWL*para->stepWavel;
                out << "Angle(deg)\t S2(real | imag)" << " @ Wavelength of " << currentWL << " nm\n";
                for (int i=0; i<para->nTheta; i++)
                {
                    theta[i] = 180.0 * i /(para->nTheta-1);
                    out << theta[i] <<"\t" << para->S2[indexWL][i].real() <<"\t"
                        << para->S2[indexWL][i].imag()<<"\n";
                }
                out << "\n";
            }
        }
    }
}


void OptionsDialog::EnableFitPara()
{
    ui->checkBox_Fitting->setDisabled(false);
}

void OptionsDialog::DisableFitPara()
{
    ui->checkBox_Fitting->setDisabled(true);
}


void OptionsDialog::on_checkBox_Csca_clicked(bool checked)
{
    cSca = checked;
    ApplyAllEnableDisable();
    CheckSelectAll();
}

void OptionsDialog::on_checkBox_Cext_clicked(bool checked)
{
    cExt = checked;
    ApplyAllEnableDisable();
    CheckSelectAll();
}

void OptionsDialog::on_checkBox_Cback_clicked(bool checked)
{
    cBack = checked;
    ApplyAllEnableDisable();
    CheckSelectAll();
}

void OptionsDialog::on_checkBox_ScatCoeff_clicked(bool checked)
{
    scatCoeff = checked;
    ApplyAllEnableDisable();
    CheckSelectAll();
}

void OptionsDialog::on_checkBox_ReducedScatCoeff_clicked(bool checked)
{
    redScatCoeff = checked;
    ApplyAllEnableDisable();
    CheckSelectAll();
}

void OptionsDialog::on_checkBox_Fitting_clicked(bool checked)
{
    fitPara = checked;
    ApplyAllEnableDisable();
    CheckSelectAll();
}

void OptionsDialog::on_checkBox_G_clicked(bool checked)
{
    g = checked;
    ApplyAllEnableDisable();
    CheckSelectAll();
}

void OptionsDialog::on_checkBox_PhaseFunc_clicked(bool checked)
{
    phaseFunc = checked;
    ApplyAllEnableDisable();
    CheckSelectAll();
}

void OptionsDialog::on_checkBox_S1_clicked(bool checked)
{
    S1 = checked;
    ApplyAllEnableDisable();
    CheckSelectAll();
}

void OptionsDialog::on_checkBox_S2_clicked(bool checked)
{
    S2 = checked;
    ApplyAllEnableDisable();
    CheckSelectAll();
}

void OptionsDialog::on_checkBox_Forward_clicked(bool checked)
{
    forward = checked;
    ApplyAllEnableDisable();
    CheckSelectAll();
}

void OptionsDialog::on_checkBox_Backward_clicked(bool checked)
{
    backward = checked;
    ApplyAllEnableDisable();
    CheckSelectAll();
}

void OptionsDialog::on_checkBox_SelectAll_clicked(bool checked)
{
    cSca = checked;
    scatCoeff = checked;
    redScatCoeff = checked;
    fitPara = checked;
    g = checked;
    phaseFunc = checked;
    S1 = checked;
    S2 = checked;
    forward = checked;
    backward = checked;
    cExt = checked;
    cBack = checked;

    ui->checkBox_Csca->setChecked(checked);
    ui->checkBox_ScatCoeff->setChecked(checked);
    ui->checkBox_ReducedScatCoeff->setChecked(checked);
    ui->checkBox_Fitting->setChecked(checked);
    ui->checkBox_G->setChecked(checked);
    ui->checkBox_PhaseFunc->setChecked(checked);
    ui->checkBox_S1->setChecked(checked);
    ui->checkBox_S2->setChecked(checked);
    ui->checkBox_Forward->setChecked(checked);
    ui->checkBox_Backward->setChecked(checked);
    ui->checkBox_Cext ->setChecked(checked);
    ui->checkBox_Cback ->setChecked(checked);

    ApplyAllEnableDisable();
}
