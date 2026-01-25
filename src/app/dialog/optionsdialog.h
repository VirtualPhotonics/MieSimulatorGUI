#ifndef OPTIONSDIALOG_H
#define OPTIONSDIALOG_H

#include <QDialog>
#include <QRadioButton>
#include <QComboBox>
#include <QSlider>
#include "qcheckbox.h"
#include "ui_optionsdialog.h"
#include "parameters.h"

namespace Ui {
class OptionsDialog;
}

class OptionsDialog : public QDialog
{
    Q_OBJECT

public:
    explicit OptionsDialog(QWidget *parent = nullptr);
    ~OptionsDialog();

    void SaveData(QRadioButton *radioButton_MonoDisperse,
                  QRadioButton *radioButton_PolyDisperse,
                  QRadioButton *radioButton_NumDen,
                  QRadioButton *radioButton_VolFrac,
                  QComboBox *comboBox_Distribution,
                  QSlider *slider_ConcPercentChange,
                  QCheckBox *checkBox_PhasePolarAve,
                  QCheckBox *checkBox_PhasePolarPara,
                  QCheckBox *checkBox_PhasePolarPerp,
                  Parameters *para);

private slots:
    void on_pushButton_ScatPara_clicked();
    void on_pushButton_PhaseFunction_clicked();
    void on_pushButton_S1_clicked();
    void on_pushButton_S2_clicked();
    void on_pushButton_Cancel_clicked();

private:
    bool mFlagScatPara;
    bool mFlagPhaseFunction;
    bool mFlagS1;
    bool mFlagS2;
    Ui::OptionsDialog *ui;

    void SaveScatPara(QRadioButton *radioButton_MonoDisperse,
                      QRadioButton *radioButton_PolyDisperse,
                      QRadioButton *radioButton_NumDen,
                      QRadioButton *radioButton_VolFrac,
                      QComboBox *comboBox_Distribution,
                      Parameters *para, QString fileName,
                      double margin);
    void SavePhaseFunction(QCheckBox *checkBox_PhasePolarAve,
                           QCheckBox *checkBox_PhasePolarPara,
                           QCheckBox *checkBox_PhasePolarPerp,
                           Parameters *para, QString fileName);
    void SaveS1(Parameters *para, QString fileName);
    void SaveS2(Parameters *para, QString fileName);
    void RememberLastDirectory(QString fileName);
};

#endif // OPTIONSDIALOG_H
