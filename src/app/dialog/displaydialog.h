#ifndef DISPLAYDIALOG_H
#define DISPLAYDIALOG_H

#include <QDialog>
#include <QRadioButton>
#include <QComboBox>
#include <QSlider>
#include "parameters.h"
#include "qcheckbox.h"

namespace Ui {
class DisplayDialog;
}

class DisplayDialog : public QDialog
{
    Q_OBJECT

public:
    explicit DisplayDialog(QWidget *parent = 0);
    ~DisplayDialog();
    void DisplayData(QRadioButton *radioButton_MonoDisperse,
                     QRadioButton *radioButton_PolyDisperse,
                     QRadioButton *radioButton_NumDen,
                     QRadioButton *radioButton_VolFrac,
                     QComboBox *comboBox_Distribution,
                     QSlider *slider_ConcPercentChange,
                     QSlider *slider_WL_PFPolar,
                     QCheckBox *checkBox_PhasePolarAve,
                     QCheckBox *checkBox_PhasePolarPara,
                     QCheckBox *checkBox_PhasePolarPerp, Parameters *);
    Ui::DisplayDialog *ui;

private slots:
    void on_pushButton_Close_clicked();

};

#endif // DISPLAYDIALOG_H
