#ifndef DISPLAYDIALOG_H
#define DISPLAYDIALOG_H

#include <QDialog>
#include <QRadioButton>
#include <QComboBox>
#include <QSlider>
#include "parameters.h"

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
                     QRadioButton *radioButton_PhaseAverage,
                     QRadioButton *radioButton_PhasePara,
                     QRadioButton *radioButton_PhasePerp, parameters *);
    Ui::DisplayDialog *ui;

private slots:
    void on_pushButton_Close_clicked();

};

#endif // DISPLAYDIALOG_H
