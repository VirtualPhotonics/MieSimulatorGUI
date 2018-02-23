#ifndef OPTIONSDIALOG_H
#define OPTIONSDIALOG_H

#include <QDialog>
#include "ui_mainwindow.h"
#include "ui_displaydialog.h"
#include "ui_optionsdialog.h"
#include "parameters.h"

namespace Ui {
class OptionsDialog;
}

class OptionsDialog : public QDialog
{
    Q_OBJECT

public:
    explicit OptionsDialog(QWidget *parent = 0);
    ~OptionsDialog();
    bool scatCross;
    bool scatCoeff;
    bool redScatCoeff;
    bool g;
    bool phaseFunc;
    bool S1;
    bool S2;
    bool forward;
    bool backward;
    bool fitPara;
    void EnableFitPara();
    void DisableFitPara();
    void SaveData(Ui_MainWindow *ui, parameters *para);
    void CheckSelectAll();
    void ApplyAllEnableDisable();


private slots:
    void on_pushButton_Apply_clicked();
    void on_pushButton_Cancel_clicked();
    void on_checkBox_ScatCross_clicked(bool checked);
    void on_checkBox_ScatCoeff_clicked(bool checked);
    void on_checkBox_ReducedScatCoeff_clicked(bool checked);
    void on_checkBox_G_clicked(bool checked);
    void on_checkBox_PhaseFunc_clicked(bool checked);
    void on_checkBox_S1_clicked(bool checked);
    void on_checkBox_S2_clicked(bool checked);
    void on_checkBox_Forward_clicked(bool checked);
    void on_checkBox_Backward_clicked(bool checked);
    void on_checkBox_Fitting_clicked(bool checked);
    void on_checkBox_SelectAll_clicked(bool checked);

private:
    Ui::OptionsDialog *ui;
};

#endif // OPTIONSDIALOG_H
