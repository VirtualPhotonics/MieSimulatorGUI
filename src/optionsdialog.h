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
    explicit OptionsDialog(QWidget *parent = nullptr);
    ~OptionsDialog();

    void SaveData(Ui_MainWindow *ui, parameters *para);

private slots:
    void on_pushButton_ScatPara_clicked();
    void on_pushButton_PhaseFunction_clicked();
    void on_pushButton_S1_clicked();
    void on_pushButton_S2_clicked();
    void on_pushButton_Cancel_clicked();

private:
    bool flagScatPara;
    bool flagPhaseFunction;
    bool flagS1;
    bool flagS2;
    Ui::OptionsDialog *ui;

    void SaveScatPara(Ui_MainWindow *ui, parameters *para, QString fileName, double margin);
    void SavePhaseFunction(Ui_MainWindow *ui, parameters *para, QString fileName);
    void SaveS1(parameters *para, QString fileName);
    void SaveS2(parameters *para, QString fileName);
    void RememberLastDirectory(QString fileName);
};

#endif // OPTIONSDIALOG_H
