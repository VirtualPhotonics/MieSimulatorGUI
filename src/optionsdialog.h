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
    void on_pushButton_scatPara_clicked();
    void on_pushButton_phaseFunction_clicked();
    void on_pushButton_s1_clicked();
    void on_pushButton_s2_clicked();
    void on_pushButton_cancel_clicked();

private:
    bool _flagScatPara;
    bool _flagPhaseFunction;
    bool _flagS1;
    bool _flagS2;
    Ui::OptionsDialog *ui;

    void SaveScatPara(Ui_MainWindow *ui, parameters *para, QString fileName, double margin);
    void SavePhaseFunction(Ui_MainWindow *ui, parameters *para, QString fileName);
    void SaveS1(parameters *para, QString fileName);
    void SaveS2(parameters *para, QString fileName);
    void RememberLastDirectory(QString fileName);
};

#endif // OPTIONSDIALOG_H
