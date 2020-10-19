#ifndef DISPLAYDIALOG_H
#define DISPLAYDIALOG_H

#include <QDialog>
#include "ui_mainwindow.h"
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
    void DisplayData(Ui_MainWindow *, parameters *);
    Ui::DisplayDialog *ui;

private slots:
    void on_pushButton_Close_clicked();

};

#endif // DISPLAYDIALOG_H
