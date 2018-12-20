#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

    double r =0;
    double K =0;
    double S =0;
    int N =0;
    double aORd =0;
    double bORu =0;

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();


private:
    Ui::MainWindow *ui;


private slots:
    void on_pushButton_clicked();
};

#endif // MAINWINDOW_H
