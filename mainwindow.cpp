#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "qmath.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->setWindowTitle("Fair price");

    ui->lineEdit_N->setText("2");
    ui->lineEdit_r->setText("0.1");
    ui->lineEdit_aORd->setText("-0.2");
    ui->lineEdit_S->setText("100");
    ui->lineEdit_K->setText("100");
    ui->lineEdit_bORu->setText("0.2");
}

MainWindow::~MainWindow()
{
    delete ui;
}

//-- max
double max(double a,double b){
    if (a>=b){
        return a;
    }
    else{
        return b;
    }
}
//-- end max

//-- CX
int CX(int n, int m)
{
    int k = n - m;
    if (m > k)
        m = k;
    if (!m)
        return 1;
    int akk = k = n - m + 1, r1, r2;
    k++;
    for (int i = 2; i <= m; i++, k++)
    {
        if (!(r1 = akk % i))
        {
            akk = akk / i * k;
            continue;
        }
        if (!(r2 = k % i))
        {
            akk = k / i * akk;
            continue;
        }
        int t1 = akk / i, t2 = k / i;
        akk = t1 * t2 * i + t1 * r2 + t2 *r1 + r1 * r2 / i;
    }
    return akk;
}

//-- end CX

//--k0
int getK0(double S, double bORu, double aORd, double N, double K){
    int k0 = 0;
    while (true){
        if (S*qPow(1+bORu, k0)*qPow(1+aORd,N - k0) > K){
            break;
        }
        k0++;
    }
    return k0;
}
//--end k0

//-- p
double getP(double r, double aORd, double bORu){
    return (r - aORd)/(bORu - aORd);
}
//-- end p

//-- C0(s)
double getC0(int k0, int N, double p, double bORu, double aORd, double S, double K, double r){
    double sum = 0;
    for (int l = k0; l <= N; l ++){
        double tmp = CX(N, l) * qPow(p,l) * qPow(1-p,N-l) * max((qPow(1 + bORu,l) * qPow(1 + aORd, N - l) * S - K), 0);
        sum = sum + tmp;
    }

    return (1/(qPow(1+r, N))) * sum;
}
//end C0(s)

//-- P0(s)
double getP0(int k0, int N, double p, double bORu, double aORd, double S, double K, double r){
    double sum = 0;
    for (int l = 0; l <= k0 - 1; l ++){
        double tmp = CX(N, l) * qPow(p,l) * qPow(1-p,N-l) * max(K - (qPow(1 + bORu,l) * qPow(1 + aORd, N - l) * S), 0);
        sum = sum + tmp;
    }

    return (1/(qPow(1+r, N))) * sum;
}
//end P0(s)

//-- F0(S)
double getF0(double S, double K, double r, double N){
    return S - (K/qPow(1+r, N));
}
//-- end F0(S)

void MainWindow::on_pushButton_clicked()
{
    N = ui->lineEdit_N->text().toInt();
    r = ui->lineEdit_r->text().toDouble();
    aORd = ui->lineEdit_aORd->text().toDouble();
    S = ui->lineEdit_S->text().toDouble();
    K = ui->lineEdit_K->text().toDouble();
    bORu = ui->lineEdit_bORu->text().toDouble();

    //*******************************************
    //*******************************************
    //---- call  V(S,N)=max(S - K, 0)
    //*******************************************
    //*******************************************
    int k0 = getK0(S, bORu,aORd,N, K);
    double p =  getP( r,  aORd,  bORu);
    double C0 = getC0( k0,  N,  p,  bORu,  aORd,  S,  K,  r);
    ui->label_answer_call->setText(QString::number(C0));

    //*******************************************
    //*******************************************
    //---- put  V(S,N)=max(K - S, 0)
    //*******************************************
    //*******************************************
    // old K0 and p
    double P0 = getP0( k0,  N,  p,  bORu,  aORd,  S,  K,  r);
    ui->label_answer_put->setText(QString::number(P0));

    //*******************************************
    //*******************************************
    //---- futures  V(S,N)= S - K
    //*******************************************
    //*******************************************
    double F0 =  getF0( S,  K, r,  N);
     ui->label_answer_futures->setText(QString::number(F0));
}
