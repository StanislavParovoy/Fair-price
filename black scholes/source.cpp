//  C=SN(d_{1})-Xe^{-rT}N(d_{2}),}
//https://www.kalkulaator.ee/ru/opcionnyj-kalkulyator-model-blehka-shoulza

#include <iostream>
#include <cmath>

using namespace std;

double pnorm(double x, double mean, double sd)
{
	// Значение функции распределения нормального распределения 
	// со средним mean и стандартным отклонением sd

	if (sd < 0)
		return 0.0; // TODO: Ошибка или NAN

					// Нормируем x:
	if (sd == 0)
		sd = 1E-10;
	/* TODO: Рассмотреть случай sd == 0. R выдает что-то разумное */
	x = (x - mean) / sd;

	// Взято из http://www.johndcook.com/cpp_phi.html
	// constants
	double a1 = 0.254829592;
	double a2 = -0.284496736;
	double a3 = 1.421413741;
	double a4 = -1.453152027;
	double a5 = 1.061405429;
	double p = 0.3275911;

	// Save the sign of x
	int sign = 1;
	if (x < 0)
		sign = -1;
	x = fabs(x) / sqrt(2.0);

	// A&S formula 7.1.26 [A&S = Abramowitz & Stegun?]
	double t = 1.0 / (1.0 + p * x);
	double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x * x);

	return 0.5*(1.0 + sign * y);
}

double getd1(double S, double sigma, double T, double r, double X) {
	return (log(S / X) + (r + ((sigma * sigma) / (2)))*T) / (sigma*sqrt(T));
}

double getd2(double d1, double sigma, double T) {
	return d1 - sigma * sqrt(T);
}


double calcCallWiki(double S, double sigma, double T, double r, double X) {
	double d1 = getd1(S, sigma, T, r, X);
	double d2 = getd2(d1, sigma, T);
	double Notd1 = pnorm(d1, 0, 1);
	cout << endl << "N(d1)= " << Notd1;
	double Notd2 = pnorm(d2, 0, 1);
	cout << endl << "N(d2)= " << Notd2;
	return S * Notd1 - X * exp(-r * T)*Notd2;
}

double calcPutWiki(double S, double sigma, double T, double r, double X) {
	double d1 = getd1(S, sigma, T, r, X);
	cout << endl << "d1= " << d1;
	double d2 = getd2(d1, sigma, T);
	cout << endl << "d2= " << d2;
	double Notmind1 = pnorm(-d1, 0, 1);
	cout << endl << "N(-d1)= " << Notmind1;
	double Notmind2 = pnorm(-d2, 0, 1);
	cout << endl << "N(-d2)= " << Notmind2;
	return X * exp(-r * T)*Notmind2 - S * Notmind1;
}

int main()
{
	setlocale(LC_ALL, "Russian");

	double S = 25.0;//цена БА
	double X = 30.0;//страйк
	double sigma = 0.15;//% волатильность
	double TT = 60;//период
	double r = 0.00;//% Безрисковая ставка доходности
	double days = 365;//дней в году 

	cout << endl << "S - Текущая цена базисной акции. S = " << S;
	cout << endl << "X - Цена исполнения опциона (strike price). X = " << X;
	cout << endl << "sigma - Волатильность базисной акции %. sigma = " << sigma;
	cout << endl << "TT - Время до истечения срока опциона (период опциона). TT = " << TT;
	cout << endl << "r - Безрисковая ставка доходности (непрерывное накоплениеx). r = " << r;
	cout << endl << "days - Дней в году = " << 365;


	cout << endl << "===========================";
	double T = TT / days;
	cout << endl << "-> T - Срок опциона (лет)= " << T;
	cout << endl << "-> C — цена опциона call. C = " << calcCallWiki(S, sigma, T, r, X);
	cout << endl << "-> P — цена опциона put . P = " << calcPutWiki(S, sigma, T, r, X);

	cin >> X;
	return 0;
}

