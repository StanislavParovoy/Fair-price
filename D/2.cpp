#include "stdafx.h"
#include <cmath>
#include <iostream>
#include "windows.h"
#include <fstream>

using namespace std;

const long double VAL_EXP = exp(1);

const int J = 20;//x
const int N = 80;//t

const long double X = 0.1;
const long double T = 2.0;

const long double _alpha = 0.01925;
const long double beta = 0.55;
const long double σ = 0.39;
const long double  γ = sqrt(beta*beta + 2.0 * σ*σ);

long double func_U(long double t) {
	long double result = 2.0 * γ *  expl(((beta + γ)*t) / (2.0));
	return result;
}

long double func_V(long double t) {
	long double result = 2.0 * γ + (beta + γ)*(expl((γ*t)) - 1);
	return result;
}

long double func_W(long double t) {
	long double result = 2.0 * (expl((γ*t)) - 1);
	return result;
}

long double func_A(long double t) {
	long double result = powl(func_U(t) / func_V(t), (2.0*_alpha) / (σ* σ));
	return result;
}

long double func_B(long double t) {
	long double result = func_W(t) / func_V(t);
	return result;
}

long double func_v(long double x, long double t) {
	long double result = func_A(t)*expl((-func_B(t)*x));
	return result;
}

long double func_u0(long double x) {
	return func_v(x, 0);
}

long double func_g0(long double t) {
	return func_v(0, t);
}

long double func_g1(long double t) {
	return func_v(X, t);
}

long double func_Mu(long double x) {
	return ((σ*σ) / (2.0))*x;
}

long double func_C(long double x) {
	return beta * x - _alpha + ((σ*σ) / (2.0));
}

long double func_f(long double x) {
	return beta - x;
}

void createScriptFile(long double u[J + 1][N + 1], long double x[J + 1], long double t[J + 1]) {//dbg
	std::ofstream out;

	out.open("script.m");
	out << "t = 0:1 : " << N;
	out << endl;
	out << "x = 0:1 : " << J;
	out << endl;
	out << "[T, X] = meshgrid(t, x)";
	out << endl;
	out << "u = [";
	for (int j = 0; j < J + 1; j++) {
		for (int n = 0; n < N + 1; n++) {
			out << u[j][n] << ' ';
		}
		out << endl;
	}
	out << ']';
	out << endl;
	out << "plot3(X, T, u)";
	out.close();
}

int main()
{
	long double   h = X / J;//шаг по x
	long double tau = T / N;//шаг по t

	long double x[J + 1];//сетка по x
	for (int j = 0; j <= J; j++) {
		x[j] = j * h;
	}

	long double t[N + 1];//сетка по t
	for (int n = 0; n <= N; n++) {
		t[n] = n * tau;
	}

	long double b[J + 1];  //коэффициенты системы уравнений
	long double a[J + 1];  //коэффициенты системы уравнений
	long double c[J + 1];  //коэффициенты системы уравнений
	long double b2[J + 1]; //коэффициенты системы уравнений
	long double a2[J + 1]; //коэффициенты системы уравнений
	long double c2[J + 1]; //коэффициенты системы уравнений

	c[0] = 1;
	c[J] = 1;
	b[0] = 0;
	a[J] = 0;

	for (int j = 1; j <= J - 1; j++) {
		a[j] = func_Mu(x[j]) / (2 * h*h) - ((0.5*σ*σ - func_C(x[j])) / (4 * h));
		b[j] = func_Mu(x[j]) / (2 * h*h) + ((σ*σ / 2 - func_C(x[j])) / (4 * h));
		c[j] = 1 / tau + func_Mu(x[j]) / (h*h) - (func_f(x[j]) - beta) / 2;

		a2[j] = a[j];
		b2[j] = b[j];
		c2[j] = 1 / tau - func_Mu(x[j]) / (h*h) + (func_f(x[j]) - beta) / 2;
	}

	long double u[J + 1][N + 1];//численное решение
	long double U[J + 1][N + 1];//точное решение 

	for (int j = 0; j < J + 1; j++) {
		for (int n = 0; n < N + 1; n++) {
			u[j][n] = -7;
			U[j][n] = -7;
		}
	}

	for (int j = 0; j <= J; j++) {
		for (int n = 0; n <= N; n++) {
			U[j][n] = func_v(x[j], t[n]);
		}
	}

	for (int j = 0; j <= J; j++) {
		u[j][0] = func_u0(x[j]);
	}

	for (int n = 1; n <= N; n++) {
		u[0][n] = func_g0(t[n]);
		u[J][n] = func_g1(t[n]);
	}

	long double alpha[J+1];
	alpha[1] = b[0] / c[0];
	for (int j = 1; j <= J - 1; j++) {
		alpha[j + 1] = b[j] / (c[j] - a[j] * alpha[j]);
	}

	long double betta[J+1];
	long double d[J+1];
	for (int n = 0; n <= N - 1; n++) {
		d[0] = u[0][n + 1];
		d[J] = u[J][n + 1];
		betta[1] = d[0] / c[0];

		for (int j = 1; j <= J - 1; j++) {
			d[j] = c2[j] * u[j][n] + b2[j] * u[j + 1][n] + a2[j] * u[j - 1][n];
			betta[j + 1] = (d[j] + a[j] * betta[j]) / (c[j] - a[j] * alpha[j]);
		}

		for (int j = J; j >= 2.0; j--) {
			u[j - 1][n + 1] = alpha[j] * u[j][n + 1] + betta[j];
		}
	}

	long double errInf = 0;
	long double errE = 0;
	long double Dinf = 0;
	long double DE = 0;

	for (int j = 1; j <= J - 1; j++) {
		long double valD = abs(u[j][N] - U[j][N]);
		if (Dinf < valD) {
			Dinf = valD;
		}

		DE = DE + valD * valD;
		for (int n = 1; n <= N + 1; n++) {
			long double val = abs(u[j][n] - U[j][n]);

			if (errInf < val) {
				errInf = val;
			}
			errE = errE + val * val;
		}
	}

	cout << endl << "D_e " << sqrt(h*DE);
	cout << endl << "D_inf " << Dinf;
	cout << endl << "Err_r " << sqrt(h*tau*errE);
	cout << endl << "Err_inf " << errInf;

	createScriptFile(u, x, t);

	cout << endl << endl;
	system("pause");

	return 0;
}
