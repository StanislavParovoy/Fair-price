#include "stdafx.h"
#include <cmath>
#include <iostream>
#include "windows.h"
#include <fstream>
//Леланд
using namespace std;

const long double M_PIl = 3.141592653589793238462643383279502884L;

const int N = 20;
const int M = 160;

const long double R = 1.0;
const long double σ = 0.2;
const long double eps = 0.0;
const long double r = 0.1;
const long double T = (σ * σ) / 2;

const long double delta_t = 0.01;
const long double kappa = 0.05;

void createScriptFile(long double u[2 * N + 1][M + 1], long double x[2 * N + 1], long double t[2 * N + 1]) {
	std::ofstream out;

	out.open("script.m");
	out << "t = 0:1 : " << M;
	out << endl;
	out << "x = 0:1 : " << 2 * N;
	out << endl;
	out << "[T, X] = meshgrid(t, x)";
	out << endl;
	out << "u = [";
	for (int j = 0; j < 2 * N + 1; j++) {
		for (int n = 0; n < M + 1; n++) {
			out << u[j][n] << ' ';
		}
		out << endl;
	}
	out << ']';
	out << endl;
	out << "plot3(X, T, u)";
	out.close();
}


long double func_u0(long double x) {
	return fmax(1 - expl(-x), 0);
}

long double func_g0(long double tau) {
	return 0;
}

long double func_g1(long double tau) {
	return 1 - expl(-((2 * r*tau / (σ*σ)) + R));
}

long double sign(long double Val) {
	if (Val == 0.0) {
		return 0;
	}
	if (Val > 0.0) {
		return 1;
	}
	else {
		return -1;
	}
}

int main()
{
	long double h = R / N;
	long double k = T / M;

	long double tau[M + 1];
	for (int n = 0; n <= M; n++) {
		tau[n] = n * k;
	}

	long double x[2 * N + 1];
	for (int i = 0; i <= 2 * N; i++) {
		x[i] = i * h - R;
	}

	long double u[2 * N + 1][M + 1];

	for (int i = 0; i < 2 * N + 1; i++) {
		for (int n = 0; n < M + 1; n++) {
			u[i][n] = -7;
		}
	}

	for (int i = 0; i <= 2 * N; i++) {
		u[i][0] = func_u0(x[i]);
	}

	for (int n = 1; n <= M; n++) {
		u[0][n] = func_g0(tau[n]);
		u[2 * N][n] = func_g1(tau[n]);
	}

	long double a = 0.0;
	long double b = 0.0;
	long double c = 1.0;
	long double a2 = 0.0;
	long double b2 = 0.0;
	long double c2 = 0.0;
	long double d = 0.0;

	long double alpha[2 * N + 1];
	alpha[1] = b / c;

	long double betta[2 * N + 1];
	long double s[2 * N + 1][M + 1];

	for (int n = 1; n <= M; n++) {
		d = u[0][n];
		betta[1] = d / c;
		for (int i = 1; i <= 2 * N - 1; i++) {
			long double D2U = (u[i + 1][n - 1] - 2 * u[i][n - 1] + u[i - 1][n - 1]) / (h*h);
			long double D1U = (u[i + 1][n - 1] - u[i - 1][n - 1]) / (2 * h);

			s[i][n - 1] = sqrt(2 / M_PIl)*(kappa / (σ*sqrt(delta_t))) * fabs(sign(D2U + D1U));

			a = (s[i][n - 1] + 1) / (h*h) - (s[i][n - 1] + 1 + 2 * r / (σ*σ)) / (2 * h);
			b = (s[i][n - 1] + 1) / (h*h) + (s[i][n - 1] + 1 + 2 * r / (σ*σ)) / (2 * h);
			c = 2 / k + (2 * s[i][n - 1] + 2) / (h * h);

			a2 = a;
			b2 = b;
			c2 = 2 / k - (2 * s[i][n - 1] + 2) / (h * h);

			d = c2 * u[i][n - 1] + b2 * u[i + 1][n - 1] + a2 * u[i - 1][n - 1];
			alpha[i + 1] = b / (c - a * alpha[i]);
			betta[i + 1] = (d + a * betta[i]) / (c - a * alpha[i]);
		}
		d = u[2 * N][n];
		for (int ii = 2 * N; ii >= 2; ii--) {
			u[ii - 1][n] = alpha[ii] * u[ii][n] + betta[ii];
		}
	}

	createScriptFile(u, x, tau);

	cout << endl << endl;
	system("pause");
	return 0;
}

