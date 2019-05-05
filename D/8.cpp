#include "stdafx.h"
#include <cmath>
#include <iostream>
#include "windows.h"
#include <fstream>

using namespace std;

const int N = 20;
const int M = 160;

const long double R = 1.0;
const long double σ = 0.2;
const long double eps = 0.0;
const long double r = 0.1;
const long double T = σ * σ / 2;
const long double σ_tsq = σ * σ * (1 + eps);

void createScriptFile(long double u[2 * N + 1][M + 1], long double x[2 * N + 1], long double t[2 * N + 1]) {//dbg
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
	return 1 - expl(-(2 * r*tau / (σ*σ) + R));
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

	long double s[2 * N + 1][M + 1];
	for (int i = 0; i < 2 * N + 1; i++) {
		for (int n = 0; n < M + 1; n++) {
			s[i][n] = σ_tsq / (σ*σ) - 1;
		}
	}

	long double a[2 * N + 1];
	long double b[2 * N + 1];
	long double c[2 * N + 1];
	long double a2[2 * N + 1];
	long double b2[2 * N + 1];
	long double c2[2 * N + 1];

	c[0] = 1;
	c[2 * N] = 1;
	b[0] = 0;
	a[2 * N] = 0;

	for (int i = 1; i <= 2 * N - 1; i++) {
		a[i] = (s[i][0] + 1) / (h*h) - (s[i][0] + 1 + 2 * r / (σ*σ)) / (2 * h);
		b[i] = (s[i][0] + 1) / (h*h) + (s[i][0] + 1 + 2 * r / (σ*σ)) / (2 * h);
		c[i] = 2 / k + (2 * s[i][0] + 2) / (h * h);

		a2[i] = a[i];
		b2[i] = b[i];
		c2[i] = 2 / k - (2 * s[i][0] + 2) / (h * h);
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

	long double alpha[2 * N + 1];
	alpha[1] = b[0] / c[0];
	for (int i = 1; i <= 2 * N - 1; i++) {
		alpha[i + 1] = b[i] / (c[i] - a[i] * alpha[i]);
	}

	long double betta[2 * N + 1];
	long double d[2 * N + 1];
	for (int n = 0; n <= M - 1; n++) {
		d[0] = u[0][n + 1];
		d[2 * N] = u[2 * N][n + 1];
		betta[1] = d[0] / c[0];

		for (int i = 1; i <= 2 * N - 1; i++) {
			d[i] = c2[i] * u[i][n] + b2[i] * u[i + 1][n] + a2[i] * u[i - 1][n];
			betta[i + 1] = (d[i] + a[i] * betta[i]) / (c[i] - a[i] * alpha[i]);
		}

		for (int i = 2 * N; i >= 2; i--) {
			u[i - 1][n + 1] = alpha[i] * u[i][n + 1] + betta[i];
		}
	}

	createScriptFile(u, x, tau);

	cout << endl << endl;
	system("pause");
	return 0;
}

