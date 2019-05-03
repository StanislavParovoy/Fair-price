void createScriptFile(long double u[J + 1][N + 1], long double x[J + 1], long double t[J + 1]) {//dbg
	std::ofstream out;

	out.open("script.m");
	out << "t = 0:1 : " << N-1;
	out << endl;
	out << "x = 0:1 : " << J-2;
	out << endl;
	out << "[T, X] = meshgrid(t, x)";
	out << endl;
	out << "u = [";
	for (int j = 1; j <= J-1 ; j++) {
		for (int n = 1; n <= N  ; n++) {
			out << u[j][n] << ' ';
		}
		out << endl;
	}
	out << ']';
	out << endl;
	out << "plot3(X, T, u)";
	out.close();
}

long double err[J + 1][N + 1];

long double val = abs(u[j][n] - U[j][n]);
err[j][n] = val;

createScriptFile(err, x, t);
