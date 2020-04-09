#include <iostream>
#include <fstream>
#include <cmath>

void coefC(double* cArr, double** mass, int n, double h) {
	int i, j;
	double** A = new double* [n - 2];
	for (i = 0; i < n - 2; i++) {
		A[i] = new double[n - 2];
	}
	for (i = 1; i < n - 3; i++) {
		j = i - 1;
		A[i][j] = h;
		A[i][j + 1] = 4 * h;
		A[i][j + 2] = h;
	}
	A[0][0] = 4 * h;
	A[0][1] = h;
	A[n - 3][n - 4] = h;
	A[n - 3][n - 3] = 4 * h;

	double** F = new double* [n - 2];
	for (i = 0; i < n - 2; i++) {
		F[i] = new double[1];
	}

	for (i = 1; i < n - 1; i++) {
		F[i - 1][0] = 3 * ((mass[i + 1][1] - mass[i][1]) / h - (mass[i][1] - mass[i - 1][1]) / h);
	}
	for (i = 1; i < n - 2; i++) {
		F[i][0] -= (A[i][i - 1] * F[i - 1][0] / A[i - 1][i - 1]);
		A[i][i] -= (A[i][i - 1] * A[i - 1][i] / A[i - 1][i - 1]);
	}

	for (i = 1; i < n - 2; i++) {
		F[i][0] -= (A[i][i - 1] * F[i - 1][0] / A[i - 1][i - 1]);
		A[i][i] -= (A[i][i - 1] * A[i - 1][i] / A[i - 1][i - 1]);
	}

	cArr[n - 2] = F[n - 3][0] / A[n - 3][n - 3];
	for (i = n - 4; i >= 0; i--) {
		cArr[i + 1] = (F[i][0] - A[i][i + 1] * cArr[i + 2]) / A[i][i];
	}
	cArr[0] = 0;
	cArr[n - 1] = 0;
	for (int i = 0; i < n - 2; i++) {
		delete[] A[i];
		delete[] F[i];
	}
	delete[] A;
	delete[] F;
}

double splineFunc(double x, double* a, double* b, double* c, double* d, double** func, int n) {
	int i = 0;
	for (int j = 1; j < n - 1; j++) {
		if (x > func[j][0])
			i++;
		else
			break;
	}
	return a[i] + b[i] * (x - func[i][0]) + c[i] * pow((x - func[i][0]), 2) + d[i] * pow((x - func[i][0]), 3);
}

double trapeze(double x0, double xn, int step, double** func, int n, double* a, double* b, double* c, double* d) {
	double sum = 0, delta = (xn - x0) / (double)step;
	while(x0 <= xn) {
		sum += (splineFunc(x0, a, b, c, d, func, n) + splineFunc(x0 + delta, a, b, c, d, func, n)) * delta / 2.;
		x0 += delta;
	}
	return sum;
}

double simpson(double x0, double xn, int step, double** func, int n, double* a, double* b, double* c, double* d) {
	double sum = 0, delta = (xn - x0) / (double)step;
	while (x0 <= xn) {
		sum += delta * (splineFunc(x0, a, b, c, d, func, n) + 4 * splineFunc(x0 + delta / 2., a, b, c, d, func, n) + splineFunc(x0 + delta, a, b, c, d, func, n)) / 6.;
		x0 += delta;
	}
	return sum;
}

double rectangle(double x0, double xn, int step, double** func, int n, double* a, double* b, double* c, double* d) {
	double sum = 0, delta = (xn - x0) / (double)step;
	while (x0 <= xn) {
		sum += splineFunc(x0 + delta / 2., a, b, c, d, func, n) * delta;
		x0 += delta;
	}
	return sum;
}

double runge(const int meth, double x0, double xn, int step, double** func, int n, double* a, double* b, double* c, double* d) {
	if (meth == 1)
		return fabs(trapeze(x0, xn, step, func, n, a, b, c, d) - trapeze(x0, xn, 2 * step, func, n, a, b, c, d)) / 3.;
	else if (meth == 2)
		return fabs(simpson(x0, xn, step, func, n, a, b, c, d) - simpson(x0, xn, 2 * step, func, n, a, b, c, d)) / 15.;
	else if (meth == 3)
		return fabs(rectangle(x0, xn, step, func, n, a, b, c, d) - rectangle(x0, xn, 2 * step, func, n, a, b, c, d)) / 3.;
}

int main(void) {
	using namespace std;
	int i, n = 9, step;
	double h = 0.25, eps = 1e-4;
	double** func = new double* [n];
	for (i = 0; i < n; i++) {
		func[i] = new double[2];
	}

	double j = 0;
	for (i = 0; i < n; i++) {
		func[i][0] = j;
		j += h;
	}
	func[0][1] = 1.;
	func[1][1] = 0.979915;
	func[2][1] = 0.927295;
	func[3][1] = 0.858001;
	func[4][1] = 0.785398;
	for (i = 5; i < n; i++) {
		func[i][1] = 0.716844;
	}

	double* a = new double[n - 1];
	double* b = new double[n - 1];
	double* c = new double[n];
	double* d = new double[n - 1];

	coefC(c, func, n, h);
	for (i = 0; i < n - 1; i++) {
		a[i] = func[i][1];
		b[i] = (func[i + 1][1] - func[i][1]) / h - c[i] * h - (c[i + 1] - c[i]) * h / 3.;
		d[i] = (c[i + 1] - c[i]) / (3 * h);
	}

	fstream fin;
	fin.open("coef.txt", ios::out);
	if (fin.is_open()) {
		for (i = 0; i < n - 1; i++) {
			fin << a[i] << " " << b[i] << " " << c[i] << " " << d[i] << endl;
		}
	}
	fin.close();
	fin.open("x.txt", ios::out);
	if (fin.is_open()) {
		for (i = 0; i < n; i++) {
			fin << func[i][0] << " " << func[i][1] << endl;
		}
	}
	fin.close();

	double integ, pogr;

	step = 50;
	do {
		integ = trapeze(func[0][0], func[n - 1][0], step, func, n, a, b, c, d);
		pogr = runge(1, func[0][0], func[n - 1][0], step, func, n, a, b, c, d);
		step *= 2;
	} while (pogr > eps);
	fin.open("ans1.dat", ios::out);
	if (fin.is_open()) {
		fin << integ << " " << pogr << endl;
	}
	fin.close();

	step = 50;
	do {
		integ = simpson(func[0][0], func[n - 1][0], step, func, n, a, b, c, d);
		pogr = runge(2, func[0][0], func[n - 1][0], step, func, n, a, b, c, d);
		step *= 2;
	} while (pogr > eps);
	fin.open("ans2.dat", ios::out);
	if (fin.is_open()) {
		fin << integ << " " << pogr << endl;
	}
	fin.close();

	step = 50;
	do {
		integ = rectangle(func[0][0], func[n - 1][0], step, func, n, a, b, c, d);
		pogr = runge(3, func[0][0], func[n - 1][0], step, func, n, a, b, c, d);
		step *= 2;
	} while (pogr > eps);
	fin.open("ans3.dat", ios::out);
	if (fin.is_open()) {
		fin << integ << " " << pogr << endl;
	}
	fin.close();

	system("python3 splineGraph.py");
	delete[] a;
	delete[] b;
	delete[] c;
	delete[] d;
	for (int i = 0; i < n; i++)
		delete[] func[i];
	delete[] func;
	return 0;
}