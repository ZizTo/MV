#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

const double A = -3;
const double B = 3;
const int N = 16; // 15 интервалов, 16 узлов

double func(double x) {
    return sin(x) * cos(x);
}

double funcpr(double x) {
    return cos(2 * x);
}

double funcprpr(double x) {
    return -2 * sin(2 * x);
}

//LAB3 from second sem
vector<double> metodProgonki(vector<double> a, vector<double> b, vector<double> c, vector<double> f) {
    vector<double> alp(N, 0);
    vector<double> bet(N, 0);

    vector<double> y(N);


    alp[0] = -b[0] / c[0];
    bet[0] = f[0] / c[0];
    for (int i = 1; i < N - 1; i++)
    {
        double delti = c[i] + a[i] * alp[i - 1];
        alp[i] = -b[i] / delti;
        bet[i] = (f[i] - a[i] * bet[i - 1]) / delti;
    }
    bet[N - 1] = (f[N - 1] - a[N - 1] * bet[N - 2]) / (c[N - 1] + a[N - 1] * alp[N - 2]);

    y[N - 1] = funcprpr(B);
    for (int i = N - 2; i >= 0; i--)
    {
        y[i] = alp[i] * y[i + 1] + bet[i];
    }

    return y;
}

double S(double x, double h, const vector<double>& M, const vector<double>& f, const vector<double>& xm) {
    int i = 0;
    while (i < N - 1 && x > xm[i + 1])
        i++;

    return M[i] * pow(xm[i + 1] - x, 3) / (6 * h) +
        M[i + 1] * pow(x - xm[i], 3) / (6 * h) +
        (f[i] - M[i] * h * h / 6) * (xm[i + 1] - x) / h +
        (f[i + 1] - M[i + 1] * h * h / 6) * (x - xm[i]) / h;
}

int main() {
    vector<double> f(N), xm(N);
    double h = (B - A) / (N - 1);

    for (int i = 0; i < N; ++i) {
        xm[i] = A + i * h;
        f[i] = func(xm[i]);
    }

    vector<double> a(N, 1), b(N, 1), c(N, 4), ff(N, 0);

    // Левое условие
    c[0] = 2;
    b[0] = 1;
    ff[0] = (f[1] - f[0]) / h - funcpr(A) + h * funcprpr(A) / 6;

    // Правое условие
    a[N - 1] = 0;
    c[N - 1] = 1;
    ff[N - 1] = funcprpr(B);

    for (int i = 1; i < N - 1; ++i) {
        ff[i] = 6 * ((f[i + 1] - 2 * f[i] + f[i - 1]) / (h * h));
    }

    vector<double> M = metodProgonki(a, b, c, ff);

    cout << fixed << setprecision(20);

    double maxC = 0;
    int imax = -1;
    ofstream out("uzl.txt");
    for (int i = 0; i < 100; i++)
    {
        double x = A + (B - A) / 100 * i;
        if (maxC < abs(S(x, h, M, f, xm) - func(x))) { 
            maxC = abs(S(x, h, M, f, xm) - func(x));
            imax = i;
        }
        out << x << " " << S(x, h, M, f, xm) << endl;
    }
    out.close();

    cout << "max: " << maxC << " - " << imax << endl;

    return 0;
}