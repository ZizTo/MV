#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

double f(double x, double y) {
    return x * x * (y * y + 1);
}

double df(double x, double y) {
    return 2 * x * x * y;
}

double exact(double x) {
    return tan(x*x*x/3);
}

double trapezoidal_step(double x, double y, double h, double tol = 1e-10) {
    double y_next = y + h * f(x, y);
    double error = 1.0;
    int max_iter = 100;

    for (int iter = 0; iter < max_iter && error > tol; ++iter) {
        double F = y_next - y - h / 2 * (f(x, y) + f(x + h, y_next));
        double dF = 1 - h / 2 * df(x + h, y_next);
        double delta = -F / dF;
        y_next += delta;
        error = fabs(delta);
    }
    return y_next;
}

int main()
{
    vector<vector<double>> rt{ 
        {0,   0,   0,   0,   0},
        {1.0/2, 1.0/2, 0,   0,   0},
        {1.0/2, 0,   1.0/2, 0,   0},
        {1.0,   0,   0,   1.0,   0},
        {0,   1.0/6, 1.0/3, 1.0/3, 1.0/6}
    };

    double a = 0, b = 1;
    double h = 0.1;
    vector<double> k(4);
    vector<double> yh((b-a)/h+1);
    vector<double> yh2((b-a)/(h/2) + 1);
    yh[0] = 0;
    yh2[0] = 0;

    for (double i = 0; i < (b - a) / h; i++)
    {
        k[0] = f(i * h, yh[i]);
        k[1] = f(i * h + rt[1][0]*h, yh[i] + h*rt[1][1]*k[0]);
        k[2] = f(i * h + rt[2][0]*h, yh[i] + h*rt[2][1]*k[0] + h * rt[2][2] * k[1]);
        k[3] = f(i * h + rt[3][0]*h, yh[i] + h*rt[3][1]*k[0] + h * rt[3][2] * k[1] + h * rt[3][3] * k[2]);

        yh[i + 1] = yh[i] + h * (rt[4][1]*k[0] + rt[4][2] * k[1] + rt[4][3] * k[2] + rt[4][4] * k[3]);
    }

    h /= 2;

    double maxnet = 0, maxnety = 0;

    ofstream out("xdhd.txt");
    for (int i = 0; i < (b - a) / h; i++)
    {
        k[0] = f(i * h, yh2[i]);
        k[1] = f(i * h + rt[1][0] * h, yh2[i] + h * rt[1][1] * k[0]);
        k[2] = f(i * h + rt[2][0] * h, yh2[i] + h * rt[2][1] * k[0] + h * rt[2][2] * k[1]);
        k[3] = f(i * h + rt[3][0] * h, yh2[i] + h * rt[3][1] * k[0] + h * rt[3][2] * k[1] + h * rt[3][3] * k[2]);

        yh2[i + 1] = yh2[i] + h * (rt[4][1] * k[0] + rt[4][2] * k[1] + rt[4][3] * k[2] + rt[4][4] * k[3]);
        out << i * h << " " << yh2[i] << endl;
        if (maxnet < fabs(exact(i * h) - yh2[i])) { maxnet = fabs(exact(i * h) - yh2[i]); }
        if (i % 2 == 0 && maxnety < fabs(yh[i / 2] - yh2[i])) { maxnety = fabs(yh[i / 2] - yh2[i]); }
    }
    out.close();

    cout << fixed << setprecision(20);
    cout << "with exact: " << maxnet << endl;
    cout << "   with h1: " << maxnety << endl;
    cout << "       Rh2: " << maxnety / (16 - 1) << endl;



    h = 0.1;
    vector<double> ytr((b - a) / h + 1);
    vector<double> ytr2((b - a) / (h/2) + 1);
    ytr[0] = 0;
    ytr2[0] = 0;

    for (int i = 0; i < (b - a) / h; i++) {
        double x = a + i * h;
        ytr[i + 1] = trapezoidal_step(x, ytr[i], h);
    }

    h /= 2;
    
    ofstream out2("xdhd2.txt");
    double maxnettr = 0, maxnettry = 0;
    for (int i = 0; i < (b - a) / h; i++) {
        ytr2[i + 1] = trapezoidal_step(i * h, ytr2[i], h);

        out2 << i * h << " " << ytr2[i] << endl;
        if (maxnettr < fabs(exact(i*h) - ytr2[i])) { maxnettr = fabs(exact(i * h) - ytr2[i]); }
        if (i % 2 == 0 && maxnety < fabs(ytr[i / 2] - ytr2[i])) { maxnettry = fabs(ytr[i / 2] - ytr2[i]); }
    }
    out2.close();


    cout << "------------------METOD TRAPECI-----------------" << endl;
    cout << "with exact: " << maxnettr << endl;
    cout << "   with h1: " << maxnettry << endl;
    cout << "       Rh2: " << maxnettry / (16 - 1) << endl;
}
