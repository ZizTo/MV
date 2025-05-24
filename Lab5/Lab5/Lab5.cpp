#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double f(double x, double y) {
    return x * x * (y * y + 1);
}

double exact(double x) {
    return tan(x*x*x/3);
}

int main()
{
    vector<vector<double>> rt{ 
        {0,   0,   0,   0,   0},
        {1/2, 1/2, 0,   0,   0},
        {1/2, 0,   1/2, 0,   0},
        {1,   0,   0,   1,   0},
        {0,   1/6, 1/3, 1/3, 1/6}
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
    for (double i = 0; i < (b - a) / h; i++)
    {
        k[0] = f(i * h, yh2[i]);
        k[1] = f(i * h + rt[1][0] * h, yh2[i] + h * rt[1][1] * k[0]);
        k[2] = f(i * h + rt[2][0] * h, yh2[i] + h * rt[2][1] * k[0] + h * rt[2][2] * k[1]);
        k[3] = f(i * h + rt[3][0] * h, yh2[i] + h * rt[3][1] * k[0] + h * rt[3][2] * k[1] + h * rt[3][3] * k[2]);

        yh2[i + 1] = yh2[i] + h * (rt[4][1] * k[0] + rt[4][2] * k[1] + rt[4][3] * k[2] + rt[4][4] * k[3]);
        if (maxnet < abs(exact(i * h) - yh2[i])) { maxnet = abs(exact(i * h) - yh2[i]); }
        if (int(i) % 2 == 0 && maxnety < abs(yh[i / 2] - yh2[i])) { maxnety = abs(yh[i / 2] - yh2[i]); }
    }


}
