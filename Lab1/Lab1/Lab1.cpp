#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

const double epsDih = 1e-2;
const double eps = 1e-7;

double func(double x) {
    return cos(x) + x / 4 - 0.5;
}

double dfunc(double x) {
    return -sin(x) + 0.25;
}

int main()
{
    cout << fixed << setprecision(10);

    double a, b;
    double dih;
    a = -5; b = 0;
    int k = 0;

    cout << "| k |      a        |       b       |       fa      |       fb      |      b-a      |" << endl;
    cout << "| " << k << " | " << a << " | " << b << " | " << func(a) << " | " << func(b) << " | " << b - a << " |" << endl;
    while (b - a > epsDih) {
        dih = (a + b) / 2;
        if (func(b) * func(dih) < 0) {
            a = dih;
        }
        else {
            b = dih;
        }
        k++;
        cout << "| " << k << " | " << a << " | " << b << " | " << func(a) << " | " << func(b) << " | " << b - a << " |" << endl;
    }
    cout << "Dihotomia: " << dih << endl << "a: " << a << endl << "b: " << b << endl << endl;

    cout << "|         Nuton with const         |" << endl;
    cout << "| k |      x        |  |x-x(-1)|   |" << endl;
    double xNC = dih;
    double dfx = dfunc(dih);
    double fx = func(dih);
    double xk1 = 1000;
    k = 0;
    cout << "| " << k << " | " << xNC << " |       -      |" << endl;
    while (abs(xNC - xk1) > eps)
    {
        xk1 = xNC;
        xNC -= fx / dfx;
        fx = func(xNC);
        k++;
        cout << "| " << k << " | " << xNC << " | " << abs(xNC - xk1) << " |" << endl;
    }

    cout << xNC << " - Nuton with const x" << endl << endl;

    cout << "|       Nuton without const        |" << endl;
    cout << "| k |      x        |  |x-x(-1)|   |" << endl;
    double xN = dih;
    fx = func(dih);
    k = 0;
    xk1 = 1000;
    cout << "| " << k << " | " << xN << " |       -      |" << endl;
    while (abs(xN - xk1) > eps)
    {
        xk1 = xN;
        xN -= fx / dfunc(dih);
        fx = func(xN);
        k++;
        cout << "| " << k << " | " << xN << " | " << abs(xN - xk1) << " |" << endl;
    }

    cout << xN << " - Nuton without const x" << endl << endl;

    cout << "|             Sekushih             |" << endl;
    cout << "| k |      x        |  |x-x(-1)|   |" << endl;
    double fxa = func(a), fxb = func(b);
    k = 0;
    cout << "| " << k << " | " << b << " |       -      |" << endl;
    while (abs(b - a) > eps) {
        double c = b - fxb * (b - a) / (fxb - fxa);
        a = b; b = c; fxa = fxb; fxb = func(c);
        k++;
        cout << "| " << k << " | " << b << " | " << abs(b - a) << " |" << endl;
    }

    cout << b << " - Metod sekushih" << endl;
}