#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;


const double exactValue = M_PI / (3.0 * sqrt(3.0));
const double eps = 1e-7;

double f(double x) {
    return 1.0 / (1.0 + 2.0 * pow(sin(x), 2));
}

double srPriamoug(double a, double b, int n) {
    double h = (b - a) / n;
    double result = 0.0;
    
    for (int i = 0; i < n; i++) {
        double x_mid = a + (i + 0.5) * h;
        result += f(x_mid);
    }
    
    return h * result;
}

double simpsonKF(double a, double b, int n) {
    double h = (b - a) / (2*n);
    double result = f(a) + f(b);
    
    for (int i = 1; i < 2*n; i += 2) {
        result += 4 * f(a + i * h);
    }
    
    for (int i = 2; i < 2*n; i += 2) {
        result += 2* f(a + i * h);
    }
    
    return h * result / 3.0;
}

double rungeError(double Qh, double Qh2, int p) {
    return abs(Qh - Qh2) / (pow(2, p) - 1.0);
}

void gaussKF4(vector<double>& t, vector<double>& A) {
    t.resize(4);
    t[0] = -sqrt((3.0 + 2.0 * sqrt(6.0/5.0)) / 7.0);
    t[1] = -sqrt((3.0 - 2.0 * sqrt(6.0/5.0)) / 7.0);
    t[2] = sqrt((3.0 - 2.0 * sqrt(6.0/5.0)) / 7.0);
    t[3] = sqrt((3.0 + 2.0 * sqrt(6.0/5.0)) / 7.0);
    
    A.resize(4);
    A[0] = (18.0 - sqrt(30.0)) / 36.0;
    A[1] = (18.0 + sqrt(30.0)) / 36.0;
    A[2] = (18.0 + sqrt(30.0)) / 36.0;
    A[3] = (18.0 - sqrt(30.0)) / 36.0;
}

double gaussKF(double a, double b, int k) {
    vector<double> t, A;
    
    switch(k) {
        case 4:
        gaussKF4(t, A);
            break;
        default:
            cout << "Unsupported number of nodes: " << k << endl;
            return 0;
    }
    
    double result = 0;

    // Преобразование отрезка [a, b] к [-1, 1]
    double mid = (b + a) / 2.0;
    double half_length = (b - a) / 2.0;
    
    for (int i = 0; i < k; i++) {
        double x = mid + half_length * t[i];
        result += A[i] * f(x);
    }
    
    return half_length * result;
}

int main() {
    double a = 0.0;
    double b = M_PI / 4.0;
    
    cout << fixed << setprecision(10);
    
    cout << "Задание 1" << endl;
    
    cout << "\nКвадратурная формула средних прямоугольников:" << endl;
    cout << "----------------------------------------------------------------------------------" << endl;
    cout << "| Число разбиений|      Шаг      |  Приближенное |   Оценка      |  Абсолютная   |" << endl;
    cout << "|       N        |       h       |   значение    |  погрешности  |   погрешность |" << endl;
    cout << "----------------------------------------------------------------------------------" << endl;

    int por1 = 2;
    int n = 2;
    double err = 1.0;
    
    while (abs(err) > eps) {
        double h = (b - a) / n;
        double Qh = srPriamoug(a, b, n/2);
        double Qh2 = srPriamoug(a, b, n);
        double error = rungeError(Qh, Qh2, por1);
        double abs_error = abs(Qh2 - exactValue);
        
        cout << "| " << setw(14) << n;
        cout << " | " << setw(13) << h;
        cout << " | " << setw(13) << Qh2;
        cout << " | " << setw(13) << error;
        cout << " | " << setw(13) << abs_error << " |" << endl;
        
        err = error;
        n *= 2;
        
        if (n > 1000000) {
            break;
        }
    }
    cout << "----------------------------------------------------------------------------------" << endl;
    
    cout << "\nКвадратурная формула Симпсона:" << endl;
    cout << "----------------------------------------------------------------------------------" << endl;
    cout << "| Число разбиений|      Шаг      |  Приближенное |   Оценка      |  Абсолютная   |" << endl;
    cout << "|       N        |       h       |   значение    |  погрешности  |   погрешность |" << endl;
    cout << "----------------------------------------------------------------------------------" << endl;
    
    int por2 = 4;
    n = 2;
    err = 1.0;
    
    while (abs(err) > eps) {
        double h = (b - a) / n;
        double Qh = simpsonKF(a, b, n/2);
        double Qh2 = simpsonKF(a, b, n);
        double error = rungeError(Qh, Qh2, por2);
        double abs_error = abs(Qh2 - exactValue);
        
        cout << "| " << setw(14) << n;
        cout << " | " << setw(13) << h;
        cout << " | " << setw(13) << Qh2;
        cout << " | " << setw(13) << error;
        cout << " | " << setw(13) << abs_error << " |" << endl;
        
        err = error;
        n *= 2;
        
        if (n > 1000000) {
            break;
        }
    }
    cout << "----------------------------------------------------------------------------------" << endl;
    
    // Задание 2: КФ НАСТ с 4 узлами
    int k = 4;
    cout << endl;
    
    double gauss_result = gaussKF(a, b, k);
    double gauss_error = abs(gauss_result - exactValue);
    
    cout << "Приближенное значение: " << gauss_result << endl;
    cout << "Абсолютная погрешность: " << gauss_error << endl;
}
