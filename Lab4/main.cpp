/*#include <iostream>
#include <cmath>

using namespace std;

float eps = 1e-7;

float func(float x) {
    
}

float QsrPr(float h, int N)

void integ(float a, float b, float(*)(float, int)) {
    int N = 1;
    while () {
        N*=2;
        float h = (b-a)/N;
    }
} 

int main() {
    integ(0, M_PI/4, );


}

*/

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
    return (Qh - Qh2) / (pow(2, p) - 1.0);
}

void gaussQuadratureNodes4(vector<double>& nodes, vector<double>& weights) {
    // Узлы на [-1, 1]
    nodes.resize(4);
    nodes[0] = -sqrt((3.0 + 2.0 * sqrt(6.0/5.0)) / 7.0);
    nodes[1] = -sqrt((3.0 - 2.0 * sqrt(6.0/5.0)) / 7.0);
    nodes[2] = sqrt((3.0 - 2.0 * sqrt(6.0/5.0)) / 7.0);
    nodes[3] = sqrt((3.0 + 2.0 * sqrt(6.0/5.0)) / 7.0);
    
    // Веса
    weights.resize(4);
    weights[0] = (18.0 - sqrt(30.0)) / 36.0;
    weights[1] = (18.0 + sqrt(30.0)) / 36.0;
    weights[2] = (18.0 + sqrt(30.0)) / 36.0;
    weights[3] = (18.0 - sqrt(30.0)) / 36.0;
}

// Квадратурная формула НАСТ с k узлами (k = 4 для варианта 4)
double gaussQuadrature(double a, double b, int k) {
    vector<double> nodes, weights;
    
    switch(k) {
        case 4:
            gaussQuadratureNodes4(nodes, weights);
            break;
        default:
            cerr << "Unsupported number of nodes: " << k << endl;
            return 0.0;
    }
    
    double result = 0.0;
    // Преобразование отрезка [a, b] к [-1, 1]
    double mid = (b + a) / 2.0;
    double half_length = (b - a) / 2.0;
    
    for (int i = 0; i < k; i++) {
        double x = mid + half_length * nodes[i];  // Преобразование узла
        result += weights[i] * f(x);
    }
    
    return half_length * result;
}

int main() {
    double a = 0.0;          // Нижний предел интегрирования
    double b = M_PI / 4.0;     // Верхний предел интегрирования
    
    cout << fixed << setprecision(10);
    
    cout << "Задание 1" << endl;
    
    cout << "\nКвадратурная формула средних прямоугольников:" << endl;
    cout << "---------------------------------------------------------------------------------" << endl;
    cout << "| Число разбиений|      Шаг      |  Приближенное |   Оценка     |  Абсолютная   |" << endl;
    cout << "|       N        |       h       |   значение    |  погрешности |   погрешность |" << endl;
    cout << "---------------------------------------------------------------------------------" << endl;

    int por1 = 2;
    int n = 2;
    double prev_error = 1.0; // Начальное значение для погрешности
    
    while (abs(prev_error) > eps) {
        double h = (b - a) / n;
        double Qh = srPriamoug(a, b, n);
        double Qh2 = srPriamoug(a, b, 2*n);
        double error = rungeError(Qh, Qh2, por1);
        double abs_error = abs(Qh - exactValue);
        
        cout << "| " << setw(14) << n;
        cout << " | " << setw(13) << h;
        cout << " | " << setw(13) << Qh;
        cout << " | " << setw(13) << error;
        cout << " | " << setw(13) << abs_error << " |" << endl;
        
        prev_error = error;
        n *= 2;
        
        if (n > 1000000) {
            break;
        }
    }
    cout << "---------------------------------------------------------------------------------" << endl;
    
    cout << "\nКвадратурная формула Симпсона:" << endl;
    cout << "---------------------------------------------------------------------------------" << endl;
    cout << "| Число разбиений|      Шаг      |  Приближенное |   Оценка     |  Абсолютная   |" << endl;
    cout << "|       N        |       h       |   значение    |  погрешности |   погрешность |" << endl;
    cout << "---------------------------------------------------------------------------------" << endl;
    
    int por2 = 4;
    n = 2;
    prev_error = 1.0;
    
    while (abs(prev_error) > eps) {
        double h = (b - a) / n;
        double Qh = simpsonKF(a, b, n);
        double Qh2 = simpsonKF(a, b, 2*n);
        double error = rungeError(Qh, Qh2, por2);
        double abs_error = abs(Qh - exactValue);
        
        cout << "| " << setw(14) << n;
        cout << " | " << setw(13) << h;
        cout << " | " << setw(13) << Qh;
        cout << " | " << setw(13) << error;
        cout << " | " << setw(13) << abs_error << " |" << endl;
        
        prev_error = error;
        n *= 2;
        
        if (n > 1000000) {
            break;
        }
    }
    cout << "---------------------------------------------------------------------------------" << endl;
    
    // Задание 2: КФ НАСТ с k узлами
    int k = 4;
    cout << "\nЗадание 2: Квадратурная формула наивысшей алгебраической степени точности (НАСТ)" << endl;
    cout << "Количество узлов k = " << k << endl;
    
    double gauss_result = gaussQuadrature(a, b, k);
    double gauss_error = abs(gauss_result - exactValue);
    
    cout << "Приближенное значение: " << gauss_result << endl;
    cout << "Абсолютная погрешность: " << gauss_error << endl;
}
