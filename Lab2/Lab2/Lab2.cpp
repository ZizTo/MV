#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

const double a = -3.0;
const double b = 3.0;

double f1(double x) {
    return sin(x) * cos(x);
}

double f2(double x) {
    return 1.0 / (1.0 + 12.0 * pow(x, 4));
}

vector<double> equidistant_nodes(int n, double (*f)(double)) {
    vector<double> nodes(n + 1);
    for (int i = 0; i <= n; ++i) {
        nodes[i] = a + i * (b - a) / n;
    }
    return nodes;
}

vector<double> chebyshev_nodes(int n, double (*f)(double)) {
    vector<double> nodes(n + 1);
    for (int i = 0; i <= n; ++i) {
        nodes[i] = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * M_PI / (2 * (n + 1)));
    }
    sort(nodes.begin(), nodes.end());
    return nodes;
}

vector<vector<double>> divided_differences(const vector<double>& nodes, double (*f)(double)) {
    int n = nodes.size();
    vector<vector<double>> table(n, vector<double>(n));

    for (int i = 0; i < n; ++i) {
        table[i][0] = f(nodes[i]);
    }

    for (int j = 1; j < n; ++j) {
        for (int i = 0; i < n - j; ++i) {
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (nodes[i + j] - nodes[i]);
        }
    }
    return table;
}

double newton_poly(double x, const vector<double>& nodes, const vector<vector<double>>& table) {
    double result = table[0][0];
    double product = 1.0;

    for (int i = 1; i < nodes.size(); ++i) {
        product *= (x - nodes[i - 1]);
        result += table[0][i] * product;
    }
    return result;
}


double printmax(vector<double> x_vals, vector<double> y_vals, double (*f)(double)) {
    double maxc = -1;
    for (int i = 0; i < x_vals.size(); ++i) {
        maxc = std::max(f(x_vals[i]) - y_vals[i], maxc);
    }
    return maxc;
}

int main() {
    vector<int> n_values = { 3, 5, 10, 15, 20, 30 };
    cout << "f1 - C" << endl;
    for (int n : n_values) {
        vector<double> nodes;
        nodes = chebyshev_nodes(n, f1);

        auto table = divided_differences(nodes, f1);

        vector<double> x_vals, y_vals;
        for (int i = 0; i <= 100; ++i) {
            double x = a + i * (b - a) / 100;
            x_vals.push_back(x);
            y_vals.push_back(newton_poly(x, nodes, table));
        }
        cout << n << " | " << printmax(x_vals, y_vals, f1) << endl;
    }
    cout << endl << endl << "f1 - P" << endl;
    for (int n : n_values) {
        vector<double> nodes;
        nodes = equidistant_nodes(n, f1);

        auto table = divided_differences(nodes, f1);

        vector<double> x_vals, y_vals;
        for (int i = 0; i <= 100; ++i) {
            double x = a + i * (b - a) / 100;
            x_vals.push_back(x);
            y_vals.push_back(newton_poly(x, nodes, table));
        }

        cout << n << " | " << printmax(x_vals, y_vals, f1) << endl;
    }
    cout << endl << endl << "f2 - C" << endl;
    for (int n : n_values) {
        vector<double> nodes;
        nodes = chebyshev_nodes(n, f2);

        auto table = divided_differences(nodes, f2);

        vector<double> x_vals, y_vals;
        for (int i = 0; i <= 100; ++i) {
            double x = a + i * (b - a) / 100;
            x_vals.push_back(x);
            y_vals.push_back(newton_poly(x, nodes, table));
        }

        cout << n << " | " << printmax(x_vals, y_vals, f2) << endl;
    }
    cout << endl << endl << "f2 - P" << endl;
    for (int n : n_values) {
        vector<double> nodes;
        nodes = equidistant_nodes(n, f2);

        auto table = divided_differences(nodes, f2);

        vector<double> x_vals, y_vals;
        for (int i = 0; i <= 100; ++i) {
            double x = a + i * (b - a) / 100;
            x_vals.push_back(x);
            y_vals.push_back(newton_poly(x, nodes, table));
        }

        cout << n << " | " << printmax(x_vals, y_vals, f2) << endl;
    }

    return 0;
}