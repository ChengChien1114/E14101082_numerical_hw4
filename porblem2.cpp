#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
using namespace std;

double f(double x) {
    return x * x * log(x);
}

double gaussQuad(double a, double b, int n) {
    vector<double> xi, ci;

    if (n == 3) {
        xi = {-sqrt(3.0/5), 0.0, sqrt(3.0/5)};
        ci = {5.0/9, 8.0/9, 5.0/9};
    } else if (n == 4) {
        xi = {-0.861136, -0.339981, 0.339981, 0.861136};
        ci = {0.347855, 0.652145, 0.652145, 0.347855};
    } else {
        cerr << "Only n = 3 or 4 supported.\n";
        return 0;
    }

    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double x = (b - a) / 2 * xi[i] + (a + b) / 2;
        sum += ci[i] * f(x);
    }

    return (b - a) / 2 * sum;
}

double exact_integral() {
    auto F = [](double x) {
        return (x * x * x * log(x)) / 3.0 - (x * x * x) / 9.0;
    };
    return F(1.5) - F(1.0);
}

int main() {
    double a = 1.0, b = 1.5;

    double result3 = gaussQuad(a, b, 3);
    double result4 = gaussQuad(a, b, 4);
    double exact = exact_integral();

    cout << fixed << setprecision(10);
    cout << "Gaussian Quadrature (n=3): " << result3 << endl;
    cout << "Gaussian Quadrature (n=4): " << result4 << endl;
    cout << "Exact Value:               " << exact << endl;
    cout << "Error (n=3):               " << fabs(result3 - exact) << endl;
    cout << "Error (n=4):               " << fabs(result4 - exact) << endl;

    return 0;
}
