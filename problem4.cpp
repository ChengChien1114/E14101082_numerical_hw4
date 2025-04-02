#include <iostream>
#include <cmath>
#include <functional>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double compositeSimpson1D(const std::function<double(double)>& func, double a, double b, int n) {
    if (n % 2 != 0 || n <= 0) return NAN;
    double h = (b - a) / n;
    double sum = func(a) + func(b);
    for (int i = 1; i < n; i += 2) sum += 4.0 * func(a + i * h);
    for (int i = 2; i < n - 1; i += 2) sum += 2.0 * func(a + i * h);
    return (h / 3.0) * sum;
}

double func_a(double x) {
    if (x == 0.0) return 0.0;
    if (x < 0.0) return NAN;
    return pow(x, -0.25) * sin(x);
}

double func_b_transformed(double t) {
    if (t == 0.0) return 0.0;
    if (t < 0.0) return NAN;
    return t * t * sin(1.0 / t);
}

int main() {
    int n = 4;
    double result_a = compositeSimpson1D(func_a, 0.0, 1.0, n);
    double result_b = compositeSimpson1D(func_b_transformed, 0.0, 1.0, n);

    if (!std::isnan(result_a)) std::cout  << "Approximate value (a): " << std::fixed << std::setprecision(10) << result_a << std::endl;
    if (!std::isnan(result_b)) std::cout << "Approximate value (b): "<< std::fixed << std::setprecision(10) << result_b << std::endl;

    return 0;
}
