#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double integrand(double x, double y) {
    return 2.0 * y * sin(x) + cos(x) * cos(x);
}

// --- Simpson's Rule Implementation ---
double simpson1D(const std::function<double(double)>& func, double a, double b, int n) {
    if (n % 2 != 0) {
        std::cerr << "Error (simpson1D): n must be even. Adjusting n to " << n + 1 << std::endl;
        n += 1; 
    }
    double h = (b - a) / n;
    double sum = func(a) + func(b);

    for (int i = 1; i < n; i += 2) { 
    }
    for (int i = 2; i < n; i += 2) { 
        sum += 2.0 * func(a + i * h);
    }
    return (h / 3.0) * sum;
}

double inner_integral_simpson(double x, int m) {
    double lower_y = sin(x);
    double upper_y = cos(x);
    std::function<double(double)> f_inner = [&](double y) {
        return integrand(x, y);
    };
     if (lower_y >= upper_y) {
         return 0.0;
     }
    return simpson1D(f_inner, lower_y, upper_y, m);
}

// 2D Simpson's Rule
double simpson2D(double a, double b, int n, int m) {
    std::function<double(double)> f_outer = [&](double x) {
        return inner_integral_simpson(x, m);
    };
    return simpson1D(f_outer, a, b, n);
}

// --- Gaussian Quadrature Implementation ---
const std::vector<double> gauss_points_3 = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
const std::vector<double> gauss_weights_3 = {5.0/9.0, 8.0/9.0, 5.0/9.0};

// 2D Gaussian Quadrature for variable limits
double gaussianQuadrature2D(double a, double b, int n, int m) {
    double integral = 0.0;

    const std::vector<double>& points_x = gauss_points_3; 
    const std::vector<double>& weights_x = gauss_weights_3;
    const std::vector<double>& points_y = gauss_points_3; 
    const std::vector<double>& weights_y = gauss_weights_3;


    double x_transform_factor = (b - a) / 2.0;
    double x_shift_factor = (b + a) / 2.0;

    for (int i = 0; i < n; ++i) { 
        double xi = points_x[i];
        double wi = weights_x[i];
        double x_val = x_transform_factor * xi + x_shift_factor;
        double g1_x = sin(x_val);
        double g2_x = cos(x_val);
         if (g1_x >= g2_x) {
             continue;
         }

        double y_transform_factor = (g2_x - g1_x) / 2.0;
        double y_shift_factor = (g2_x + g1_x) / 2.0;

        double inner_sum = 0.0;
        for (int j = 0; j < m; ++j) { 
            double nj = points_y[j];
            double wj = weights_y[j];
            double y_val = y_transform_factor * nj + y_shift_factor;

            inner_sum += wj * integrand(x_val, y_val);
        }
        integral += wi * y_transform_factor * inner_sum;
    }

    return x_transform_factor * integral;
}


int main() {
    double a = 0.0;
    double b = M_PI / 4.0;
    int n_simpson = 4;
    int m_simpson = 4;
    int n_gauss = 3; 
    int m_gauss = 3; 
    // a. Simpson's Rule
    double simpson_result = simpson2D(a, b, n_simpson, m_simpson);
    // b. Gaussian Quadrature
    double gauss_result = gaussianQuadrature2D(a, b, n_gauss, m_gauss);
    // c. Exact Value
    double exact_value = (5.0 * sqrt(2.0) - 4.0) / 6.0;
    // Output Results
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Approximating integral of (2y*sin(x) + cos^2(x)) dy dx" << std::endl;
    std::cout << "Limits: x from 0 to pi/4, y from sin(x) to cos(x)" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "a. Simpson's Rule (n=" << n_simpson << ", m=" << m_simpson << "): " << simpson_result << std::endl;
    std::cout << "b. Gaussian Quadrature (n=" << n_gauss << ", m=" << m_gauss << "):    " << gauss_result << std::endl;
    std::cout << "c. Exact Value:                         " << exact_value << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Comparison:" << std::endl;
    std::cout << "   Simpson Error:   " << std::abs(simpson_result - exact_value) << std::endl;
    std::cout << "   Gaussian Error:  " << std::abs(gauss_result - exact_value) << std::endl;


    return 0;
}