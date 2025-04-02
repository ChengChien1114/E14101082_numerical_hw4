#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>


double func(double x) {
    return exp(x) * sin(4.0 * x);
}

// Composite Trapezoidal Rule
double compositeTrapezoidal(double a, double b, int n, double h) {
    double sum = func(a) + func(b);
    for (int i = 1; i < n; ++i) {
        double x_i = a + i * h;
        sum += 2.0 * func(x_i);
    }
    return (h / 2.0) * sum;
}

// Composite Simpson's Rule
double compositeSimpson(double a, double b, int N, double h) {
    if (N % 2 != 0) {
        std::cerr << "Error: N must be even for Composite Simpson's Rule." << std::endl;
        return NAN; 
    }
    double sum = func(a) + func(b);
    int n_pdf = N / 2; 

    for (int i = 1; i <= n_pdf; ++i) {
        double x_odd = a + (2 * i - 1) * h;
        sum += 4.0 * func(x_odd);
    }

    for (int i = 1; i < n_pdf; ++i) { 
        double x_even = a + (2 * i) * h;
        sum += 2.0 * func(x_even);
    }

    return (h / 3.0) * sum;
}


// Composite Midpoint Rule
double compositeMidpoint(double a, double b, int n, double h) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double midpoint = a + (i + 0.5) * h;
        sum += func(midpoint);
    }
    return h * sum;
}


int main() {
    double a = 1.0;
    double b = 2.0;
    double h = 0.1;
    int n = static_cast<int>(round((b - a) / h));
    if (std::abs(n * h - (b - a)) > 1e-9) {
         n = static_cast<int>((b - a) / h + 0.5);
         h = (b - a) / n; 
         std::cout << "Adjusted n to " << n << " and h to " << h << std::endl;
    } else {
        std::cout << "Using n = " << n << " and h = " << h << std::endl;
    }



    double trapezoidal_result = compositeTrapezoidal(a, b, n, h);
    double simpson_result = NAN;
     if (n % 2 == 0) {
        simpson_result = compositeSimpson(a, b, n, h);
    } else {
        std::cout << "N (" << n << ") is odd, Composite Simpson's rule requires an even number of intervals." << std::endl;
    }

    double midpoint_result = compositeMidpoint(a, b, n, h);

    // Output the results
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Integral of e^x * sin(4x) from " << a << " to " << b << " with h = " << h << ":" << std::endl;
    std::cout << "a. Composite Trapezoidal Rule Result: " << trapezoidal_result << std::endl;
     if (!std::isnan(simpson_result)) {
        std::cout << "b. Composite Simpson's Rule Result:   " << simpson_result << std::endl;
    }
    std::cout << "c. Composite Midpoint Rule Result:    " << midpoint_result << std::endl;

    return 0;
}