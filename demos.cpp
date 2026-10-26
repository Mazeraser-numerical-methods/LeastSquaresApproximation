#include "demos.h"
#include "LeastSquaresApproximation.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

// ============================ Пример 8 ============================

void demoExample8() {
    std::cout << "================================" << std::endl;
    std::cout << " Пример 8: Аппроксимация полиномами" << std::endl;
    std::cout << "================================" << std::endl;

    std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<double> y = {0.9, 4.2, 5.9, 7.2, 7.6, 7.4, 6.9, 6.2, 4.5, 2.4};

    LeastSquaresApproximation ls(x, y);

    // Линейная аппроксимация
    ls.linearFit();
    ls.printResults("Линейная аппроксимация (Пример 8)");
    std::cout << "Ожидаемые коэффициенты: a0 = 4.787, a1 = 0.097"
              << std::endl << std::endl;

    // Квадратичная аппроксимация
    ls.quadraticFit();
    ls.printResults("Квадратичная аппроксимация (Пример 8)");
    std::cout << "Ожидаемые коэффициенты: a0 = -1.6467, a1 = 3.3136, a2 = -0.2924"
              << std::endl << std::endl;
}

// ============================ Пример 9 ============================

void demoExample9() {
    std::cout << "================================" << std::endl;
    std::cout << " Пример 9: Аппроксимация суммой Фурье" << std::endl;
    std::cout << "================================" << std::endl;

    std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<double> y = {1.3, 2.8, 2.5, 1.1, -0.8, -1.9, -1.1, 0.9, 3.85, 0.0};

    LeastSquaresApproximation ls(x, y);
    int m = 6;
    double a = 0.0, L = 10.0;
    ls.fourierFit(m, a, L);
    ls.printResults("Сумма Фурье m = 6 (Пример 9)", a, L);

    std::cout << "Ожидаемые коэффициенты (Таблица 11):" << std::endl;
    std::cout << " a0 = 0.86" << std::endl;
    std::cout << " a1 = 1.257, b1 = 1.105" << std::endl;
    std::cout << " a2 = -0.719, b2 = -1.246" << std::endl;
    std::cout << " a3 = -0.577, b3 = -0.551" << std::endl;
    std::cout << " a4 = -0.541, b4 = -0.323" << std::endl;
    std::cout << " a5 = -0.56, b5 = 0.0" << std::endl;
    std::cout << " a6 = -0.541, b6 = 0.551" << std::endl;
    std::cout << std::endl;

    // Сравнение с m = 2
    LeastSquaresApproximation ls2(x, y);
    ls2.fourierFit(2, a, L);

    std::cout << "Для сравнения – сумма Фурье порядка m = 2:" << std::endl;
    std::cout << " Q(m=2) = " << ls2.residualSumSquares(a, L) << std::endl;
    std::cout << " Q(m=6) = " << ls.residualSumSquares(a, L) << std::endl;
    std::cout << std::endl;
}

// ============================ Экспоненциальная зависимость ============================

void demoExponential() {
    std::cout << "================================" << std::endl;
    std::cout << " Демонстрация: A(t) = A0 * exp(-k*t)" << std::endl;
    std::cout << "================================" << std::endl;

    double A0_true = 10.0, k_true = 0.5;
    std::vector<double> t, A;
    for (int i = 0; i <= 10; ++i) {
        double ti = 0.2 * i;
        t.push_back(ti);
        A.push_back(A0_true * std::exp(-k_true * ti));
    }

    double A0, k;
    LeastSquaresApproximation::exponentialFit(t, A, A0, k);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Истинные значения: A0 = " << A0_true << ", k = " << k_true << std::endl;
    std::cout << "Найденные значения: A0 = " << A0 << ", k = " << k << std::endl;
    std::cout << std::endl;
}

// ============================ Пример 10 ============================

void demoExample10() {
    std::cout << "==============================================" << std::endl;
    std::cout << " Пример 10: y = A0*(1 - exp(-b*x))" << std::endl;
    std::cout << "==============================================" << std::endl;

    // Синтетические данные: A0 = 5, b = 0.8
    double A0_true = 5.0, b_true = 0.8;
    std::vector<double> x, y;
    for (int i = 0; i <= 10; ++i) {
        double xi = 0.2 * i;
        x.push_back(xi);
        y.push_back(A0_true * (1.0 - std::exp(-b_true * xi)));
    }

    LeastSquaresApproximation ls(x, y);
    ls.expSaturationFit();
    ls.printResults("Нелинейная аппроксимация A0*(1-exp(-b*x)) (Пример 10)");

    double A0, b;
    ls.getExpSaturationParams(A0, b);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Истинные значения: A0 = " << A0_true << ", b = " << b_true << std::endl;
    std::cout << "Найденные значения: A0 = " << A0 << ", b = " << b << std::endl;
    std::cout << std::endl;
}

// ============================ Пример 11 ============================

void demoExample11() {
    std::cout << "==============================================" << std::endl;
    std::cout << " Пример 11: y = A*exp(-((x-b)/c)^2)" << std::endl;
    std::cout << "==============================================" << std::endl;

    // Синтетические данные: A = 3, b = 2, c = 1
    double A_true = 3.0, b_true = 2.0, c_true = 1.0;
    std::vector<double> x, y;
    for (int i = 0; i <= 20; ++i) {
        double xi = -1.0 + 0.2 * i;
        x.push_back(xi);
        double d = (xi - b_true) / c_true;
        y.push_back(A_true * std::exp(-d * d));
    }

    LeastSquaresApproximation ls(x, y);
    ls.gaussianFit(A_true, b_true, c_true);
    ls.printResults("Нелинейная аппроксимация Гаусса (Пример 11)");

    double A, b, c;
    ls.getGaussianParams(A, b, c);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Истинные значения: A = " << A_true
              << ", b = " << b_true << ", c = " << c_true << std::endl;
    std::cout << "Найденные значения: A = " << A
              << ", b = " << b << ", c = " << c << std::endl;
    std::cout << std::endl;
}