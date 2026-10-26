#include "LeastSquaresApproximation.h"

#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================ Конструктор ============================

LeastSquaresApproximation::LeastSquaresApproximation(
        const std::vector<double>& xValues,
        const std::vector<double>& yValues)
    : x(xValues), y(yValues)
{
    if (x.size() != y.size()) {
        throw std::invalid_argument("Размеры векторов x и y должны совпадать");
    }
    if (x.size() < 2) {
        throw std::invalid_argument("Требуется минимум 2 точки");
    }
}

// ============================ Вспомогательные ============================

bool LeastSquaresApproximation::isEquidistant(double& h) const {
    if (x.size() < 2) return false;
    h = x[1] - x[0];
    if (std::abs(h) < 1e-12) return false;
    for (size_t i = 1; i < x.size(); ++i) {
        if (std::abs((x[i] - x[i - 1]) - h) > 1e-10) return false;
    }
    return true;
}

double LeastSquaresApproximation::det3(double a11, double a12, double a13,
                                       double a21, double a22, double a23,
                                       double a31, double a32, double a33) {
    return a11 * (a22 * a33 - a23 * a32)
         - a12 * (a21 * a33 - a23 * a31)
         + a13 * (a21 * a32 - a22 * a31);
}

// ============================ Линейная аппроксимация ============================

void LeastSquaresApproximation::linearFit() {
    size_t n = x.size();
    double Sx = 0.0, Sy = 0.0, Sxy = 0.0, Sx2 = 0.0;
    for (size_t i = 0; i < n; ++i) {
        Sx  += x[i];
        Sy  += y[i];
        Sxy += x[i] * y[i];
        Sx2 += x[i] * x[i];
    }
    double denom = static_cast<double>(n) * Sx2 - Sx * Sx;
    if (std::abs(denom) < 1e-12) {
        throw std::runtime_error("Вырожденная система: невозможно построить линейный аппроксимант");
    }
    a1 = (static_cast<double>(n) * Sxy - Sx * Sy) / denom;
    a0 = (Sx2 * Sy - Sx * Sxy) / denom;
    currentFit = FitType::Linear;
}

// ============================ Квадратичная аппроксимация ============================

void LeastSquaresApproximation::quadraticFit() {
    size_t n = x.size();
    double Sx = 0.0, S2x = 0.0, S3x = 0.0, S4x = 0.0;
    double Sy = 0.0, Sxy = 0.0, S2xy = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double xi = x[i];
        double xi2 = xi * xi;
        Sx   += xi;
        S2x  += xi2;
        S3x  += xi2 * xi;
        S4x  += xi2 * xi2;
        Sy   += y[i];
        Sxy  += xi * y[i];
        S2xy += xi2 * y[i];
    }

    double nn = static_cast<double>(n);

    // Определитель основной матрицы системы (13)
    double Dz = det3(nn,  Sx,  S2x,
                     Sx,  S2x, S3x,
                     S2x, S3x, S4x);
    if (std::abs(Dz) < 1e-12) {
        throw std::runtime_error("Вырожденная система: невозможно построить квадратичный аппроксимант");
    }

    // Определители Крамера (формула (14))
    double D0 = det3(Sy,  Sx,  S2x,
                     Sxy, S2x, S3x,
                     S2xy, S3x, S4x);
    double D1 = det3(nn,  Sy,  S2x,
                     Sx,  Sxy, S3x,
                     S2x, S2xy, S4x);
    double D2 = det3(nn,  Sx,  Sy,
                     Sx,  S2x, Sxy,
                     S2x, S3x, S2xy);

    a0 = D0 / Dz;
    a1 = D1 / Dz;
    a2 = D2 / Dz;
    currentFit = FitType::Quadratic;
}

// ============================ Аппроксимация Фурье ============================

void LeastSquaresApproximation::fourierFit(int m, double a, double L) {
    size_t n = x.size();
    if (m < 1) {
        throw std::invalid_argument("Порядок суммы Фурье должен быть >= 1");
    }
    if (n < static_cast<size_t>(2 * m + 1)) {
        throw std::invalid_argument("Недостаточно узлов для выбранного порядка Фурье");
    }

    double h;
    if (!isEquidistant(h)) {
        throw std::invalid_argument("Для формулы (24)-(26) узлы должны быть равноотстоящими");
    }

    if (L <= 0.0) {
        L = x.back() - x.front();
    }

    // Сохраняем параметры преобразования аргумента
    fourierA = a;
    fourierL = L;
    fourierParamsSet = true;
    fourierOrder = m;

    fourierAk.assign(m + 1, 0.0);
    fourierBk.assign(m + 1, 0.0);

    // Постоянная составляющая a0 (24)
    double sumY = 0.0;
    for (size_t i = 0; i < n; ++i) sumY += y[i];
    fourierA0 = sumY / static_cast<double>(n);

    // Коэффициенты ak, bk (25)-(26)
    for (int k = 1; k <= m; ++k) {
        double sumC = 0.0, sumS = 0.0;
        for (size_t i = 0; i < n; ++i) {
            double t = 2.0 * M_PI * (x[i] - a) / L;
            sumC += y[i] * std::cos(k * t);
            sumS += y[i] * std::sin(k * t);
        }
        fourierAk[k] = 2.0 * sumC / static_cast<double>(n);
        fourierBk[k] = 2.0 * sumS / static_cast<double>(n);
    }

    currentFit = FitType::Fourier;
}

// ============================ Нелинейная: A0*(1-exp(-b*x)) ============================

void LeastSquaresApproximation::expSaturationFit() {
    size_t n = x.size();
    if (n < 2) {
        throw std::invalid_argument("Требуется минимум 2 точки");
    }

    // Вспомогательная функция F(b), корень которой – искомая b (31).
    auto F = [this, n](double b) -> double {
        double s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0;
        for (size_t i = 0; i < n; ++i) {
            double e = std::exp(-b * x[i]);
            double one_minus_e = 1.0 - e;
            s1 += y[i] * one_minus_e;
            s2 += one_minus_e * x[i] * e;
            s3 += x[i] * y[i] * e;
            s4 += one_minus_e * one_minus_e;
        }
        return s1 * s2 - s3 * s4;
    };

    // Поиск интервала смены знака для b > 0
    double bLo = 0.0, bHi = 1.0;
    double fLo = F(bLo), fHi = F(bHi);
    int iter = 0;
    while (fLo * fHi > 0.0 && iter < 100) {
        bHi *= 2.0;
        fHi = F(bHi);
        ++iter;
    }
    if (fLo * fHi > 0.0) {
        throw std::runtime_error("Не удалось локализовать корень для b (Пример 10)");
    }

    // Бисекция
    for (int it = 0; it < 200; ++it) {
        double bMid = 0.5 * (bLo + bHi);
        double fMid = F(bMid);
        if (std::abs(fMid) < 1e-14) { bLo = bHi = bMid; break; }
        if (fLo * fMid <= 0.0) { bHi = bMid; fHi = fMid; }
        else                    { bLo = bMid; fLo = fMid; }
    }
    double b = 0.5 * (bLo + bHi);

    // A0 из первого уравнения (30)
    double num = 0.0, den = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double one_minus_e = 1.0 - std::exp(-b * x[i]);
        num += y[i] * one_minus_e;
        den += one_minus_e * one_minus_e;
    }
    if (std::abs(den) < 1e-14) {
        throw std::runtime_error("Вырожденная система для A0 (Пример 10)");
    }

    expSatA0 = num / den;
    expSatB = b;
    currentFit = FitType::ExpSaturation;
}

// ============================ Нелинейная: функция Гаусса ============================

void LeastSquaresApproximation::gaussianFit(double A_init, double b_init, double c_init) {
    size_t n = x.size();
    if (n < 3) {
        throw std::invalid_argument("Для функции Гаусса требуется минимум 3 точки");
    }

    // Начальные приближения
    if (A_init <= 0.0) {
        double ymax = -std::numeric_limits<double>::infinity();
        for (double yi : y) if (yi > ymax) ymax = yi;
        A_init = (ymax > 0.0) ? ymax : 1.0;
    }
    if (b_init == 0.0) {
        double ymax = -std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < n; ++i) {
            if (y[i] > ymax) { ymax = y[i]; b_init = x[i]; }
        }
    }
    if (c_init <= 0.0) {
        double xmin = x.front(), xmax = x.back();
        c_init = (xmax - xmin) / 2.0;
        if (c_init <= 0.0) c_init = 1.0;
    }

    auto Q = [this, n](double A, double b, double c) -> double {
        double sum = 0.0;
        for (size_t i = 0; i < n; ++i) {
            double d = (x[i] - b) / c;
            double f = A * std::exp(-d * d);
            double diff = f - y[i];
            sum += diff * diff;
        }
        return sum;
    };

    double A = A_init, b = b_init, c = c_init;
    double Qcur = Q(A, b, c);

    // Покоординатный спуск с адаптивным шагом
    for (int outer = 0; outer < 200; ++outer) {
        double stepA = std::max(1e-3, std::abs(A) * 0.1);
        double stepB = std::max(1e-3, std::abs(c) * 0.1);
        double stepC = std::max(1e-3, std::abs(c) * 0.1);
        bool improved = false;

        // A
        for (int s = 0; s < 20; ++s) {
            double An = A + stepA;
            double Qn = Q(An, b, c);
            if (Qn < Qcur) { A = An; Qcur = Qn; improved = true; }
            else { stepA *= 0.5; if (stepA < 1e-12) break; }
        }
        // b
        for (int s = 0; s < 20; ++s) {
            double bn = b + stepB;
            double Qn = Q(A, bn, c);
            if (Qn < Qcur) { b = bn; Qcur = Qn; improved = true; }
            else { stepB *= 0.5; if (stepB < 1e-12) break; }
        }
        // c (должен оставаться > 0)
        for (int s = 0; s < 20; ++s) {
            double cn = c + stepC;
            if (cn <= 1e-9) { stepC *= 0.5; if (stepC < 1e-12) break; continue; }
            double Qn = Q(A, b, cn);
            if (Qn < Qcur) { c = cn; Qcur = Qn; improved = true; }
            else { stepC *= 0.5; if (stepC < 1e-12) break; }
        }

        if (!improved) break;
    }

    gaussA = A;
    gaussB = b;
    gaussC = c;
    currentFit = FitType::Gaussian;
}

// ============================ evaluate ============================

double LeastSquaresApproximation::evaluate(double xVal, double a, double L) const {
    switch (currentFit) {
        case FitType::Linear:
            return a1 * xVal + a0;

        case FitType::Quadratic:
            return a2 * xVal * xVal + a1 * xVal + a0;

        case FitType::Fourier: {
            double aa = std::isnan(a) ? fourierA : a;
            double LL = std::isnan(L) ? fourierL : L;

            if (!fourierParamsSet && (std::isnan(a) || std::isnan(L))) {
                throw std::runtime_error("Параметры Фурье-преобразования не заданы");
            }
            if (LL <= 0.0) LL = x.back() - x.front();

            double t = 2.0 * M_PI * (xVal - aa) / LL;
            double result = fourierA0;
            for (int k = 1; k <= fourierOrder; ++k) {
                result += fourierAk[k] * std::cos(k * t)
                        + fourierBk[k] * std::sin(k * t);
            }
            return result;
        }

        case FitType::ExpSaturation:
            return expSatA0 * (1.0 - std::exp(-expSatB * xVal));

        case FitType::Gaussian:
            return gaussA * std::exp(-std::pow((xVal - gaussB) / gaussC, 2.0));

        default:
            throw std::runtime_error("Аппроксимант ещё не построен");
    }
}

// ============================ residualSumSquares ============================

double LeastSquaresApproximation::residualSumSquares(double a, double L) const {
    double Q = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        double diff = evaluate(x[i], a, L) - y[i];
        Q += diff * diff;
    }
    return Q;
}

// ============================ exponentialFit ============================

void LeastSquaresApproximation::exponentialFit(const std::vector<double>& t,
                                               const std::vector<double>& A,
                                               double& A0, double& k) {
    if (t.size() != A.size() || t.size() < 2) {
        throw std::invalid_argument("Некорректные размеры векторов t и A");
    }

    std::vector<double> z(A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        if (A[i] <= 0.0) {
            throw std::invalid_argument("Амплитуды должны быть положительными");
        }
        z[i] = std::log(A[i]);
    }

    // Линейная аппроксимация z = a - k*t
    LeastSquaresApproximation ls(t, z);
    ls.linearFit();

    double aa, slope;
    ls.getLinearCoefficients(aa, slope);

    A0 = std::exp(aa);
    k  = -slope;
}

// ============================ Геттеры ============================

void LeastSquaresApproximation::getLinearCoefficients(double& a0_out, double& a1_out) const {
    if (currentFit != FitType::Linear) {
        throw std::runtime_error("Линейный аппроксимант не построен");
    }
    a0_out = a0;
    a1_out = a1;
}

void LeastSquaresApproximation::getQuadraticCoefficients(double& a0_out, double& a1_out,
                                                         double& a2_out) const {
    if (currentFit != FitType::Quadratic) {
        throw std::runtime_error("Квадратичный аппроксимант не построен");
    }
    a0_out = a0;
    a1_out = a1;
    a2_out = a2;
}

void LeastSquaresApproximation::getExpSaturationParams(double& A0_out, double& b_out) const {
    if (currentFit != FitType::ExpSaturation) {
        throw std::runtime_error("Аппроксимант A0*(1-exp(-b*x)) не построен");
    }
    A0_out = expSatA0;
    b_out  = expSatB;
}

void LeastSquaresApproximation::getGaussianParams(double& A_out, double& b_out, double& c_out) const {
    if (currentFit != FitType::Gaussian) {
        throw std::runtime_error("Аппроксимант Гаусса не построен");
    }
    A_out = gaussA;
    b_out = gaussB;
    c_out = gaussC;
}

// ============================ printResults ============================

void LeastSquaresApproximation::printResults(const std::string& title, double a, double L) const {
    std::cout << "=== " << title << " ===" << std::endl;
    std::cout << std::fixed << std::setprecision(4);

    std::cout << "Исходные данные:" << std::endl;
    for (size_t i = 0; i < x.size(); ++i) {
        std::cout << " (" << x[i] << ", " << y[i] << ")" << std::endl;
    }
    std::cout << std::endl;

    switch (currentFit) {
        case FitType::Linear:
            std::cout << "Линейный аппроксимант: P1(x) = a1*x + a0" << std::endl;
            std::cout << " a0 = " << a0 << std::endl;
            std::cout << " a1 = " << a1 << std::endl;
            break;

        case FitType::Quadratic:
            std::cout << "Квадратичный аппроксимант: P2(x) = a2*x^2 + a1*x + a0" << std::endl;
            std::cout << " a0 = " << a0 << std::endl;
            std::cout << " a1 = " << a1 << std::endl;
            std::cout << " a2 = " << a2 << std::endl;
            break;

        case FitType::Fourier:
            std::cout << "Сумма Фурье порядка m = " << fourierOrder << std::endl;
            std::cout << " a0 = " << fourierA0 << std::endl;
            for (int k = 1; k <= fourierOrder; ++k) {
                std::cout << " a" << k << " = " << fourierAk[k]
                          << ", b" << k << " = " << fourierBk[k] << std::endl;
            }
            break;

        case FitType::ExpSaturation:
            std::cout << "Аппроксимант A0*(1-exp(-b*x))" << std::endl;
            std::cout << " A0 = " << expSatA0 << std::endl;
            std::cout << " b  = " << expSatB  << std::endl;
            break;

        case FitType::Gaussian:
            std::cout << "Аппроксимант Гаусса A*exp(-((x-b)/c)^2)" << std::endl;
            std::cout << " A = " << gaussA << std::endl;
            std::cout << " b = " << gaussB << std::endl;
            std::cout << " c = " << gaussC << std::endl;
            break;

        default:
            std::cout << "Аппроксимант не построен" << std::endl;
            return;
    }

    std::cout << std::endl;
    std::cout << "Значения аппроксиманта в узлах:" << std::endl;
    std::cout << " x\t\ty\t\tF(x)\t\t|F(x)-y|" << std::endl;
    for (size_t i = 0; i < x.size(); ++i) {
        double fx = evaluate(x[i], a, L);
        std::cout << " " << x[i] << "\t\t" << y[i] << "\t\t" << fx
                  << "\t\t" << std::abs(fx - y[i]) << std::endl;
    }
    std::cout << "Сумма квадратов отклонений Q = "
              << residualSumSquares(a, L) << std::endl;
    std::cout << std::endl;
}