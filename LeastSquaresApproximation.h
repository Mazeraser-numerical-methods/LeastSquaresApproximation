#ifndef LEAST_SQUARES_APPROXIMATION_H
#define LEAST_SQUARES_APPROXIMATION_H

#include <vector>
#include <string>
#include <limits>

/**
 * @brief Класс для аппроксимации табличных данных методом наименьших квадратов.
 *
 * Реализует:
 *  - линейную аппроксимацию P1(x) = a1*x + a0;
 *  - квадратичную аппроксимацию P2(x) = a2*x^2 + a1*x + a0;
 *  - аппроксимацию суммами Фурье (для равноотстоящих узлов);
 *  - сведение экспоненциальной зависимости к линейной;
 *  - нелинейную аппроксимацию A0*(1 - exp(-b*x)) (Пример 10);
 *  - нелинейную аппроксимацию функцией Гаусса A*exp(-((x-b)/c)^2) (Пример 11).
 */
class LeastSquaresApproximation {
public:
    /**
     * @brief Конструктор по таблично заданной функции.
     * @param xValues узлы (аргументы)
     * @param yValues значения функции в узлах
     * @throws std::invalid_argument при несовпадении размеров или < 2 точек
     */
    LeastSquaresApproximation(const std::vector<double>& xValues,
                              const std::vector<double>& yValues);

    /// Линейная аппроксимация P1(x) = a1*x + a0 (формулы (10)-(11)).
    void linearFit();

    /// Квадратичная аппроксимация P2(x) = a2*x^2 + a1*x + a0 (формулы (13)-(14)).
    void quadraticFit();

    /**
     * @brief Аппроксимация суммой Фурье T_m(t) = a0 + sum[ak*cos(k*t) + bk*sin(k*t)].
     * @param m порядок суммы Фурье (m >= 1)
     * @param a начало интервала (по умолчанию 0)
     * @param L длина интервала (по умолчанию определяется по узлам)
     */
    void fourierFit(int m, double a = 0.0, double L = 0.0);

    /**
     * @brief Нелинейная аппроксимация y = A0*(1 - exp(-b*x)) (Пример 10).
     * Уравнение для b (31) решается методом бисекции; A0 — из (30).
     */
    void expSaturationFit();

    /**
     * @brief Нелинейная аппроксимация функцией Гаусса (Пример 11).
     * Покоординатный спуск с адаптивным шагом.
     */
    void gaussianFit(double A_init = 0.0, double b_init = 0.0, double c_init = 1.0);

    /**
     * @brief Вычисляет значение аппроксиманта в точке xVal.
     * Для Фурье-аппроксиманта параметры a и L берутся из сохранённых
     * при вызове fourierFit, если не переданы явно.
     */
    double evaluate(double xVal,
                    double a = std::numeric_limits<double>::quiet_NaN(),
                    double L = std::numeric_limits<double>::quiet_NaN()) const;

    /// Остаточная сумма квадратов Q (формула (3)).
    double residualSumSquares(double a = std::numeric_limits<double>::quiet_NaN(),
                              double L = std::numeric_limits<double>::quiet_NaN()) const;

    /**
     * @brief Сведение экспоненциальной зависимости A(t) = A0*exp(-k*t) к линейной.
     * ln A = ln A0 - k*t.
     * @param t вектор моментов времени
     * @param A вектор измеренных амплитуд (должны быть > 0)
     * @param[out] A0 начальная амплитуда
     * @param[out] k коэффициент затухания
     */
    static void exponentialFit(const std::vector<double>& t,
                               const std::vector<double>& A,
                               double& A0, double& k);

    /// Возвращает коэффициенты линейного аппроксиманта.
    void getLinearCoefficients(double& a0_out, double& a1_out) const;

    /// Возвращает коэффициенты квадратичного аппроксиманта.
    void getQuadraticCoefficients(double& a0_out, double& a1_out, double& a2_out) const;

    /// Возвращает параметры аппроксиманта A0*(1-exp(-b*x)).
    void getExpSaturationParams(double& A0_out, double& b_out) const;

    /// Возвращает параметры аппроксиманта Гаусса.
    void getGaussianParams(double& A_out, double& b_out, double& c_out) const;

    /// Выводит результаты аппроксимации.
    void printResults(const std::string& title,
                      double a = std::numeric_limits<double>::quiet_NaN(),
                      double L = std::numeric_limits<double>::quiet_NaN()) const;

private:
    // Исходные данные
    std::vector<double> x;
    std::vector<double> y;

    // Коэффициенты полиномиальных аппроксимантов
    double a0 = 0.0, a1 = 0.0, a2 = 0.0;

    // Коэффициенты Фурье
    double fourierA0 = 0.0;
    std::vector<double> fourierAk;
    std::vector<double> fourierBk;
    int fourierOrder = 0;

    // Параметры нелинейных аппроксимантов
    double expSatA0 = 0.0, expSatB = 0.0;     // Пример 10
    double gaussA = 0.0, gaussB = 0.0, gaussC = 1.0; // Пример 11

    // Параметры преобразования аргумента для Фурье
    double fourierA = 0.0;
    double fourierL = 1.0;
    bool fourierParamsSet = false;

    enum class FitType {
        None, Linear, Quadratic, Fourier, ExpSaturation, Gaussian
    };
    FitType currentFit = FitType::None;

    /// Проверяет, являются ли узлы равноотстоящими; возвращает шаг h.
    bool isEquidistant(double& h) const;

    /// Определитель 3x3 (для метода Крамера в quadraticFit).
    static double det3(double a11, double a12, double a13,
                       double a21, double a22, double a23,
                       double a31, double a32, double a33);
};

#endif // LEAST_SQUARES_APPROXIMATION_H