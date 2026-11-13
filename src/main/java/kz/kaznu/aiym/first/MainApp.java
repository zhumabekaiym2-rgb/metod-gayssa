package kz.kaznu.aiym.first;

/**
 * Метод Гаусса с частичным выбором главного элемента (Java 17, без библиотек).
 * Решает квадратную СЛАУ A·x = b размера n.
 */
public class MainApp {

    /** Порог: ведущий элемент по модулю меньше EPS считаем нулевым. */
    private static final double EPS = 1e-12;

    /** Результат: верхнетреугольная расширенная матрица и решение. */
    public record Result(double[][] upper, double[] x) { }

    /** Строит расширенную матрицу [A | b] (копия, исходные данные не портим). */
    static double[][] augment(double[][] a, double[] b) {
        int n = a.length;
        double[][] m = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) m[i][j] = a[i][j];
            m[i][n] = b[i];
        }
        return m;
    }

    static Result solve(double[][] a, double[] b) {
        int n = a.length;
        double[][] m = augment(a, b);

        // ---------- ПРЯМОЙ ХОД ----------
        for (int k = 0; k < n; k++) {                 // k — номер шага (ведущий столбец)
            int p = k;                                // поиск главного элемента в столбце k
            for (int i = k + 1; i < n; i++)
                if (Math.abs(m[i][k]) > Math.abs(m[p][k])) p = i;

            if (Math.abs(m[p][k]) < EPS)              // весь столбец ниже диагонали ≈ 0
                throw new ArithmeticException("Матрица вырождена или почти вырождена: "
                        + "ведущий элемент в столбце " + (k + 1) + " равен нулю");

            if (p != k) {                             // R_k <-> R_p (обмен ссылок на строки)
                double[] tmp = m[k]; m[k] = m[p]; m[p] = tmp;
            }

            for (int i = k + 1; i < n; i++) {         // строки под ведущей
                double mult = m[i][k] / m[k][k];      // m_ik = a_ik / a_kk
                for (int j = k; j <= n; j++)          // R_i := R_i - m_ik * R_k
                    m[i][j] -= mult * m[k][j];
                m[i][k] = 0.0;                        // убираем «хвост» округления
            }
        }

        // ---------- ОБРАТНЫЙ ХОД ----------
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = m[i][n];                     // b_i
            for (int j = i + 1; j < n; j++)
                sum -= m[i][j] * x[j];                // минус вклад найденных x_j
            x[i] = sum / m[i][i];                     // деление на a_ii
        }
        return new Result(m, x);
    }

    /** Невязка r = b - A·x и её максимальная по модулю компонента. */
    static double residual(double[][] a, double[] b, double[] x) {
        double max = 0.0;
        for (int i = 0; i < a.length; i++) {
            double ax = 0.0;
            for (int j = 0; j < a.length; j++) ax += a[i][j] * x[j];
            max = Math.max(max, Math.abs(b[i] - ax));
        }
        return max;
    }

    static void print(Result r, double[][] a, double[] b, String[] names) {
        int n = a.length;
        System.out.println("Верхнетреугольная расширенная матрица:");
        for (double[] row : r.upper()) {
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < n; j++) sb.append(String.format("%10.6f", row[j]));
            sb.append(String.format(" | %10.6f", row[n]));
            System.out.println(sb);
        }
        System.out.println("Решение:");
        for (int i = 0; i < n; i++)
            System.out.printf("%s = %.6f%n", names[i], r.x()[i]);
        System.out.printf("Невязка max|b - A·x| = %.3e%n", residual(a, b, r.x()));
    }

    public static void main(String[] args) {
        // Основная задача из раздела IV
        double[][] a = { { 2, 1, -1 }, { -2, 2, 1 }, { 1, -1, 1 } };
        double[] b = { 8, 1, -2 };
        print(solve(a, b), a, b, new String[] { "x", "y", "z" });

        // Нулевой ведущий элемент: y = 2, x + y = 3
        System.out.println();
        double[][] a2 = { { 0, 1 }, { 1, 1 } };
        double[] b2 = { 2, 3 };
        print(solve(a2, b2), a2, b2, new String[] { "x", "y" });

        // Вырожденная матрица: строки пропорциональны
        System.out.println();
        try {
            solve(new double[][] { { 1, 2 }, { 2, 4 } }, new double[] { 3, 6 });
        } catch (ArithmeticException e) {
            System.out.println("Ошибка: " + e.getMessage());
        }
    }
}