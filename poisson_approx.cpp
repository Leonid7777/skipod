#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <string>

class PoissonSolver
{
    int M, N;
    double h1, h2, eps, delta;
    double h1_square, h2_square;
    double A1, A2, B1, B2;
    double alpha, betta;
    double zr_k, zr_k_prev, app;
    std::vector<double> r, z, p, w, Ap, Q_arr, F_arr;
    
    double S_ij (int i, int j)
    {
        double x = std::abs(A1 + i * h1);
        double y = std::abs(A2 + j * h2);

        if(x >= h1 / 2 && x <= 2 - h1 / 2) 
        {
            if(y - h2 / 2 >= 1 + x + h1 / 2) {
                return 0;
            } else if(y + h2 / 2 <= 1 + x - h1 / 2) {
                return h1 * h2;
            } else if(y - h2 / 2 <= 1 + x - h1 / 2) {
                double l = y + h2 / 2 - 1 - x + h1 / 2;
                return h1 * h2 - l * l / 2;
            } else if(y + h2 / 2 <= 1 + x + h1 / 2) {
                return (x - y + 1 + h1 / 2) * h2;
            } else {
                double l = x + h1 / 2 - y + h2 / 2 + 1;
                return l * l / 2;
            }
        }
        else if(x >= 2 + h1 / 2 && x <= 5 - h1 / 2)
        {
            if(y - h2 / 2 >= 5 - x + h1 / 2) {
                return 0;
            } else if(y + h2 / 2 <= 5 - x - h1 / 2) {
                return h1 * h2;
            } else if(y - h2 / 2 <= 5 - x - h1 / 2) {
                double l = y + h2 / 2 - 5 + x + h1 / 2;
                return h1 * h2 - l * l / 2;
            } else if(y + h2 / 2 <= 5 - x + h1 / 2) {
                return (5 - x - y + h1 / 2) * h2;
            } else {
                double l = 5 - x + h1 / 2 - y + h2 / 2;
                return l * l / 2;
            }
        }
        else if(x < h1 / 2) 
        {
            if(y + h2 / 2 <= 1) {
                return h1 * h2;
            } else if (y >= 1 + h1 / 2 + h2 / 2) {
                return 0;
            } else if (y + h2 / 2 >= 1 && y + h2 / 2 <= 1 + h1 / 2) {
                double l = y + h2 / 2 - 1;
                return h1 * h2 - l * l;
            } else if (y - h2 / 2 <= 1) {
                double l1 = 1 - y + h2 / 2;
                double l2 = 1 + h1 / 2 + h2 / 2 - y;
                return (l1 + l2) * h1 / 2;
            } else {
                double l = 1 + h1 / 2 + h2 / 2 - y;
                return l * l;
            }
        }
        else if(x < 2 + h1 / 2)
        {
            if(y + h2 / 2 <= 3 - h1 / 2) {
                return h1 * h2;
            } else if (y - h2 / 2 >= 3) {
                return 0;
            } else if (y + h2 / 2 >= 3 - h1 / 2 && y + h2 / 2 <= 3) {
                double l = y + h2 / 2 - 3 + h1 / 2;
                return h1 * h2 - l * l;
            } else if (y - h2 / 2 <= 3 - h1 / 2) {
                double l1 = 3 - y + h2 / 2;
                double l2 = 3 - h1 / 2 - y + h2 / 2;
                return (l1 + l2) * h1 / 2;
            } else {
                double l = 3 - y + h2 / 2;
                return l * l;
            }
        }
        return 0;
    }

    void values_Q_F() {
        for(int j = 1; j < N; j++) {
            for(int i = 1; i < M; i++) {
                double S = S_ij(i, j);
                Q_arr[j * (M + 1) + i] = (h1 * h2 - S) / (h1 * h2 * eps * eps);
                F_arr[j * (M + 1) + i] = S / (h1 * h2);
            }
        }
    }

    double mat_vec_Ap_ij(const std::vector<double>& p, int i, int j)
    {
        int abs_index = j * (M + 1) + i;
        double Q = Q_arr[abs_index];
        return -p[abs_index + 1] / h1_square + p[abs_index] * (2 / h1_square + 2 / h2_square + Q) - p[abs_index - 1] / h1_square - p[abs_index + M + 1] / h2_square - p[abs_index - M - 1] / h2_square;
    }

    void mat_vec_Ap(const std::vector<double>& p)
    {
        for(int j = 1; j < N; j++) {
            for(int i = 1; i < M; i++) {
                Ap[j * (M + 1) + i] = mat_vec_Ap_ij(p, i, j);
            }
        }
    }

    double right_part_ij(int i, int j)
    {
        int abs_index = j * (M + 1) + i;
        return F_arr[abs_index] + Q_arr[abs_index] * std::sin(A1 + i * h1 + A2 + j * h2);
    }

    void right_part()
    {
        for(int j = 1; j < N; j++) {
            for(int i = 1; i < M; i++) {
                r[j * (M + 1) + i] = right_part_ij(i, j);
            }
        }
    }

    void z_from_r()
    {
        for(int j = 1; j < N; j++) {
            for(int i = 1; i < M; i++) {
                int abs_index = j * (M + 1) + i;
                z[abs_index] = r[abs_index] / (2 / h1_square + 2 / h2_square + Q_arr[abs_index]);
            }
        }
    }

    void update_r() 
    {
        int abs_index;
        for(int j = 1; j < N; j++) {
            for(int i = 1; i < M; i++) {
                abs_index = j * (M + 1) + i;
                r[abs_index] -= alpha * Ap[abs_index];
            }
        }
    }

    void dot_prod_Ap_p()
    {
        app = 0;
        int abs_index;
        for(int j = 1; j < N; j++) {
            for(int i = 1; i < M; i++) {
                abs_index = j * (M + 1) + i;
                app += Ap[abs_index] * p[abs_index];
            }
        }
    }

    void dot_prod_z_r()
    {
        zr_k = 0;
        int abs_index;
        for(int j = 1; j < N; j++) {
            for(int i = 1; i < M; i++) {
                abs_index = j * (M + 1) + i;
                zr_k += z[abs_index] * r[abs_index];
            }
        }
    }

    void update_p()
    {
        int abs_index;
        for(int j = 1; j < N; j++) {
            for(int i = 1; i < M; i++) {
                abs_index = j * (M + 1) + i;
                p[abs_index] = z[abs_index] + betta * p[abs_index];
            }
        }
    }

    void update_w()
    {
        int abs_index;
        for(int j = 1; j < N; j++) {
            for(int i = 1; i < M; i++) {
                abs_index = j * (M + 1) + i;
                w[abs_index] += alpha * p[abs_index];
            }
        }
    }

    double error_w()
    {
        double error = 0;
        int abs_index;
        for(int j = 1; j < N; j++) {
            for(int i = 1; i < M; i++) {
                abs_index = j * (M + 1) + i;
                error += p[abs_index] * p[abs_index];
            }
        }
        return alpha * std::sqrt(error * h1 * h2);
    }

public:
    PoissonSolver(int M_, int N_, double A1_, double B1_, double A2_, double B2_, double delta_)
    {
        M = M_;
        N = N_;
        A1 = A1_;
        B1 = B1_;
        A2 = A2_;
        B2 = B2_;
        delta = delta_;
        h1 = (B1 - A1) / M;
        h2 = (B2 - A2) / N;
        h1_square = h1 * h1;
        h2_square = h2 * h2;
        eps = std::max(h1_square, h2_square);
        r.resize((M + 1) * (N + 1), 0);
        z.resize((M + 1) * (N + 1), 0);
        p.resize((M + 1) * (N + 1), 0);
        Ap.resize((M + 1) * (N + 1), 0);
        w.resize((M + 1) * (N + 1), 0);
        Q_arr.resize((M + 1) * (N + 1), 0);
        F_arr.resize((M + 1) * (N + 1), 0);

        values_Q_F();
    }

    int CG()
    {
        right_part();

        for(int i = 0; i <= M; i++) {
            w[i] = std::sin(A1 + i * h1 + A2);
            w[N * (M + 1) + i] = std::sin(A1 + i * h1 + B2);
        }
        for(int j = 0; j <= N; j++) {
            w[j * (M + 1)] = std::sin(A1 + A2 + j * h2);
            w[j * (M + 1) + M] = std::sin(B1 + A2 + j * h2);
        }

        mat_vec_Ap(w);

        for(int j = 1; j < N; j++) {
            for(int i = 1; i < M; i++) {
                int abs_index = j * (M + 1) + i;
                r[abs_index] -= Ap[abs_index];
            }
        }

        z_from_r();

        for(int j = 1; j < N; j++) {
            for(int i = 1; i < M; i++) {
                int abs_index = j * (M + 1) + i;
                p[abs_index] = z[abs_index];
            }
        }
        dot_prod_z_r();
        mat_vec_Ap(p);
        dot_prod_Ap_p();
        alpha = zr_k / app;

        for(int j = 1; j < N; j++) {
            for(int i = 1; i < M; i++) {
                int abs_index = j * (M + 1) + i;
                w[abs_index] = alpha * p[abs_index];
            }
        }

        zr_k_prev = zr_k;

        int iter = 1;
        while (error_w() >= delta && iter < (N - 1) * (M - 1))
        {
            update_r();
            z_from_r();
            dot_prod_z_r();
            betta = zr_k / zr_k_prev;
            update_p();
            mat_vec_Ap(p);
            dot_prod_Ap_p();
            alpha = zr_k / app;
            update_w();
            zr_k_prev = zr_k;
            ++iter;
        }
        return iter;
    }

    const std::vector<double>& get_w() const
    {
        return w;
    }

    double relative_error(const std::vector<double>& v)
    {
        double res = 0;
        double norm_B = 0;
        for(int j = 1; j < N; j++) {
            for(int i = 1; i < M; i++) {
                double B = right_part_ij(i, j);
                double diff = B - mat_vec_Ap_ij(v, i, j);
                res += diff * diff;
                norm_B += B * B;
            }
        }
        return std::sqrt(res / norm_B);
    }

    void save_w(const std::string& filename) const
    {
        std::ofstream out(filename);
        for(int j = 0; j <= N; j++) {
            for(int i = 0; i <= M; i++) {
                out << A1 + i * h1 << " " << A2 + j * h2 << " " << w[j * (M + 1) + i] << "\n";
            }
            out << "\n";
        }
    }
};

int main(void)
{
    int M, N;
    std::cin >> M >> N;

    PoissonSolver Ps(M, N, -5, 5, -3, 3, 1e-6);

    auto start = std::chrono::steady_clock::now();
    int iter = Ps.CG();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "N_iter = " << iter << std::endl;
    std::cout << "Time elapsed: " << elapsed.count() << " ms\n";

    const std::vector<double>& w = Ps.get_w();
    std::cout << "Error = " << Ps.relative_error(w) << std::endl;

    Ps.save_w("solution_" + std::to_string(M) + "x" + std::to_string(N) + ".txt");
    return 0;
}