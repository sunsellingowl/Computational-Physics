#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>

#define PI 3.14159265358
using namespace std;

typedef complex<double> Complex;
// 一维FFT
vector<Complex> FFT(vector<Complex>& x) {
    int n = x.size();
    if (n <= 1) return x;

    vector<Complex> even, odd;
    for (int i = 0; i < n; i++) {
        if (i % 2 == 0) {
            even.push_back(x[i]);
        }
        else {
            odd.push_back(x[i]);
        }
    }

    even = FFT(even);
    odd = FFT(odd);

    vector<Complex> X(n);
    for (int k = 0; k < n / 2; k++) {
        Complex t = polar(1.0, -2 * PI * k / n) * odd[k];
        X[k] = even[k] + t;
        X[k + n / 2] = even[k] - t;
    }

    return X;
}
// 二维FFT
vector<vector<Complex>> FFT2D(vector<vector<Complex>>& f) {
    int row = f.size();    // 行数
    int col = f[0].size(); // 列数

    // 对每一行进行一维FFT
    for (int i = 0; i < row; i++) {
        f[i] = FFT(f[i]);
    }

    // 对每一列进行一维FFT
    vector<vector<Complex>> result(col, vector<Complex>(row));
    for (int j = 0; j < col; j++) {
        vector<Complex> column0(row);
        vector<Complex> column(row);
        for (int i = 0; i < row; i++) {
            column0[i] = f[i][j];
        }
        column = FFT(column0);
        for (int i = 0; i < row; i++) {
            result[i][j] = column[i];
        }
    }

    return result;
}

// 输出
void printMatrix(vector<vector<Complex>> matrix, int N, int M) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            std::cout << "(" << matrix[i][j].real() << ", " << matrix[i][j].imag() << ") " << '\t';
        }
        std::cout << std::endl;
    }
}
void printMatrix(vector<vector<Complex>> matrix, int N, int M, std::ofstream& outfile) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            outfile << matrix[i][j].real() << '\t';
        }
        outfile << std::endl;
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            outfile << matrix[i][j].imag() << '\t';
        }
        outfile << std::endl;
    }
}
int main() {
    const int N = 32; // 矩阵的大小
    const int M = 32; // 矩阵的大小
    double h = 2 * PI / N;

    // 定义二维矩阵

    vector<vector<Complex>> f(N, vector<Complex>(M));
    // 函数矩阵赋值
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            f[i][j] = sin(5 * h * i + 7 * h * j);

        }
    }
    std::ofstream outfile1("data1.txt"); //原始的数值矩阵
    printMatrix(f, N, M);
    printMatrix(f, N, M, outfile1);
    // 定义输出矩阵:频域
    vector<vector<Complex>>F = FFT2D(f);
    std::ofstream outfile2("data2.txt"); //变换后的矩阵
    std::cout << "变换前矩阵" << std::endl;
    std::cout << '\n' << '\n' << "变换后的频域信号矩阵 F(u, v):" << std::endl;
    printMatrix(F, N, M);
    printMatrix(F, N, M, outfile2);

    return 0;
}
