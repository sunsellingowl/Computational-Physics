#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>

#define PI 3.14159265358
using namespace std;

typedef complex<double> Complex;
// һάFFT
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
// ��άFFT
vector<vector<Complex>> FFT2D(vector<vector<Complex>>& f) {
    int row = f.size();    // ����
    int col = f[0].size(); // ����

    // ��ÿһ�н���һάFFT
    for (int i = 0; i < row; i++) {
        f[i] = FFT(f[i]);
    }

    // ��ÿһ�н���һάFFT
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

// ���
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
    const int N = 32; // ����Ĵ�С
    const int M = 32; // ����Ĵ�С
    double h = 2 * PI / N;

    // �����ά����

    vector<vector<Complex>> f(N, vector<Complex>(M));
    // ��������ֵ
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            f[i][j] = sin(5 * h * i + 7 * h * j);

        }
    }
    std::ofstream outfile1("data1.txt"); //ԭʼ����ֵ����
    printMatrix(f, N, M);
    printMatrix(f, N, M, outfile1);
    // �����������:Ƶ��
    vector<vector<Complex>>F = FFT2D(f);
    std::ofstream outfile2("data2.txt"); //�任��ľ���
    std::cout << "�任ǰ����" << std::endl;
    std::cout << '\n' << '\n' << "�任���Ƶ���źž��� F(u, v):" << std::endl;
    printMatrix(F, N, M);
    printMatrix(F, N, M, outfile2);

    return 0;
}
