#define _CRT_SECURE_NO_WARNINGS //忽略安全检测
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <fftw3.h>
#define PI 3.14159265358

typedef std::complex<double> Complex; //复数模板

int main() {
    const int N = 50; // 矩阵的大小
    const int M = 50; // 矩阵的大小
    double h = 2 * PI / N;

    // 定义二维矩阵
    Complex** f = new Complex * [N];
    for (int i = 0; i < N; ++i) {
        f[i] = new Complex[M];
    }
    // 函数矩阵赋值
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            f[i][j] = sin(5 * h * i + 7 * h * j);
        }
    }
    // 输出矩阵:频域
    Complex** F = new Complex * [N];
    for (int i = 0; i < N; ++i) {
        F[i] = new Complex[M];
    }

    // 创建 2D FFT
    fftw_plan plan_forward = fftw_plan_dft_2d(N, N, f, F, FFTW_FORWARD, FFTW_ESTIMATE);

    // 执行 2D FFT
    fftw_execute(plan_forward);

    // 输出 2D FFT 结果
    std::cout << "2D FFT 结果：" << std::endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << F[i][j].real() << " + " << F[i][j].imag() << "i\t";
        }
        std::cout << std::endl;
    }

    // 销毁计划和释放内存
    fftw_destroy_plan(plan_forward);

    return 0;
}
