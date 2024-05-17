#define _CRT_SECURE_NO_WARNINGS //���԰�ȫ���
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <fftw3.h>
#define PI 3.14159265358

typedef std::complex<double> Complex; //����ģ��

int main() {
    const int N = 50; // ����Ĵ�С
    const int M = 50; // ����Ĵ�С
    double h = 2 * PI / N;

    // �����ά����
    Complex** f = new Complex * [N];
    for (int i = 0; i < N; ++i) {
        f[i] = new Complex[M];
    }
    // ��������ֵ
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            f[i][j] = sin(5 * h * i + 7 * h * j);
        }
    }
    // �������:Ƶ��
    Complex** F = new Complex * [N];
    for (int i = 0; i < N; ++i) {
        F[i] = new Complex[M];
    }

    // ���� 2D FFT
    fftw_plan plan_forward = fftw_plan_dft_2d(N, N, f, F, FFTW_FORWARD, FFTW_ESTIMATE);

    // ִ�� 2D FFT
    fftw_execute(plan_forward);

    // ��� 2D FFT ���
    std::cout << "2D FFT �����" << std::endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << F[i][j].real() << " + " << F[i][j].imag() << "i\t";
        }
        std::cout << std::endl;
    }

    // ���ټƻ����ͷ��ڴ�
    fftw_destroy_plan(plan_forward);

    return 0;
}
