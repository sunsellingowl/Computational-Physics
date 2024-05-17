#define _CRT_SECURE_NO_WARNINGS //忽略安全检测
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>

#define PI 3.14159265358

//f(x,y)= sin(5x+7y)

typedef std::complex<double> Complex; //复数模板

void DFT2D(Complex** f, Complex** F, int N, int M) {
    for (int u = 0; u < N; ++u) {
        for (int v = 0; v < M; ++v) {
            F[u][v] = 0;
            for (int x = 0; x < N; ++x) {
                for (int y = 0; y < M; ++y) {
                    Complex expo = exp(Complex(0, -2 * PI * (double)(u * x) / N - 2 * PI * (double)(v * y) / M));
                    F[u][v] += f[x][y] * expo;
                }
            }
        }
    }
}

// 输出
void printMatrix(Complex** matrix, int N, int M) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            std::cout << "(" << matrix[i][j].real() << ", " << matrix[i][j].imag() << ") "<<'\t';
        }
        std::cout << std::endl;
    }
}
//重载
void printMatrix(Complex** matrix, int N, int M, std::ofstream & outfile) {
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

//内存释放
void DDDelete(Complex** f, Complex** F, int N){
	// 释放内存
	for (int i = 0; i < N; ++i) {
		delete[] f[i];
		delete[] F[i];
	}
	delete[] f;
	delete[] F;


}
int main() {
    const int N = 50;
    const int M = 50;
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

    // 计算二维离散傅立叶变换
    DFT2D(f, F, N, M);

    std::ofstream outfile1("data1.txt"); //原始的数值矩阵
    std::ofstream outfile2("data2.txt"); //变换后的矩阵


    // 输出
    std::cout << "变换前矩阵" << std::endl;
    printMatrix(f, N, M);
    printMatrix(f, N, M, outfile1);

    std::cout << '\n' << '\n' << "变换后的频域信号矩阵 F(u, v):" << std::endl;
    printMatrix(F, N, M);
    printMatrix(F, N, M, outfile2);

    DDDelete(f, F, N);//释放内存
    return 0;
}