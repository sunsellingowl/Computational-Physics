#define _CRT_SECURE_NO_WARNINGS //���԰�ȫ���
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>

#define PI 3.14159265358

//f(x,y)= sin(5x+7y)

typedef std::complex<double> Complex; //����ģ��

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

// ���
void printMatrix(Complex** matrix, int N, int M) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            std::cout << "(" << matrix[i][j].real() << ", " << matrix[i][j].imag() << ") "<<'\t';
        }
        std::cout << std::endl;
    }
}
//����
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

//�ڴ��ͷ�
void DDDelete(Complex** f, Complex** F, int N){
	// �ͷ��ڴ�
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

    // �����ά��ɢ����Ҷ�任
    DFT2D(f, F, N, M);

    std::ofstream outfile1("data1.txt"); //ԭʼ����ֵ����
    std::ofstream outfile2("data2.txt"); //�任��ľ���


    // ���
    std::cout << "�任ǰ����" << std::endl;
    printMatrix(f, N, M);
    printMatrix(f, N, M, outfile1);

    std::cout << '\n' << '\n' << "�任���Ƶ���źž��� F(u, v):" << std::endl;
    printMatrix(F, N, M);
    printMatrix(F, N, M, outfile2);

    DDDelete(f, F, N);//�ͷ��ڴ�
    return 0;
}