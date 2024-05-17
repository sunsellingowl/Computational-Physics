#define _CRT_SECURE_NO_WARNINGS //���԰�ȫ���
#include<iostream>
#include<cmath>
#include <fstream>  // �����ļ���ͷ�ļ�
#include <cstdlib>  // ���ڷ��ʻ�������

#define PI 3.1415926

int main(void)
{
	const int N = 500;
	const double dt = 4 * PI / N;
	double x[N]; //positon
	double v[N]; //velocity

	x[0] = 0;
	v[0] = 1;

    std::ofstream outfile;  // ��������ļ�������
    // ����·��
    std::string filePath = "damped_oscillation.dat";

    outfile.open(filePath);  // ����Ϊ "damped_oscillation.dat" ���ļ�

    // F = -k*x - b*v
    outfile << "time" << "\t" << "x" << "\t" << 'v' << std::endl;
    for (int i = 1; i < N; i++)
    {
        x[i] = x[i - 1] + v[i - 1] * dt;
        v[i] = v[i - 1] - x[i - 1] * dt - 0.5 * v[i - 1] * dt;	// k/m=1, b/m=0.5
        
        // д�����ݵ��ļ�
        outfile << (i - 1) * dt << "\t" << x[i - 1] << "\t" << v[i - 1] << std::endl;

    }
    // �ر��ļ�
    outfile.close();

	return 0;
}