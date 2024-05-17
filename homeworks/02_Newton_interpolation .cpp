#define _CRT_SECURE_NO_WARNINGS //���԰�ȫ���
#include <iostream>
#include <cmath>
#include <fstream>  // �����ļ���ͷ�ļ�

#define PI 3.14159265358

//ϵ������
void caculate_coefficients_a(double a[], double x[], double y[], int size);
//���
double interpolation(double a[], double x[], double x_inter, int length);


int main(void)
{
	const int N = 10;
	double x[N];
	double y[N];

	//����cosֵ
	for (int i = 0; i < N; i++)
	{
		x[i] = PI / (N - 1) * i;
		y[i] = cos( x[i] );
	}

	double a[N]; //ϵ��
	caculate_coefficients_a(a, x, y, N);

	//���100����
	const int NN = 100;
	double x_inter[NN];
	double f[NN];
	for (int i = 0; i < NN; i++)
	{
		x_inter[i] = PI / (NN - 1) * i;
		f[i] = interpolation(a, x, x_inter[i], N);
	}

	std::ofstream outfile1;
	outfile1.open("data_ori.dat");
	outfile1 << "ԭʼx" << '\t' << "ԭʼy" << std::endl;
	for (int i = 0; i < N; i++)
	{
		outfile1 << x[i] << '\t' << y[i] << std::endl;
	}
	outfile1.close();

	std::ofstream outfile2;
	outfile2.open("data_interpolated.dat");
	outfile2 << "���x" << '\t' << "���y" << std::endl;
	for (int i = 0; i < NN; i++)
	{
		outfile2 << x_inter[i] << '\t' << f[i] << std::endl;
	}
	outfile2.close();
	return 0;
}

//ϵ������
void caculate_coefficients_a(double a[], double x[], double y[], int size)
{
	int n = 0;
	a[0] = y[0]; //a0 ��������
	//a1 �����
	for (int n = 1; n < size; n++)
	{
		double sum = 0;
		for (int i = 0; i <= n; i++)
		{	
			double multi = 1;
			for (int j = 0; j <= n ; j++) 
			{
				if (i != j)
				{
					multi *= (x[i] - x[j]);
				}
			}
			sum += y[i] / multi;
		}
		a[n] = sum;
	}

}
//���
double interpolation(double a[], double x[], double x_inter, int length_of_poly)
{
	double f = 0;
	for (int j = 0; j < length_of_poly; j++)
	{
		double multiplicate = 1;
		for (int k = 0; k < j; k++)
		{
			multiplicate *= (x_inter - x[k]);
		}
		f += a[j] * multiplicate;
	}
	return f;
}


