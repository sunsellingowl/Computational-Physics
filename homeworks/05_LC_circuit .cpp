#define _CRT_SECURE_NO_WARNINGS //忽略安全检测
#include <iostream>
#include <cmath>
#include <fstream> 

#define PI 3.14159265358

void RungeKutta(double I[], 
				double q[], 
				double E[],
				double R = 1.5, //电阻
				double C = 2, //电容
				double L = 3,  //电感
				double t = 0.01,  //时间间隔
				double A = 3);//振幅

int main()
{
	double E[5000] = {};  //激励电压（这里考虑使用正弦波）
	double q[5000] = {};  //电容电荷
	double I[5000] = {};  //电路电流   

	I[0] = 0;
	q[0] = 0;
	RungeKutta(I, q, E);
	return 0;
}


void RungeKutta(double I[], double q[], double E[], double R, double C, double L, double t, double A)
{
	double E1, E2, E3, E4;  //电压1至4阶导数
	double I1, I2, I3, I4;  //电流1至4阶导数
	double t2 = t * t / 2;
	double t3 = t * t2 / 3;
	double t4 = t * t3 / 4;
	std::ofstream file_out("RungeKutta.dat");
	file_out << "时间" << '\t' << "电流" << '\t' << "电感电压" << '\t' << "电阻电压" << '\t' << "电容电压" << '\t' << "电压" << std::endl;
	for (int i = 0; i < 5000 - 1; i++)
	{
		E[i] = A * sin(i * t);
		E1 = A * cos(i * t);
		E2 = -A * sin(i * t);
		E3 = -A * cos(i * t);
		E4 = A * sin(i * t);
		I1 = E[i] / L - I[i] * R / L - q[i] / C / L;
		I2 = E1 / L - R * I1 / L - I[i] / L / C;
		I3 = E2 / L - I2 * R / L - I1 / C / L;
		I4 = E3 / L - I3 * R / L - I2 / C / L;
		q[i + 1] += q[i] + I[i] * t + I1 * t2 + I2 * t3 + I3 * t4;
		I[i + 1] += I[i] + I1 * t + I2 * t2 + I3 * t3 + I4 * t4;
		
		file_out << i * t << '\t' << I[i] << '\t' << L * I1 << '\t' << I[i] * R << '\t' << q[i] / C << '\t' << E[i] << std::endl;
	}
	file_out.close();
}
