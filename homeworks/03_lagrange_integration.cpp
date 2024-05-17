#define _CRT_SECURE_NO_WARNINGS //忽略安全检测
#include <iostream>
#include <cmath>
#include <fstream> 

#define PI 3.14159265358

double TrapezoidalIntegrate(const int n);//Trapezoidal rules 
double LagrangeIntegrate(const int n);//Lagrange rules

int main(void)
{
	std::ofstream outfile;
	outfile.open("data.dat");

	outfile << "间距个数" << '\t' << "TInte" << '\t' << "LInte" << '\t' << "T2one" << '\t' << "L2one" << std::endl;
	for (int n = 2; n < 30; n++)
	{
		double Trap = TrapezoidalIntegrate(n);
		double Lag = LagrangeIntegrate(n);
		outfile << n<<'\t' << Trap << '\t' << Lag << '\t' << fabs(1 - Trap) << '\t' << fabs(1 - Lag) <<std::endl ;
	}
	outfile.close();

	return 0;
}

double TrapezoidalIntegrate(const int n)
{
	const double h = PI / 2 / n;

	double sum = 0;

	for (int i = 0; i < n; i++)
	{
		sum += sin(i * h) + sin(i * h + h);	
	}

	return sum * h / 2;
}

double LagrangeIntegrate(const int n)
{
	const double h = PI / 2 / n;
	double sum = 0;
	if (n % 2 == 0) //奇数个点，间距偶数个
	{
		for (int i = 0; i <= n / 2 - 1; i++)
		{
			sum += sin(2 * i * h) + 4 * sin(2 * i * h + h) + sin(2 * i * h + 2 * h);
		}

		return sum * h / 3;
	}
	else
	{
		for (long i = 0; i <= (n - 1) / 2 - 1; i++)
		{
			sum += sin(2 * i * h) + 4 * sin(2 * i * h + h) + sin(2 * i * h + 2 * h);
		}
		sum *= h / 3;
		return sum + (sin(n * h - h) * h + sin(n * h) * h) / 2; //最后一个线性直接算
	}
}