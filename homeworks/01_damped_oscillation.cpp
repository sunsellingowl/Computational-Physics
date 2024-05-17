#define _CRT_SECURE_NO_WARNINGS //忽略安全检测
#include<iostream>
#include<cmath>
#include <fstream>  // 包含文件流头文件
#include <cstdlib>  // 用于访问环境变量

#define PI 3.1415926

int main(void)
{
	const int N = 500;
	const double dt = 4 * PI / N;
	double x[N]; //positon
	double v[N]; //velocity

	x[0] = 0;
	v[0] = 1;

    std::ofstream outfile;  // 创建输出文件流对象
    // 完整路径
    std::string filePath = "damped_oscillation.dat";

    outfile.open(filePath);  // 打开名为 "damped_oscillation.dat" 的文件

    // F = -k*x - b*v
    outfile << "time" << "\t" << "x" << "\t" << 'v' << std::endl;
    for (int i = 1; i < N; i++)
    {
        x[i] = x[i - 1] + v[i - 1] * dt;
        v[i] = v[i - 1] - x[i - 1] * dt - 0.5 * v[i - 1] * dt;	// k/m=1, b/m=0.5
        
        // 写入数据到文件
        outfile << (i - 1) * dt << "\t" << x[i - 1] << "\t" << v[i - 1] << std::endl;

    }
    // 关闭文件
    outfile.close();

	return 0;
}