#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <random>

#define PI 3.14159265358
#define N 5 // five masses 

std::vector<double> verlet(std::vector<double> Position, std::vector<double> previousPosition, double dt, double l0, double k, double m) {
	std::vector<double> acceleration(N);
	std::vector<double> newPosition(N);
	acceleration[0] = (Position[1] - Position[0] - l0) * k / m;
	acceleration[4] = (l0 - Position[4] + Position[3]) * k / m;
	for (int i = 1; i < 4; i++)
	{
		acceleration[i] = (Position[i + 1] - 2 * Position[i] + Position[i - 1]) * k / m;
	}
	for (int i = 0; i < N; i++)
	{
		newPosition[i] = 2 * Position[i] - previousPosition[i] + acceleration[i] * dt * dt;
	}
	return newPosition;
}
void printLine(std::vector<double> Position, std::ofstream & outfile) {
	for (size_t j = 0; j < Position.size(); j++) {
		outfile << Position[j] << '\t';
	}
	outfile << '\n';
}
void printLine(std::vector<double> Position) {
	for (size_t j = 0; j < Position.size(); j++) {
		std::cout << Position[j] << '\t';
	}
	std::cout << '\n';
}

int main() {
	const int n = 1000;
	const double dt = 0.01;
	const double k = 1; //弹簧系数
	const double l0 = 1; //初始长度
	const double m = 1; //重量
	std::ofstream outfile("data.txt");
	std::vector<double> Position = { 0.0, 1.0 ,2.0, 3.0, 4.0 };
	std::vector<double> newPosition(N);

	//随机初始化 Position[1]
	std::random_device rd;  // 随机设备，用于生成种子
	std::mt19937 gen(rd()); // Mersenne Twister 19937 生成器，用随机设备生成种子
	std::uniform_real_distribution<double> dis(-0.02, 0.02);	// 定义实数分布，范围是 [-0.2, 0.2)
	for (size_t i = 0; i < Position.size(); i++)
	{
		newPosition[i] = Position[i] +  dis(gen);
		//newPosition[i] = Position[i] + 0.1;
	}


	for (size_t i = 0; i < n; i++){
		printLine(Position);
		printLine(Position,outfile);
		std::vector<double> temp(N);
		temp = verlet(newPosition, Position, dt, l0, k, m);
		Position = newPosition;
		newPosition = temp;
	}
	
}
