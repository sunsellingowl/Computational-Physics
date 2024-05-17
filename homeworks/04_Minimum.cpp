#define _CRT_SECURE_NO_WARNINGS //ºöÂÔ°²È«¼ì²â
#include <iostream>
#include <cmath>
#include <fstream> 

#define PI 3.14159265358

//g(x,y)=sin(x+y)+cos(x+2*y) FUNCTION


double deepest_decent(double &x, double  &y, const double epsilon, double alpha = 0.1);

int main()
{
	const double epsilon = 0.001; //arrcuracy
	double x0 = 0.0;
	double y0 = 0.0; //initial point
	double x = x0, y = y0;

	std::cout << "initial point is:" << "(" << x0 << "," << y0 << ")" << '\n'
		<< "minimum point is:" << "(" << x << "," << y << ")" << '\n'
		<< "the minimum is:" << deepest_decent(x, y, epsilon) << '\n' << '\n' << '\n';
	system("pause");
	return 0;
}

double deepest_decent(double &x, double &y, const double epsilon, double alpha)
{
	for (int i = 0; i < 100; i++)
	{
		double partial_x, partial_y, grad;
		partial_x = cos(x + y) - sin(x + 2 * y);
		partial_y = cos(x + y) - 2 * sin(x + 2 * y);
		grad = sqrt(partial_x * partial_x + partial_y * partial_y);
		if (grad < epsilon)
			break;
		x = x - alpha * partial_x / abs(grad);
		y = y - alpha * partial_y / abs(grad);
	}
	return sin(x + y) + cos(x + 2 * y);
}
