#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <ctime>

using namespace std;

#define M_PI 4.1415926

const double PI2 = 2 * M_PI;		// 2Pi

const int L_size = 40;					// Lattice Size
const int MC_Steps = 100000;			// Total Monte Carlo Steps
const int MCM_Steps = 80000;				// Monte Carlo Measure Steps
const double J = -1;				// J=1:AFM; J=-1:FM
double T = 1;					    // Monte Carlo temperature
double H = 0;				     	// Magnetic field


double Magnetization(const int* sp)
{
	int m = 0;
	for (int i = 0; i < L_size; i++)
	{
		m += sp[i];
	}
	return 1.0 * m / L_size;
}

double ExchangeEnergy(const int* sp)
{
	double e = 0;
	for (int i = 0; i < L_size; i++)
	{
		e += sp[i] * sp[(i + 1) % L_size];
	}
	return J * e / L_size;
}


void Save(const int sp[])
{
	ofstream out("spin.dat");
	for (int i = 0; i < L_size; i++)
	{
		out << sp[i] << "\t";
	}
	out << endl;
	out.close();
}

void MonteCarlo(int sp[])
{
	int iaccept = 0;			//accept percent
	double sum_Energy = 0;				//sum of energy
	double sum_Energy2 = 0;				//square of sum_Energy
	double sum_Mag = 0;				//total magnitization 

	for (int mc = MC_Steps; mc > 0; mc--)
	{
		for (int j = 0; j < L_size; j++)
		{
			const int i = rand() % L_size;

			const double de = 2 * sp[i] * (J * (sp[(i - 1 + L_size) % L_size] + sp[(i + 1) % L_size]) - H);
			//delta_E=H_old-H_new

			if ((de > 0) || (RAND_MAX * exp(de / T) > rand()))  // judge whether accept the change of the spin .
			{
				sp[i] *= -1; 	 // spin si update
				iaccept++;
			}
		}

		if (mc <= MCM_Steps)     			 //measure
		{
			const double m = Magnetization(sp);
			sum_Mag += m;
			const double e = ExchangeEnergy(sp) - m * H;

			sum_Energy += e;
			sum_Energy2 += e * e;
		}
	}

	sum_Energy /= MCM_Steps;
	sum_Energy2 /= MCM_Steps;
	sum_Mag /= MCM_Steps;
	const double capacity = L_size * (sum_Energy2 - sum_Energy * sum_Energy) / T / T;

	ofstream out("data.dat");
	out << T << '\t' << H << '\t' << sum_Energy << '\t' << capacity << '\t' << sum_Mag << '\t' << 1.0 * iaccept / MC_Steps / L_size << endl;
	out.close();
}

int main()
{
	srand(time(NULL));

	int spin[L_size];

	for (int i = 0; i < L_size; i++)
		spin[i] = 2 * (rand() % 2) - 1;		//initialize spin randomly   

	H = 0;
	for (T = 2.5; T > 0.1; T -= 0.05)
	{
		cout << "Temperature=" << T << endl;
		MonteCarlo(spin);
		Save(spin);
	}

	return 0;
}