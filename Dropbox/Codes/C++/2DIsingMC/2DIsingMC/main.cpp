//
//  2D_Ising.cpp
//  2D_Ising_MC
//
//  Created by Xianglin Liu on 1/16/19.
//  Copyright © 2019 Xianglin Liu. All rights reserved.
//

#include "Array3d.hpp"
#include <cmath>
//#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>


using std::cout;
using std::string;
using std::vector;
using Array3d::Array2d;

//	Use the Array2d class because it is easy to work with for a 2d matrix.
//	A Mersenne twister random number generator is used. Better than rand()
//	in cstdlib (commented below).

/*
 double getRandom_cstd()
 {
	int randomNumber = std::rand();
	return float(randomNumber)/float(RAND_MAX);
 }
 */

class Random_mt
{
public:
	Random_mt():
	mt(1000),dist(0.0, 1.0)
	{;}
	
	double getRandom() { return dist(mt); }
	
private:
	std::mt19937 mt;
	std::uniform_real_distribution<double> dist;
};

Random_mt random_mt;

//	a Timer class to record the calculation time.
class Timer
{
public:
	Timer(){c_start = std::clock();}
	
	double getTime()
	{
		c_end = std::clock();
		time_elapsed_s = double(c_end-c_start) / CLOCKS_PER_SEC;
		return time_elapsed_s;
	}
	
	void printTime()
	{
		double time = getTime();
		if (time < 1.)
		{
			cout << "CPU time used: " << 1000.*time << " ms\n";
		}
		else
		{
			cout << "CPU time used: " << time << " s\n";
		}
		
	}
	
	void reset_time()
	{
		c_start = std::clock();
	}
	
private:
	std::clock_t c_start;
	std::clock_t c_end;
	double time_elapsed_s;
};

// Sweep_bin is used to store the information of multiple sweeps.
struct Sweep_bin
{
public:
	vector<double> energy_ave;
	vector<double> magnetization_ave;
	vector<double> accept_ratio;
	void reset() // erase all elements and set the length to zero;
	{
		energy_ave.resize(0);
		magnetization_ave.resize(0);
		accept_ratio.resize(0);
	}
	
};

//	A Lattice class to set the lattice parameters and store the spins.
//	To facilitate the test, all methods are set public.
class Lattice
{
public:
	Lattice(int lx_in, int ly_in, bool random_in=false,
			double h_in = 0., double J_in = 1.):
	lx(lx_in), ly(ly_in), sigma(lx_in, ly_in), random(random_in),
	h(h_in), J(J_in), kb(1.0)
	{
		if (random)
			initRandom();
		else
			sigma = 1;// if random = false, set the spin of all sites to 1.
	}
	int lx;
	int ly;
	int n_step;
	bool random;
	double h;
	double J;
	const double kb;
	Array2d<int> sigma;
	
	void printSigma()
	{
		for (int j=0; j<ly; j++)
		{
			for (int i=0; i<lx; i++)
			{
				printf( "%2d",sigma(i,j));
			}
			cout << endl;
		}
	}
	
	void thermalize(int n, double T, Sweep_bin& sweep_bin)
	{
		
		sweep_bin.reset();
		for (int i=0; i<n; i++)
		{
			updateSweep_Metro(T, sweep_bin);
			calObservables(sweep_bin);
		}
	}

    void thermalize(int n, double T, Sweep_bin& sweep_bin, std::ofstream& file)
    {

        sweep_bin.reset();
        for (int i=0; i<n; i++)
        {
            updateSweep_Metro(T, sweep_bin);
            calObservables(sweep_bin);
            for (int j=0; j<ly; j++)
            {
                for (int i=0; i<lx; i++)
                {
                    file << sigma(i,j) << " ";
                }
            }
            file << endl;
        }
    }
	
	// spins initialized randomly
	void initRandom()
	{
		for (int j=0; j<ly; j++)
			for (int i=0; i<lx; i++)
			{
				if (random_mt.getRandom() < 0.5)
				{
					sigma(i,j) = 1;
				}
				else
				{
					sigma(i,j) = -1;
				}
			}
	}
	
	//	calculate energy corresponds to site (i,j).
	double Hamiltonian(int i, int j)
	{
		double H = 0;
		H = -sigma(i,j)*(h + J * ( sigma((i+1)%lx,j) + sigma(i,(j+1)%ly ) ) );
		return H;
	}
	
	// calculate the change of energy due to flip of the spin at site (i,j)
	double changeE(int i, int j)
	{
		double deltaE = 0;
		double change = sigma((i+1)%lx,j) + sigma((i-1)%lx,j)
		+ sigma(i,(j+1)%ly) + sigma(i,(j-1)%ly);
		deltaE = 2*sigma(i,j)*(h + J*change);
		return deltaE;
	}
	
	//	updata one MC sweep with the Heat-Bath algorithm.
	void updateSweep_HB(double T, Sweep_bin &sweep_bin)
	{
		double beta = 1./(kb*T);
		double accept_tot = 0;
		for (int j=0; j<ly; j++)
		{
			for (int i=0; i<lx; i++)
			{
				double rd = random_mt.getRandom();
				double deltaE = changeE(i, j); // change of E due to flip of the spin.
				if (1./(1.+exp(beta*deltaE))>rd)
				{
					sigma(i,j) = -sigma(i,j);
					accept_tot += 1.0;
				}
			}
		}
		sweep_bin.accept_ratio.push_back( accept_tot/(lx*ly) );
	}
	
	//	updata one MC sweep with the Metropolis alogrithm.
	void updateSweep_Metro(double T, Sweep_bin &sweep_bin)
	{
		double beta = 1./(kb*T);
		double accept_tot = 0;
		for (int j=0; j<ly; j++)
		{
			for (int i=0; i<lx; i++)
			{
				double rd = random_mt.getRandom();
				double deltaE = changeE(i, j); // change of E due to flip of the spin.
				if (deltaE > 0)
				{
					if (exp(-beta*deltaE) > rd)
					{
						sigma(i,j) = -sigma(i,j);
						accept_tot += 1.0;
					}
				}
				else
				{
					sigma(i,j) = -sigma(i,j);
					accept_tot += 1.0;
				}
			}
		}
		sweep_bin.accept_ratio.push_back( accept_tot/(lx*ly) );
	}
	
	//	calculate the observables in the Sweep_bin class
	void calObservables(Sweep_bin &sweep_bin)
	{
		double en_tot = 0;
		double mag_tot = 0; // double type magnetization.
		for (int j=0; j<ly; j++)
		{
			for (int i=0; i<lx; i++)
			{
				mag_tot += sigma(i,j);
				en_tot += Hamiltonian(i, j); // record the energy before updating
			}
		}
		sweep_bin.energy_ave.push_back( en_tot/(lx*ly) );
		sweep_bin.magnetization_ave.push_back( mag_tot/(lx*ly) );
	}
	
};

void printSweep_bin(const Sweep_bin& sweep_bin)
{
	printf("     step     energy     magnetization     accept_ratio\n");
	for (int i=0; i<sweep_bin.energy_ave.size(); i++)
	{
		printf("%8d %14.7f %14.7f %14.7f\n", i+1,
			   sweep_bin.energy_ave[i], sweep_bin.magnetization_ave[i],
			   sweep_bin.accept_ratio[i]);
	}
}

double mean( const vector<double>& list )
{
	double sum = 0;
	for (int i=0; i<list.size(); i++){ sum += list[i]; }
	return sum/double(list.size());
}

double variance(const vector<double>& list )
{
	double var = 0;
	double average = mean(list);
	for (int i=0; i<list.size(); i++)
	{
		var += (list[i]-average)*(list[i]-average);
	}
	var = var/double(list.size());
	return var;
}

//==========================================================
//	these functions are only for test purpose only.
void TEST_PRINT_LINE(string title)
{
	printf("============================================\n");
	cout << title << endl;
	printf("============================================\n");
}

//
void TEST_getRandom()
{
	TEST_PRINT_LINE("TEST_getRrandom");
	for (int i=0; i<12; i++)
	{
		printf("%.8f  ", random_mt.getRandom());
		if ( i % 4 == 3)
			std::cout << "\n";
	}
}

void TEST()
{
	Timer time;
	TEST_PRINT_LINE("Test begins");
	time.printTime();
	time.reset_time();
	TEST_getRandom();
	Lattice lattice_test(30, 30, false, 0.0, 1.0);
	TEST_PRINT_LINE("Initial state of spins on lattice");
	lattice_test.printSigma();
	Sweep_bin sweep_bin_test;
	lattice_test.thermalize(100, 1., sweep_bin_test);
	TEST_PRINT_LINE("After thermalize");
	lattice_test.printSigma();
}

//#define TEST
#ifdef TEST
	TEST();
#endif

//==========================================================
// This is the main function
//==========================================================

int main(int argc, const char * argv[])
{
	Timer time; //record time
	//	set up parameters
	double temp_ini = 4.00; // initial temperature (T)
	double temp_end = 0.01; // endding T
	int n_temp = 20;   // number of T
	int n_warm = 2000; // number of MC steps discarded for the first T
	int n_skip = 2000; // number of MC steps discarded when T changes
	int n_measure = 4000; // number of MC steps to measure the observables
    int n_output = 100; // number of MC samples to be written at each T
	
	//	calculate temperatures
	vector<double> temp;
	for (int i=0; i<n_temp; i++)
	{
		double delta_temp = 0;
		if (n_temp > 1)
		{
			delta_temp= (temp_end-temp_ini)/(n_temp-1);
		}
		
		temp.push_back( temp_ini+i*delta_temp );
	}
	
	//	set up lattice
	int lx = 20;
	int ly = 20; // lattice size (lx*ly)
	double h = 0.00; // magnetic field
	double J = 1.0;  // coupling constant
	bool randomInit = true; // random initialization or not for the spin configuration
	Lattice lattice(lx, ly, randomInit, h, J);
	
	//	print the parameters
	printf("lx = %3d\n", lx);
	printf("ly = %3d\n", ly);
	printf("magnetic field h = %12.8f\n", h);
	printf("exchange interaction J = %12.8f\n", J);
	printf("n_warm = %8d\n", n_warm);
	printf("n_skip = %8d\n", n_skip);
	printf("n_measure = %8d\n", n_measure);
	
	//	storage for the ouput results
	vector<double> energy_T;
	vector<double> magnetization_T;
	vector<double> specificHeat_T;
	vector<double> magneticSusceptibility_T;
	
	//	MC simulation
	Sweep_bin sweep_bin;
	lattice.thermalize(n_warm, temp_ini, sweep_bin); // initial warm up
	
    for (int i = 0; i < n_temp; i++ )
	{
		double T = temp[i];
		double kb = lattice.kb;
		double beta  = 1./(kb*T);
		lattice.thermalize(n_skip, T, sweep_bin);	// discard the nonequilibrium MC steps
        
        std::ofstream file;
        string file_name = "spins_";
        file_name += std::to_string(i);
        file.open(file_name);
		lattice.thermalize(n_measure, T, sweep_bin);  // make measurements
		lattice.thermalize(n_output, T, sweep_bin, file);  // write the spin configurations
        file.close();
		
		int NN = lx*ly;
		double en_ave = mean(sweep_bin.energy_ave);
		double mag_ave = mean(sweep_bin.magnetization_ave);
		// times NN to accout for the extra NN divided when (en/NN)**2 is applied
		// should have been en**2/NN
		double cv = NN*variance(sweep_bin.energy_ave)/(kb*T*T);
		double chi = NN*beta*variance(sweep_bin.magnetization_ave);

		energy_T.push_back( en_ave );
		magnetization_T.push_back( mag_ave );
		specificHeat_T.push_back( cv );
		magneticSusceptibility_T.push_back( chi );
	}
	
	// output of the results
	std::ofstream file;
	file.open("output.txt");
	file << "temperature " << "energy " << "magnetization "
	<< "specific_heat " << "magnetic_susceptibility " << endl;
	for (int i=0; i < n_temp; i++ )
	{
		file << std::setw(12) << temp[i] << std::setw(12) <<
		energy_T[i] << std::setw(12) << magnetization_T[i] <<
		std::setw(16) << specificHeat_T[i] << std::setw(16) <<
		magneticSusceptibility_T[i] << "\n";
	}
	file.close();
	
	//	print the results on screen
	time.printTime();
	cout << "temperature " << "energy " << "magnetization "
	<< "specific_heat " << "magnetic_susceptibility " << endl;
	for (int i=0; i < n_temp; i++ )
	{
		cout << std::setw(16) << temp[i] << std::setw(16) <<
		energy_T[i] << std::setw(16) << magnetization_T[i] <<
		std::setw(16) << specificHeat_T[i] << std::setw(16) <<
		magneticSusceptibility_T[i] << "\n";
	}
	
	return 0;
}

