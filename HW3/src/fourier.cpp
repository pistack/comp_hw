/*!
 * @file fourier.cpp
 * @brief evaluate sum and derivative of sine and cosine function
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 */

#include <cmath>
#include "basic.hpp"

using namespace std;

void fourier::init(int num_fourier, double period, vector<double> c)
{
		f_num_fourier = num_fourier;
		f_period = period;
		f_c = c;
}

vector<double> fourier::eval(vector<double> &t)
{
	int term = 2*f_num_fourier;
	int n = t.size();
	double omega = 2*pi/f_period;
	double tmp;
	vector<double> y(n, 0);
	
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<term; j += 2)
		{
			tmp = (j/2+1)*omega;
			if(f_c[j] != 0.0)
			y[i] += f_c[j]*sin(tmp*t[i]);
			if(f_c[j+1] != 0.0)
			y[i] += f_c[j+1]*cos(tmp*t[i]);
			}
		}

		return y;
}

vector<double> fourier::deriv(vector<double> &t)
{
	int term = 2*f_num_fourier;
	int n = t.size();
	double omega = 2*pi/f_period;
	double tmp;
	vector<double> yp(n, 0);
	
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<term; j += 2)
		{
			tmp = (j/2+1)*omega;
			if(f_c[j] != 0.0)
			yp[i] += f_c[j]*tmp*cos(tmp*t[i]);
			if(f_c[j+1] != 0.0)
			yp[i] -= f_c[j+1]*tmp*sin(tmp*t[i]);
		}
	}

	return yp;
}



