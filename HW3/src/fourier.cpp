/*!
 * @file fourier.cpp
 * @ingroup libfourier
 * @brief evaluate sum and derivative of sine and cosine function
 * @author pistack (Junho Lee)
 * @date 2021. 10. 24.
 */

#include <algorithm>
#include <cmath>
#include "fourier.hpp"

using namespace std;

fourier& fourier::operator=(const fourier &copy)
{
	f_num_fourier = copy.f_num_fourier;
	f_period = copy.f_period;
	f_c = copy.f_c;
	return *this;
}

void fourier::update(vector<double> c){ f_c = c;}

double fourier::eval(double t)
{
	int term = 2*f_num_fourier;
	double omega = 2*pi/f_period;
	double tmp;
	double y=0;
	
	for(int j=0; j<term; j += 2)
	{
		tmp = (j/2+1)*omega;
		if(f_c[j] != 0.0)
		y += f_c[j]*sin(tmp*t);
		if(f_c[j+1] != 0.0)
		y += f_c[j+1]*cos(tmp*t);
	}

	return y;
}

vector<double> fourier::eval(vector<double> t)
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

double fourier::deriv(double t)
{
	int term = 2*f_num_fourier;
	double omega = 2*pi/f_period;
	double tmp;
	double yp=0;
	
	for(int j=0; j<term; j += 2)
	{
		tmp = (j/2+1)*omega;
		if(f_c[j] != 0.0)
		yp += f_c[j]*tmp*cos(tmp*t);
		if(f_c[j+1] != 0.0)
		yp -= f_c[j+1]*tmp*sin(tmp*t);
	}

	return yp;
}

vector<double> fourier::deriv(vector<double> t)
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

double fourier::nderiv(int n, double t)
{
	int term = 2*f_num_fourier;
	double omega = 2*pi/f_period;
	double tmp;
	vector<double> c = f_c;
	double yp=0;

	if(n % 2 == 0)
	{}
	else
	{
		for(int j=0; j<term; j+=2)
		{
			iter_swap(c.begin()+j, c.begin()+j+1);
			c[j] *= -1.0;
		}
	}
	for(int j=0; j<term; j += 2)
	{
		tmp = pow((j/2+1)*omega, double(n));
		if(c[j] != 0.0)
		yp += c[j]*tmp*sin(tmp*t);
		if(f_c[j+1] != 0.0)
		yp += f_c[j+1]*tmp*cos(tmp*t);
	}
	if((n%2 == 0 && n/2 % 2 != 0) ||
	(n%2 != 0 && (n-1)/2 % 2 != 0))
	yp *= -1.0;
	return yp;
}

vector<double> fourier::nderiv(int n, vector<double> t)
{
	int term = 2*f_num_fourier;
	int n_t = t.size();
	double omega = 2*pi/f_period;
	double tmp;
	vector<double> c = f_c;
	vector<double> yp(n, 0);

	if(n % 2 == 0)
	{
		if(n/2 %2 == 0){}
		else
		{
			transform(c.begin(), c.end(), c.begin(), 
			[](double &x){return x *= -1.0;});
		}
	}
	else
	{
		for(int j=0; j<term; j+=2)
		{
			iter_swap(c.begin()+j, c.begin()+j+1);
			c[j] *= -1.0;
		}
		if((n-1)/2 %2 == 0){}
		else
		{
			transform(c.begin(), c.end(), c.begin(), 
			[](double &x){return x *= -1.0;});
		}
	}
	
	for(int i=0; i<n_t; i++)
	{
		for(int j=0; j<term; j += 2)
		{
			tmp = pow((j/2+1)*omega, double(n));
			if(c[j] != 0.0)
			yp[i] += c[j]*tmp*sin(tmp*t[i]);
			if(f_c[j+1] != 0.0)
			yp[i] += f_c[j+1]*tmp*cos(tmp*t[i]);
		}
	}

	return yp;
}


