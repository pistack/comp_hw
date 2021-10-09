/*
HW3
C++ code for homework3 in computer1 class in Yonsei University
Using least action principle with Markov Chain Monte Carlo Method
*/

#include <cmath>
#include <random>
#include <tuple>
#include <vector>
#include <iostream>

using namespace std;

const double pi = 3.141592653589793;

default_random_engine generator;
uniform_real_distribution<double> distribution(-1, 1);

vector<double>
scale_and_add_vector(vector<double> &v, double scale, double add)
{
  int n = v.size();
  vector<double> avpb(n, 0);

  for(int i=0; i<n; i++)
    {
      avpb[i] = scale*v[i]+add;
    }

  return avpb;
}

tuple<vector<double>, vector<double>>
sum_of_sine(vector<double> &t, vector<double> c, int num_sine)
{
  int n = t.size();
  double period = 4*(t[n-1]-t[0]);
  double omega = 2*pi/period;
  double tmp = 0;
  vector<double> y(n, 0);
  vector<double> deriv(n, 0);

  for(int i = 0; i<n; i++)
    {
      for(int j = 0; j<num_sine; j++)
	{
	  if(c[j] != 0.0)
	    {
	      tmp = (j+1)*omega;
	      y[i] += c[j]*sin(tmp*t[i]);
	      deriv[i] += c[j]*tmp*cos(tmp*t[i]);
	    }
	}
    }

  return make_tuple(y, deriv);
}

vector<double> move_step(vector<double> &init_guess, double step)
{
  int n = init_guess.size();
  vector<double> move(n, 0);
  for(int i=0; i<n; i++)
    {
      move[i] = init_guess[i] + step*distribution(generator);
      if(move[i] > 1)
	move[i] = 1;
      if(move[i] < -1)
	move[i] = -1;
    }
  return move;
}

double eval_action(vector<double> &t, vector<double> zeta,
		   vector<double> deriv_zeta, vector<double> theta,
		   vector<double> deriv_theta)
{
  int n = t.size();
  double action = 0;
  for(int i=1; i<n; i++)
    {
      double dt = t[i]-t[i-1];
      double zeta_inv = 1/abs(zeta[i]);
      action += dt*(0.5*(pow(deriv_zeta[i],2.0)+
			 pow((zeta[i]*deriv_theta[i]),2.0))+
		    zeta_inv);
    }
  return action;
}


tuple<double,vector<double>,vector<double>,vector<double>>
HW3(double zeta_min, double t0, int n, int num_sine, int num_iter, double step)
{
  // int n_accept = 0;
  double zeta_max = zeta_min/(2*zeta_min-1);
  double a = 0.5*(zeta_min+zeta_max);
  double tmax = pi*pow(a,1.5);
  double scale_zeta;
  double scale_theta;
  double min_action;
  vector<double> t(n, 0);
  vector<double> min_c_zeta(num_sine, 0);
  vector<double> min_c_theta(num_sine, 0);
  vector<double> tmp_c_zeta(num_sine, 0);
  vector<double> tmp_c_theta(num_sine, 0);
  vector<double> min_zeta(n, 0);
  vector<double> min_theta(n, 0);
  vector<double> min_deriv_zeta(n, 0);
  vector<double> min_deriv_theta(n, 0);

  for(int i = 0; i<n; i++)
    {
      t[i] = double(i)/double(n-1)*tmax;
    }

  min_c_zeta[0] = 1;
  min_c_theta[0] = 1;

  tie(min_zeta, min_deriv_zeta) = sum_of_sine(t, min_c_zeta, num_sine);
  tie(min_theta, min_deriv_theta) = sum_of_sine(t, min_c_theta, num_sine);
  scale_zeta = (zeta_max-zeta_min)/min_zeta[n-1];
  scale_theta = pi/min_theta[n-1];
  min_zeta = scale_and_add_vector(min_zeta, scale_zeta, zeta_min);
  min_theta = scale_and_add_vector(min_theta, scale_theta, 0);
  min_deriv_zeta = scale_and_add_vector(min_deriv_zeta, scale_zeta, 0);
  min_deriv_theta = scale_and_add_vector(min_deriv_theta, scale_theta, 0);
  min_action = eval_action(t, min_zeta, min_deriv_zeta,
			   min_theta, min_deriv_theta);

  for(int i = 0; i < num_iter; i++)
    {
      double scale_zeta, scale_theta;
      double tmp_action;
      vector<double> tmp_zeta(n, 0);
      vector<double> tmp_theta(n, 0);
      vector<double> tmp_deriv_zeta(n, 0);
      vector<double> tmp_deriv_theta(n, 0);
      do
	{
	  tmp_c_zeta = move_step(min_c_zeta, step);
          tmp_c_theta = move_step(min_c_theta, step);
          tie(tmp_zeta, tmp_deriv_zeta) = \
	    sum_of_sine(t, tmp_c_zeta, num_sine);
	  tie(tmp_theta, tmp_deriv_theta) = \
	    sum_of_sine(t, tmp_c_theta, num_sine);
	}
      while(abs(tmp_zeta[n-1]) < 1e-8 || abs(tmp_theta[n-1]) < 1e-8);

      scale_zeta = (zeta_max-zeta_min)/tmp_zeta[n-1];
      scale_theta = pi/tmp_theta[n-1];
      tmp_zeta = scale_and_add_vector(tmp_zeta, scale_zeta, zeta_min);
      tmp_deriv_zeta = scale_and_add_vector(tmp_deriv_zeta, scale_zeta, 0);
      tmp_theta = scale_and_add_vector(tmp_theta, scale_theta, 0);
      tmp_deriv_theta = scale_and_add_vector(tmp_deriv_theta, scale_theta, 0);
      tmp_action = eval_action(t, tmp_zeta, tmp_deriv_zeta,
			       tmp_theta, tmp_deriv_theta);

      if(tmp_action < min_action)
	{
	  min_action = tmp_action;
	  min_c_zeta = tmp_c_zeta;
	  min_c_theta = tmp_c_theta;
	  min_zeta = tmp_zeta;
	  min_theta = tmp_theta;
	  // n_accept++;
	}
    }

  // cout << double(n_accept)/double(num_iter)*100 << endl;

  return make_tuple(min_action,
		    scale_and_add_vector(t,1, t0),
		    min_zeta,
		    min_theta);
}

int main(void)
{
  int n = 100;
  int num_sine = 4;
  int num_iter = pow(2,24);
  double zeta_min = 0.9;
  double t0 = 0;
  double step = 0.25;
  double min_action;
  vector<double> t(n,0);
  vector<double> zeta(n,0);
  vector<double> theta(n,0);

  tie(min_action, t, zeta, theta) = HW3(zeta_min,
					t0, n, num_sine, num_iter, step);
  cout << min_action << endl;
  /*
  for(int i=0; i<n; i++)
    {
      cout << t[i] << '\t' << zeta[i] << endl;
    }
  */
  
  return 0;
}
      
      

  
						  
