/*!
 * @file support.cpp
 * @brief support functions for 
 *        homework3 of Computer1 class in Yonsei University
 *        scale and add vector, evaluate sum and derivative of sine function,
 *        randomly move initial guess by step and evaluate the action of
 *        given path.
 * @author pistack (Junho Lee)
 * @date 2021. 10. 10.
 */

#include "hw3.hpp"

using namespace std;

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

vector<double>
move_step(vector<double> &init_guess, double step, mt19937 &gen,
	  uniform_real_distribution<double> &dist)
{
  int n = init_guess.size();
  vector<double> move(n, 0);
  for(int i=0; i<n; i++)
    {
      move[i] = init_guess[i] + step*dist(gen);
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
