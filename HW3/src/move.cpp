/*!
 * @file move.cpp
 * @brief randomly move initial guess, at most step size.
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 */

#include <vector>
#include <random>

using namespace std;

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
