/*!
 * @file vector.cpp
 * @brief scale and add vector by scalar
 * @author pistack (Junho Lee)
 * @date 2021. 10. 12.
 */

#include <vector>

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
