/*!
 * @file hw2.tpp
 * @brief template for homework2 of Computer1 class in Yonsei University
 * Use numerical integration to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 10. 28.
 */

template<typename T>
std::tuple<std::vector<T>, std::vector<T>>
HW2(T zeta_min, T t0, int n)
{
  //initialize variable
  T grid_space; // grid spacing
  T a; 
  T zeta_max; // maximum value of zeta
  T tmp, tmp2;
  T c1, c2;
  std::vector<T> t(n+1, 0);
  std::vector<T> zeta(n+1, 0);
  std::vector<T> u(n+1, 0);
  std::vector<T> u_mid(n, 0);
  std::vector<T> integral(n+1, 0);

  // Note a = zeta_min**(-2)-2*zeta_min**(-1)
  //        = zeta_max**(-2)-2*zeta_max**(-1)

  tmp = 1/zeta_min;
  a = tmp*(tmp-2);
  tmp = 1.0+sqrt(1.0+a);
  zeta_max = -tmp/a;

  // calculates c1 and c2
  tmp = zeta_max - zeta_min;
  c1 = 2*zeta_min/tmp;
  c2 = 2*zeta_max/tmp;
  tmp = std::sqrt(-a);
  c1 /= tmp;
  c2 /= tmp;

  // use uniform n pts between zeta_min and zeta_max
  // additional one point needed for end point.

  grid_space = (zeta_max-zeta_min)/T(n);

  // change of variable
  // and seperate integral
  // Note v_i == u_(n-i)
  // v_(mid,i) == u_(mid,n-1-i)
  for(int i=1; i<n; i++)
    {
      zeta[i] = zeta[i-1] + grid_space;
      u[i] = std::sqrt(zeta[i]);
      u_mid[i-1] = std::sqrt(0.5*(zeta[i]+zeta[i-1]));
    }

  // end point
  zeta[n] = zeta_max-zeta_min;
  u[n] = std::sqrt(zeta[n]);
  u_mid[n-1] = std::sqrt(0.5*(zeta[n]+zeta[n-1]));

  // numerical integration using 
  // non equi-spacing Simpson's rule
  // Quite fun, integrand_1(u) = v and
  // integrand_2(v) = u
  // integral_1[0] = integral_2[0]
  // integral_1[n+1-i] = integral_2[i]
  for(int i=1; i<n+1; i++)
    {
      T h0 = u_mid[i-1]-u[i-1];
      T h1 = u[i] - u_mid[i-1];
      T d = u[i] - u[i-1];

      tmp = d/6.0*((2.0-h1/h0)*u[n+1-i]+
      (d/h0)*(d/h1)*u_mid[n-i]+
      (2.0-h0/h1)*u[n-i]);
      integral[i] = tmp;
    }

  for(int i=1; i<n+1; i++)
  t[i] = t[i-1] + c1*integral[i] + c2*integral[n+1-i];

  // initial condition
  // 1. t_0 = t0 
  // 2. zeta_0 = zeta_min
  // 3. zeta_n = zeta_max
  std::transform(t.begin(), t.end(), t.begin(),
  [t0](T &x){x += t0;});
  std::transform(zeta.begin(), zeta.end(), zeta.begin(),
  [zeta_min](T &x){x += zeta_min;});
  zeta[n] = zeta_max;

  return std::make_tuple(t, zeta);
}
