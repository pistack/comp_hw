/*!
 * @file hw2.tpp
 * @brief template for homework2 of Computer1 class in Yonsei University
 * Use numerical integration to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 11. 10.
 */

template<typename T>
std::tuple<std::vector<T>, std::vector<T>>
HW2(T zeta_min, T t0, std::size_t n)
{
  //initialize variable
  T grid_space; // grid spacing
  T a; 
  T zeta_max; // maximum value of zeta
  T tmp;
  T c1, c2;
  std::vector<T> t(n+1, 0);
  std::vector<T> zeta(n+1, 0);
  std::vector<T> u(n+1, 0);

  // Note.
  // zeta_max = zeta_min/(2*zeta_min-1)
  // a = - 1/(zeta_max*zeta_min)
  // by Vieta's Formula

  zeta_max = zeta_min/(2*zeta_min-1);
  a = 1/(zeta_min*zeta_max); // use -a instead

  // calculates c1 and c2
  tmp = zeta_max - zeta_min;
  c1 = 2*zeta_min/tmp;
  c2 = 2*zeta_max/tmp;
  tmp = std::sqrt(a);
  c1 /= tmp;
  c2 /= tmp;

  // use uniform n pts between zeta_min and zeta_max
  // additional one point needed for end point.

  grid_space = (zeta_max-zeta_min)/T(n);

  // change of variable
  // and seperate integral
  // Note v_i == u_(n-i)
  // v_(mid,i) == u_(mid,n-1-i)
  for(std::size_t i=1; i<n; ++i)
    {
      zeta[i] = zeta[i-1] + grid_space;
      u[i] = std::sqrt(zeta[i]);
    }

  // end point
  zeta[n] = zeta_max-zeta_min;
  u[n] = std::sqrt(zeta[n]);

  // numerical integration using 
  // non equi-spacing trapezoidal rule
  // Quite fun, integrand_1(u) = v and
  // integrand_2(v) = u
  // integral_1[0] = integral_2[0]
  // integral_1[n+1-i] = integral_2[i]
  for(std::size_t i=1; i<n+1; ++i)
  {
    T integral_part1;
    T integral_part2;
    integral_part1 = (u[i]-u[i-1])/2.0*(u[n-i]+u[n+1-i]);
    integral_part2 = (u[n+1-i]-u[n-i])/2.0*(u[i]+u[i-1]);
    t[i] = t[i-1] + c1*integral_part1 + c2*integral_part2;
  }

  // initial condition
  // 1. t_0 = t0 
  // 2. zeta_0 = zeta_min
  std::transform(t.begin(), t.end(), t.begin(),
  [t0](T &x){return x += t0;});
  std::transform(zeta.begin(), zeta.end(), zeta.begin(),
  [zeta_min](T &x){return x += zeta_min;});

  return std::make_tuple(t, zeta);
}
