/*!
 * @file hw1.tpp
 * @brief template for homework1 of Computer1 class in Yonsei University
 * Use finite difference method to solve Kepler problem
 * @author pistack (Junho Lee)
 * @date 2021. 11. 10.
 */

template<typename T>
std::tuple<std::vector<T>, std::vector<T>, std::vector<T>> 
HW1(T t0, T t1,
std::size_t n, T y0, T y0p, T theta0)
{ 
    // initialize the variable
    T spacing;
    T y_inv;
    T tmp;
    std::vector<T> t(n+1, 0);
    std::vector<T> y(n+1, 0);
    std::vector<T> theta(n+1, 0);

    spacing = (t1 - t0) / T(n);

    // use uniform n points bewteen t0 and t1
    // additional one point needed for end point
    // to avoid numerical problem add t0 later
    for(std::size_t i=1; i< n+1; ++i)
      t[i] = t[i-1] + spacing;

    // to avoid numerical problem add initial condition
    // later
    // estimate y1 using 2nd order talyor expension
    y_inv = 1/y0;
    spacing = t[1] - t[0];
    y[1] = spacing*(y0p+y_inv*(0.5*spacing*y_inv*
    (y_inv - 1.0)));
    tmp = 1.0/(1.0+y_inv*y[1]);
    // estimate theta1 using trapezoid rule
    theta[1] = theta[0] + \
    0.5*spacing*y_inv*y_inv * \
    (1.0+tmp*tmp);

    // Solve 2nd order ODE using
    // Explict Euler Method

    for(std::size_t i = 2; i < n+1; ++i)
    {
      spacing = (t[i] - t[i-2])/2.0; // spacing for zeta
      y_inv = 1 /(y0+y[i-1]);
      tmp = y_inv*spacing;
      y[i] = tmp*tmp * \
      (y_inv - 1.0) + \
      (2*y[i-1] - y[i-2]);
      spacing = t[i] - t[i-1]; // spacing for theta
      tmp = (y0+y[i-1])/(y0+y[i]);
      theta[i] = theta[i-1] + \
      0.5*spacing*y_inv*y_inv * \
      (1.0+tmp*tmp);
    }

    std::transform(t.begin(), t.end(), t.begin(),
    [t0](T &x){return x += t0;});

    std::transform(y.begin(), y.end(), y.begin(),
    [y0](T &x){return x += y0;});

    std::transform(theta.begin(), theta.end(), theta.begin(),
    [theta0](T &x){return x += theta0;});
    
    return std::make_tuple(t, y, theta);
}
