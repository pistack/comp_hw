/*!
 * @file fourier_path.tpp
 * @ingroup libpath
 * @brief evaluate path and derivative
 * @author pistack (Junho Lee)
 * @date 2021. 11. 12.
 */

namespace libpath {
template<typename T>
void fourier_path<T>::init_helper()
{
  // init variable
  T tstp, tstf;
  std::vector<T> t_ends = {p_t0, p_tf};
  std::vector<T> f_ends(2, 0);
  // practical mechine epsilon
  const T eps = 100*std::numeric_limits<T>::epsilon();

  f_ends = p_func.eval(t_ends);
  tstp = p_f - p_0;
  tstf = f_ends[1]-f_ends[0];

  T deltap = std::abs(tstp);
  T deltaf = std::abs(tstf);
  T pmean = (std::abs(p_f) + std::abs(p_0))/2;
  T fmean = (std::abs(f_ends[1])+std::abs(f_ends[0]))/2;

  if(deltap < eps*(1+pmean) && 
  deltaf < eps*(1+fmean))
  {
    p_vaild = true;
    scale = 1;
    add = p_0 - f_ends[0];
  }
  else if(deltap > eps*(1+pmean) 
  && deltaf > eps*(1+fmean))
  {
    p_vaild = true;
    scale = tstp/tstf;
    add = p_0 - scale*f_ends[0];
  }
  else
  {
    p_vaild = false;
    scale = 0;
    add = 0;
  }
}

template<typename T>
std::vector<T> fourier_path<T>::eval(std::vector<T> t)
{
  // init variable
  T scaler = scale;
  T adder = add;
  std::vector<T> result = t;
  
  result = p_func.eval(t);

  // scale and add
  std::transform(result.begin(), result.end(), result.begin(),
  [scaler, adder](T &x){x *= scaler; 
      x += adder;
      return x;});

  return result;
}

template<typename T>
std::vector<T> fourier_path<T>::deriv(std::vector<T> t)
{
  // init variable
  T scaler = scale;
  std::vector<T> result = t;
  
  result = p_func.deriv(t);

  std::transform(result.begin(), result.end(), result.begin(),
  [scaler](T &x){return x *= scaler;});

  return result;
}

template<typename T>
std::vector<T> fourier_path<T>::nderiv(unsigned int n, std::vector<T> t)
{
  // init variable
  T scaler = scale;
  std::vector<T> result = t;
  
  result = p_func.nderiv(n, t);

  std::transform(result.begin(), result.end(), result.begin(),
  [scaler](T &x){return x *= scaler;});

  return result;
}
}
  
