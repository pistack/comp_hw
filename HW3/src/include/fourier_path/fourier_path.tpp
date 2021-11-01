/*!
 * @file fourier_path.tpp
 * @ingroup libfourier
 * @brief evaluate path and derivative
 * @author pistack (Junho Lee)
 * @date 2021. 10. 30.
 */

template<typename T>
fourier_path<T> & fourier_path<T>::operator=(const fourier_path<T> & copy)
{
  // copy initial conditions
  p_t0 = copy.p_t0; p_tf = copy.p_tf; p_0 = copy.p_0; p_f = copy.p_f; 
  // copy fourier function, scaler and adder
	p_func = copy.p_func; p_vaild = copy.p_vaild; 
  scale = copy.scale; add = copy.add;
  return *this;
}

template<typename T>
void fourier_path<T>::init_helper()
{
  // init variable
  T tstp, tstf;
  std::vector<T> t_ends = {p_t0, p_tf};
  std::vector<T> f_ends(2, 0.0);
  // practical mechine epsilon
  T eps = 100.0*std::numeric_limits<T>::epsilon();

  f_ends = p_func.eval(t_ends);
  tstp = p_f - p_0;
  tstf = f_ends[1]-f_ends[0];

  T deltap = std::abs(tstp);
  T deltaf = std::abs(tstf);
  T pmean = 0.5*(std::abs(p_f) + std::abs(p_0));
  T fmean = 0.5*(std::abs(f_ends[1])+std::abs(f_ends[0]));

  if(deltap < eps*(1.0+pmean) && 
  deltaf < eps*(1.0+fmean))
  {
    p_vaild = true;
    scale = 1.0;
    add = p_0 - f_ends[0];
  }
  else if(deltap > eps*(1.0+pmean) 
  && deltaf > eps*(1.0+fmean))
  {
    p_vaild = true;
    scale = tstp/tstf;
    add = p_0 - scale*f_ends[0];
  }
  else
  {
    p_vaild = false;
    scale = 0.0;
    add = 0.0;
  }
}

template<typename T>
void fourier_path<T>::update(fourier<T> fourier)
{
  p_func = fourier;
  init_helper();
}

template<typename T>
bool fourier_path<T>::is_vaild()
{
  return p_vaild;
}

template<typename T>
T fourier_path<T>::get_adder()
{
  return add;
}

template<typename T>
T fourier_path<T>::get_scaler()
{
  return scale;
}

template<typename T>
std::tuple<T, T> fourier_path<T>::get_endtimes()
{
  return std::make_tuple(p_t0, p_tf);
}

template<typename T>
T fourier_path<T>::eval(T t)
{
  return scale*p_func.eval(t)+add;
}

template<typename T>
std::vector<T> fourier_path<T>::eval(std::vector<T> t)
{
  // init variable
  int n = t.size();
  T scaler = scale;
  T adder = add;
  std::vector<T> result(n, 0);
  
  result = p_func.eval(t);

  // scale and add
  std::transform(result.begin(), result.end(), result.begin(),
  [scaler, adder](T &x){x *= scaler; 
      x += adder;
      return x;});

  return result;
}

template<typename T>
T fourier_path<T>::deriv(T t)
{
  return scale*p_func.deriv(t);
}

template<typename T>
std::vector<T> fourier_path<T>::deriv(std::vector<T> t)
{
  // init variable
  int n = t.size();
  T scaler = scale;
  std::vector<T> result(n, 0);
  
  result = p_func.deriv(t);

  std::transform(result.begin(), result.end(), result.begin(),
  [scaler](T &x){return x *= scaler;});

  return result;
}

template<typename T>
T fourier_path<T>::nderiv(int n, T t)
{
  return scale*p_func.nderiv(n, t);
}

template<typename T>
std::vector<T> fourier_path<T>::nderiv(int n, std::vector<T> t)
{
  // init variable
  int n_t = t.size();
  T scaler = scale;
  std::vector<T> result(n_t, 0);
  
  result = p_func.nderiv(n, t);

  std::transform(result.begin(), result.end(), result.begin(),
  [scaler](T &x){return x *= scaler;});

  return result;
}
  
