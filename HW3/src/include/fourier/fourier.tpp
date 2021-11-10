/*!
 * @file fourier.tpp
 * @ingroup libpath
 * @brief evaluate sum and derivative of sine and cosine function
 * @author pistack (Junho Lee)
 * @date 2021. 11. 10.
 */

namespace libpath {

template<typename T>
T fourier<T>::eval(T t)
{
	int term = 2*f_num_fourier;
	T omega = 2*pi<T>/f_period;
	T tmp=0;
	T y=0;
	
	for(int j=0; j<term; ++++j)
	{
		tmp += omega;
		y += f_c[j]*std::sin(tmp*t);
		y += f_c[j+1]*std::cos(tmp*t);
	}

	return y;
}

template<typename T>
std::vector<T> fourier<T>::eval(std::vector<T> t)
{
	int term = 2*f_num_fourier;
	int n = t.size();
	T omega = 2*pi<T>/f_period;
	std::vector<T> y(n, 0);
	
	for(int i=0; i<n; ++i)
	{
		T tmp=0;
		for(int j=0; j<term; ++++j)
		{
			tmp += omega;
			y[i] += f_c[j]*std::sin(tmp*t[i]);
			y[i] += f_c[j+1]*std::cos(tmp*t[i]);
		 }
	}

	return y;
}

template<typename T>
T fourier<T>::deriv(T t)
{
	int term = 2*f_num_fourier;
	T omega = 2*pi<T>/f_period;
	T tmp=0;
	T yp=0;
	
	for(int j=0; j<term; ++++j)
	{
		tmp += omega;
		yp += f_c[j]*tmp*std::cos(tmp*t);
		yp -= f_c[j+1]*tmp*std::sin(tmp*t);
	}

	return yp;
}

template<typename T>
std::vector<T> fourier<T>::deriv(std::vector<T> t)
{
	int term = 2*f_num_fourier;
	int n = t.size();
	T omega = 2*pi<T>/f_period;
	std::vector<T> yp(n, 0);
	
	for(int i=0; i<n; ++i)
	{
		T tmp = 0;
		for(int j=0; j<term; ++++j)
		{
			tmp += omega;
			yp[i] += f_c[j]*tmp*std::cos(tmp*t[i]);
			yp[i] -= f_c[j+1]*tmp*std::sin(tmp*t[i]);
		}
	}
	return yp;
}

template<typename T>
T fourier<T>::nderiv(int n, T t)
{
	int term = 2*f_num_fourier;
	T omega = 2*pi<T>/f_period;
	T tmp=0;
	T yp=0;
	std::vector<T> c(term, 0);
	c = f_c;

	if(n % 2 == 0)
	{}
	else
	{
		for(typename std::vector<T>::iterator iter = c.begin(); 
		iter != c.end(); ++++iter)
		{
			std::iter_swap(iter, iter+1);
			*iter *= -1;
		}
	}
	for(int j=0; j<term; ++++j)
	{
		T tmpower;
		tmp += omega;
		tmpower = std::pow(tmp, T(n));
		yp += c[j]*tmpower*std::sin(tmp*t);
		yp += c[j+1]*tmpower*std::cos(tmp*t);
	}
	if((n%2 == 0 && n/2 % 2 != 0) ||
	(n%2 != 0 && (n-1)/2 % 2 != 0))
	yp *= -1;
	return yp;
}

template<typename T>
std::vector<T> fourier<T>::nderiv(int n, std::vector<T> t)
{
	int term = 2*f_num_fourier;
	int n_t = t.size();
	T omega = 2*pi<T>/f_period;
	std::vector<T> c(term, 0);
	std::vector<T> yp(n_t, 0);

	c = f_c;

	if(n % 2 == 0)
	{
		if(n/2 %2 == 0){}
		else
		{
			std::transform(c.begin(), c.end(), c.begin(), 
			[](T &x){return x *= -1;});
		}
	}
	else
	{
		for(typename std::vector<T>::iterator iter = c.begin(); 
		iter != c.end(); ++++iter)
		{
			std::iter_swap(iter, iter+1);
			*iter *= -1;
		}
		if((n-1)/2 %2 == 0){}
		else
		{
			std::transform(c.begin(), c.end(), c.begin(), 
			[](T &x){return x *= -1;});
		}
	}
	
	for(int i=0; i<n_t; ++i)
	{
		T tmp = 0;
		for(int j=0; j<term; ++++j)
		{
			T tmpower;
			tmp += omega;
			tmpower = std::pow(tmp, T(n));
			yp[i] += c[j]*tmpower*std::sin(tmp*t[i]);
			yp[i] += c[j+1]*tmpower*std::cos(tmp*t[i]);
		}
	}
	return yp;
}
}

