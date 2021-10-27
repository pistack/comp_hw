/*!
 * @file fourier.tpp
 * @ingroup libfourier
 * @brief evaluate sum and derivative of sine and cosine function
 * @author pistack (Junho Lee)
 * @date 2021. 10. 28.
 */

template<typename T>
fourier<T>& fourier<T>::operator=(const fourier<T> &copy)
{
	f_num_fourier = copy.f_num_fourier;
	f_period = copy.f_period;
	f_c = copy.f_c;
	return *this;
}

template<typename T>
void fourier<T>::update(std::vector<T> c){ f_c = c;}

template<typename T>
T fourier<T>::eval(T t)
{
	int term = 2*f_num_fourier;
	T omega = 2*pi/f_period;
	T tmp;
	T y=0;
	
	for(int j=0; j<term; j += 2)
	{
		tmp = (j/2+1)*omega;
		if(f_c[j] != 0.0)
		y += f_c[j]*std::sin(tmp*t);
		if(f_c[j+1] != 0.0)
		y += f_c[j+1]*std::cos(tmp*t);
	}

	return y;
}

template<typename T>
std::vector<T> fourier<T>::eval(std::vector<T> t)
{
	int term = 2*f_num_fourier;
	int n = t.size();
	T omega = 2*pi/f_period;
	T tmp;
	std::vector<T> y(n, 0);
	
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<term; j += 2)
		{
			tmp = (j/2+1)*omega;
			if(f_c[j] != 0.0)
			y[i] += f_c[j]*std::sin(tmp*t[i]);
			if(f_c[j+1] != 0.0)
			y[i] += f_c[j+1]*std::cos(tmp*t[i]);
		 }
	}

	return y;
}

template<typename T>
T fourier<T>::deriv(T t)
{
	int term = 2*f_num_fourier;
	T omega = 2*pi/f_period;
	T tmp;
	T yp=0;
	
	for(int j=0; j<term; j += 2)
	{
		tmp = (j/2+1)*omega;
		if(f_c[j] != 0.0)
		yp += f_c[j]*tmp*std::cos(tmp*t);
		if(f_c[j+1] != 0.0)
		yp -= f_c[j+1]*tmp*std::sin(tmp*t);
	}

	return yp;
}

template<typename T>
std::vector<T> fourier<T>::deriv(std::vector<T> t)
{
	int term = 2*f_num_fourier;
	int n = t.size();
	T omega = 2*pi/f_period;
	T tmp;
	std::vector<T> yp(n, 0);
	
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<term; j += 2)
		{
			tmp = (j/2+1)*omega;
			if(f_c[j] != 0.0)
			yp[i] += f_c[j]*tmp*std::cos(tmp*t[i]);
			if(f_c[j+1] != 0.0)
			yp[i] -= f_c[j+1]*tmp*std::sin(tmp*t[i]);
		}
	}

	return yp;
}

template<typename T>
T fourier<T>::nderiv(int n, T t)
{
	int term = 2*f_num_fourier;
	T omega = 2*pi/f_period;
	T tmp;
	std::vector<T> c = f_c;
	T yp=0;

	if(n % 2 == 0)
	{}
	else
	{
		for(int j=0; j<term; j+=2)
		{
			std::iter_swap(c.begin()+j, c.begin()+j+1);
			c[j] *= -1.0;
		}
	}
	for(int j=0; j<term; j += 2)
	{
		tmp = std::pow((j/2+1)*omega, T(n));
		if(c[j] != 0.0)
		yp += c[j]*tmp*std::sin(tmp*t);
		if(f_c[j+1] != 0.0)
		yp += f_c[j+1]*tmp*std::cos(tmp*t);
	}
	if((n%2 == 0 && n/2 % 2 != 0) ||
	(n%2 != 0 && (n-1)/2 % 2 != 0))
	yp *= -1.0;
	return yp;
}

template<typename T>
std::vector<T> fourier<T>::nderiv(int n, std::vector<T> t)
{
	int term = 2*f_num_fourier;
	int n_t = t.size();
	T omega = 2*pi/f_period;
	T tmp;
	std::vector<T> c = f_c;
	std::vector<T> yp(n, 0);

	if(n % 2 == 0)
	{
		if(n/2 %2 == 0){}
		else
		{
			std::transform(c.begin(), c.end(), c.begin(), 
			[](T &x){return x *= -1.0;});
		}
	}
	else
	{
		for(int j=0; j<term; j+=2)
		{
			std::iter_swap(c.begin()+j, c.begin()+j+1);
			c[j] *= -1.0;
		}
		if((n-1)/2 %2 == 0){}
		else
		{
			std::transform(c.begin(), c.end(), c.begin(), 
			[](T &x){return x *= -1.0;});
		}
	}
	
	for(int i=0; i<n_t; i++)
	{
		for(int j=0; j<term; j += 2)
		{
			tmp = std::pow((j/2+1)*omega, T(n));
			if(c[j] != 0.0)
			yp[i] += c[j]*tmp*std::sin(tmp*t[i]);
			if(f_c[j+1] != 0.0)
			yp[i] += f_c[j+1]*tmp*std::cos(tmp*t[i]);
		}
	}

	return yp;
}
