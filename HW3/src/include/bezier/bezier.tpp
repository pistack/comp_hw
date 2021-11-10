/*!
 * @file bezier.tpp
 * @ingroup libpath
 * @brief define and evaluates bezier curve
 * @author pistack (Junho Lee)
 * @date 2021. 11. 10.
 */

namespace libpath {

/// evaluates bezier curve of order n with 
/// De Casteljau's Algorithm
template<typename T>
T bezier<T>::eval(T t)
{
    std::vector<T> result = c; // store beta_i^{0}
    for(int j=0; j < n; j++)
    {
        for(int i=0; i<n-j; i++)
        result[i] = result[i]*(1-t) + result[i+1]*t;
    }
    return result[0];
}

template<typename T>
std::vector<T> bezier<T>::eval(std::vector<T> t)
{
    std::size_t dim = t.size();
    std::vector<T> result(dim, 0);
    for(std::size_t i=0; i<dim; i++)
    result[i] = eval(t[i]);

    return result;
}

/// evaluates derivative of bezier curve of order n with 
/// De Casteljau's Algorithm
/// Note derivative of bezier curve of order n is
/// n \times Bezier curve of order n-1 with
///  control points {c_{i+1}-c_i}.
template<typename T>
T bezier<T>::deriv(T t)
{
    std::vector<T> result(n, 0); 
    // store beta_i^{0}
    for(int i=0; i<n; i++)
    result[i] = c[i+1] - c[i];

    for(int j=0; j < n-1; j++)
    {
        for(int i=0; i<n-j-1; i++)
        result[i] = result[i]*(1-t) + result[i+1]*t;
    }
    return result[0];
}

template<typename T>
std::vector<T> bezier<T>::deriv(std::vector<T> t)
{
    std::size_t dim = t.size();
    std::vector<T> result(dim, 0);
    for(std::size_t i=0; i<dim; i++)
    result[i] = deriv(t[i]);
    return result;
}
}
