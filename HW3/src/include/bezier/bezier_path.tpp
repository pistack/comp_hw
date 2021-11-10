/*!
 * @file bezier_path.tpp
 * @ingroup libpath
 * @brief define and evaluates path approximated by bezier curve
 * @author pistack (Junho Lee)
 * @date 2021. 11. 10.
 */

namespace libpath {

template<typename T>
void bezier_path<T>::init_helper()
{
    T c0=1, cn=1; // first and last control points
    const T eps = 100*std::numeric_limits<T>::epsilon();
    vaildity = false; // initialize vaildity to false
}

template<typename T>
std::vector<T> bezier_path<T>::eval(std::vector<T> t)
{
    T diff = t_final - t_init;
    std::vector<T> result = t;
    std::transform(result.begin(), result.end(), result.begin(),
    [t_init, diff](T &x){return (x-t_init)/diff;});
    result = B.eval(result);
    std::transform(result.begin(), result.end(),
    result.begin(), [scale2](T &x){return x *= scale2;});
    return result;
}

template<typename T>
std::vector<T> bezier_path<T>::deriv(std::vector<T> t)
{
    T diff = t_final - t_init;
    T scale_deriv = scale2/diff;
    std::vector<T> result = t;
    std::transform(result.begin(), result.end(), result.begin(),
    [t_init, diff](T &x){return (x-t_init)/diff;});
    result = B.deriv(result);
    std::transform(result.begin(), result.end(),
    result.begin(), [scale_deriv](T &x){return x *= scale_deriv;});
    return result;
}

}