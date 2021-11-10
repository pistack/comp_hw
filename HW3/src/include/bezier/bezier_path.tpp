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
    // variable to check initial and final value of path
    // should be zero
    bool c0=false, cn=false; 
    const T eps = 100*std::numeric_limits<T>::epsilon();
    scale2 = 0; // initialize second scale parameter to 0
    vaildity = false; // initialize vaildity to false

    // check initial or final value of path
    // should be zero
    if(std::abs(p_init)< eps*(1+std::abs(p_final)))
    c0 = true;
    if(std::abs(p_final)< eps*(1+std::abs(p_init)))
    cn = true;

    // if initial value of path should be zero
    // set first control points to zero
    if(c0)
    B.c[0] = 0;
    // if finial value of path should be zero
    // set last control points to zero
    if(cn)
    B.c[B.n] = 0;

    // bezier curve is invalid only when
    // 1. initial value of path should not be zero but
    // first control point is zero.
    // 2. finial value of path should not be zero but
    // last control points is zero.
    if(!c0 && std::abs(B.c[0]) < eps*(1+std::abs(B.c[B.n])))
    return;
    if(!cn && std::abs(B.c[B.n]) < eps*(1+std::abs(B.c[0])))
    return;

    vaildity = true;
    // set first scaling parameter which modifies control points
    // and set second scaling parameter
    if(c0 && cn)
    scale2 = 1;
    else if(c0 && !cn)
    scale2 = p_final/B.c[B.n];
    else if(!c0 && cn)
    scale2 = p_init/B.c[0];
    else
    {
        B.c[B.n] = p_final/p_init*B.c[0];
        scale2 = p_init/B.c[0];
    }
    return;
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