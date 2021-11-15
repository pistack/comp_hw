/*!
 * @file math_const.tpp
 * @brief template implementation of math_const
 * @author pistack (Junho Lee)
 * @date 2021. 11. 15.
*/

namespace libpath {
template<typename T>
constexpr T PI()
{
    unsigned int i = 0;
    T h = 1;
    T value = 0;
    T term = 0;
    // use Simon Plouffe's infinite series on pi to
    // evaluate pi.
    while(h >std::numeric_limits<T>::epsilon()*std::numeric_limits<T>::epsilon())
    do
    {
        term = h*(T(4)/(T(8*i+1))-T(2)/(T(8*i+4))-T(1)/(T(8*i+5))-T(1)/T(8*i+6));
        value += term;
        h /= 16;
        i++;
    } 
    while(term > value*std::numeric_limits<T>::epsilon());
    // we know term > 0;
    return value;
}
}
