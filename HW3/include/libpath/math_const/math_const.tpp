/*!
 * @file math_const.tpp
 * @brief template implementation for mathematical constant
 * @author pistack (Junho Lee)
 * @date 2021. 11. 15.
*/

namespace libpath {
template<typename T>
constexpr T PI()
{
    unsigned int i = 1;
    T h = 1/T(16);
    T value = T(47)/T(15);
    T corr = 0;
    T term = 0;
    // use Simon Plouffe's infinite series on pi to
    // evaluate pi.
    do
    {
        term = h*(T(4)/(T(8*i+1))-T(2)/(T(8*i+4))-T(1)/(T(8*i+5))-T(1)/T(8*i+6));
        corr += term;
        h /= 16;
        ++i;
    } 
    while(term > corr*std::numeric_limits<T>::epsilon());
    // we know term > 0;
    return value+corr;
}

template<typename T>
constexpr T EXP1()
{
    unsigned int i = 7;
    T value = T(1957)/T(720);
    T corr = 0;
    T term = T(1)/T(720);
    do
    {
        term /= i;
        ++i;
        corr += term;
    } 
    while (term > corr*std::numeric_limits<T>::epsilon());
    return value+corr;
}
}
