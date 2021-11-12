/*!
 * @file mcm_fourier.tpp
 * @ingroup libmcm
 * @brief template implementation for mcm_fourier class
 * which is derivated by libmcm::mcm class
 * @author pistack (Junho Lee)
 * @date 2021. 11. 12.
 */

namespace libmcm {

// get_init_coeff
template<typename T, typename Lag>
std::tuple<std::vector<T>, std::vector<T>, std::vector<std::vector<T>>>
mcm_fourier<T, Lag>::get_init_coeff()
{
    std::vector<libpath::fourier_path<T>> init_path_ = \
    mcm<T, libpath::fourier<T>, libpath::fourier_path<T>, Lag>::get_init_path();
    unsigned int dim1 = init_path_.size();
    unsigned int dim2 = init_path_[0].get_coeff().size();
    std::vector<T> adder(dim1, 0);
    std::vector<T> scaler(dim1, 0);
    std::vector<std::vector<T>> coeff(dim1, std::vector<T>(dim2, 0));
    for(unsigned int i=0; i<dim1; ++i)
    {
        adder[i] = init_path_[i].get_adder();
        scaler[i] = init_path_[i].get_scaler();
        coeff[i] = init_path_[i].get_coeff(); 
    }
    return std::make_tuple(adder, scaler, coeff);
}

// get_min_coeff
template<typename T, typename Lag>
std::tuple<std::vector<T>, std::vector<T>, std::vector<std::vector<T>>>
mcm_fourier<T, Lag>::get_min_coeff()
{
    std::vector<libpath::fourier_path<T>> min_path_ = \
    mcm<T, libpath::fourier<T>, libpath::fourier_path<T>, Lag>::get_min_path();
    unsigned int dim1 = min_path_.size();
    unsigned int dim2 = min_path_[0].get_coeff().size();
    std::vector<T> adder(dim1, 0);
    std::vector<T> scaler(dim1, 0);
    std::vector<std::vector<T>> coeff(dim1, std::vector<T>(dim2, 0));
    for(unsigned int i=0; i<dim1; ++i)
    {
        adder[i] = min_path_[i].get_adder();
        scaler[i] = min_path_[i].get_scaler();
        coeff[i] = min_path_[i].get_coeff(); 
    }
    return std::make_tuple(adder, scaler, coeff);
}

}