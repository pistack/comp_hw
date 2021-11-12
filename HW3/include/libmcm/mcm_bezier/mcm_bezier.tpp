/*!
 * @file mcm_bezier.tpp
 * @ingroup libmcm
 * @brief template implementation for mcm_bezier class
 * which is derivated by libmcm::mcm class
 * @author pistack (Junho Lee)
 * @date 2021. 11. 12.
 */

namespace libmcm {

// get_init_coeff
template<typename T, typename Lag>
std::tuple<std::vector<T>, std::vector<std::vector<T>>>
mcm_bezier<T, Lag>::get_init_coeff()
{
    std::vector<libpath::bezier_path<T>> init_path_ = \
    mcm<T, libpath::bezier<T>, libpath::bezier_path<T>, Lag>::get_init_path();
    unsigned int dim1 = init_path_.size();
    unsigned int dim2 = init_path_[0].get_ctrl_pts().size();
    std::vector<T> scaler(dim1, 0);
    std::vector<std::vector<T>> ctrl_pts(dim1, std::vector<T>(dim2, 0));
    for(unsigned int i=0; i<dim1; ++i)
    {
        scaler[i] = init_path_[i].get_scaler();
        ctrl_pts[i] = init_path_[i].get_ctrl_pts(); 
    }
    return std::make_tuple(scaler, ctrl_pts);
}

// get_min_coeff
template<typename T, typename Lag>
std::tuple<std::vector<T>, std::vector<std::vector<T>>>
mcm_bezier<T, Lag>::get_min_coeff()
{
    std::vector<libpath::bezier_path<T>> min_path_ = \
    mcm<T, libpath::bezier<T>, libpath::bezier_path<T>, Lag>::get_min_path();
    unsigned int dim1 = min_path_.size();
    unsigned int dim2 = min_path_[0].get_ctrl_pts().size();
    std::vector<T> scaler(dim1, 0);
    std::vector<std::vector<T>> ctrl_pts(dim1, std::vector<T>(dim2, 0));
    for(unsigned int i=0; i<dim1; ++i)
    {
        scaler[i] = min_path_[i].get_scaler();
        ctrl_pts[i] = min_path_[i].get_ctrl_pts(); 
    }
    return std::make_tuple(scaler, ctrl_pts);
}

}