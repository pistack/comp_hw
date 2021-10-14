///////////////////////////////////////////////////////////////////////
/*! @page theta Theta
 * @tableofcontents
 * By conservation of angular momentum,
 * Angle \f$ \theta \f$ satisfies following relation
 * \f{equation}{\label{eq:theta_ode}
 * \frac{\mathrm{d} \theta}{\mathrm{d} t} = \frac{1}{\zeta^2}
 * \f}
 * Integrate both side then we can deduce
 * \f{equation}{\label{eq:theta_int}
 * \theta(t) = \theta_0 + \int_{t_0}^t \frac{1}{\zeta^2} \mathrm{d} t
 * \f}
 * Let \f$ \theta_i = \theta(t_i) \f$ as in \ref fdm, then
 * for \f$ 1 \leq i \f$,
 * \f{equation}{\label{eq:theta_recurr}
 * \theta_i = \theta_{i-1} + \int_{t_{i-1}}^{t_i} \frac{1}{\zeta^2} \mathrm{d} t
 * \f}
 * Next approximate the integral using trapezoidal rule then
 * \f{equation}{\label{eq:theta_approx} 
 * \theta_i \approx \theta_{i-1} + 
 * \frac{t_i-t_{i-1}}{2}\left(\frac{1}{\zeta_{i-1}^2}+\frac{1}{\zeta_i^2}\right)
 * \f}
 * \f$\theta_i \f$ has \f$ O(n^{-3}) \f$ local turncation error for 
 * trapezoidal rule and additional \f$ O(n^{-3}) \f$ for the global turncation error of
 * \f$ \zeta \f$ (see \ref fdm \ref accurcay).
 * So the global turncation error of \f$ \theta \f$ can be estimated to
 * \f$ O(n^{-2}) \f$
 * @section complexity_theta Complexity
 * Clearly
 * \f{equation}{\label{eq:complex_theta}
 * O(n)
 * \f}
 * @section accuracy_theta Accuracy
 * The global turncation error of \f$ \theta \f$ is roughtly estimated to
 * \f{equation}{\label{error_theta}
 * O(n^{-2}) \f}
 * @section conv Convergence
 * - Initial Condition
 * \f{align*}{
 * \zeta(0) &= 0.9 \\
 * \zeta'(0) &= 0 \\
 * \theta(0) &= 0
 * \f}
 * - Initial time: 0
 * - Final time: 10
 * @image html plot_converge_theta.png "Convergence plot"
 * @image latex plot_converge_theta.eps "Convergence plot" width=10cm
 */
///////////////////////////////////////////////////////////////////////////////////////////