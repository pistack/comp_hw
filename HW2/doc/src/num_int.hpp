////////////////////////////////////////////////////////////////////////
/// @page num_int Numerical Integration
/// @tableofcontents
/// By conservation of energy, we can derive following integral
/// equation.
/// \f{equation}{\label{eq:int_eq}
/// t - t_0 = \int_{\zeta_{min}}^{\zeta(t)}
/// \frac{\zeta'}{\sqrt{-a(\zeta'-\zeta_{min})(\zeta_{max}-\zeta')}}
/// \mathrm{d} \zeta'
/// \f}
/// ,where 
/// - \f$\zeta_{min}\f$ is periapsis (minimum value of \f$ \zeta\f$),
/// - \f$\zeta_{max}\f$ is apoapsis (maximum value of \f$ \zeta \f$)
/// - \f$ \alpha \f$ is parameter defined by following relation
/// \f{align*}{
/// \alpha &= \frac{1}{\zeta_{min}^2} - \frac{2}{\zeta_{min}} \\
///   &= \frac{1}{\zeta_{max}^2} - \frac{2}{\zeta_{max}}
/// \f}
/// To solve above integral equation 
/// \latexonly \eqref{eq:int_eq} \endlatexonly
/// we need to view time \f$ t \f$ as a function of \f$ \zeta \f$ with
/// domain \f$ D = [\zeta_{min}, \zeta_{max}] \f$.
/// Now uniformly divide the domain \f$ D \f$ into \f$ n \f$ sub intervals.
/// Let \f$ \zeta_i \f$ be end points of the sub intervals then
/// for \f$ 0 \leq \zeta \leq n \f$,
/// \f{equation}{\label{eq:zeta}
/// \zeta_i = \zeta_{min} + i \frac{\zeta_{max}-\zeta_{min}}{n}
/// \f}
/// Define \f$ t_i = t(\zeta_i) \f$ then we have following recurrence relation
/// for \f$ i \geq 1\f$,
/// \f{equation}{\label{eq:t_recurr}
/// t_i = t_{i-1} + 
/// \int_{\zeta_{i-1}}^{\zeta_i}
/// \frac{\zeta'}{\sqrt{-\alpha(\zeta'-\zeta_{min})(\zeta_{max}-\zeta')}}
/// \mathrm{d} \zeta'
/// \f}
/// However, due to the divergence feature of the integrand at 
/// \f$ \zeta_0 = \zeta_{min} \f$ and \f$ \zeta_n = \zeta_{max} \f$,
/// It is hard to approximate such integral directly.
/// To remove singularity, first consider the following equation.
/// \f{align}{
/// \frac{\zeta'}{\sqrt{-\alpha(\zeta'-\zeta_{min})(\zeta_{max}-\zeta')}} &=
/// \frac{\zeta_{min}}{\sqrt{-\alpha}(\zeta_{max}-\zeta_{min})}
/// \frac{\sqrt{\zeta_{max}-\zeta'}}{\sqrt{\zeta'-\zeta_{min}}} 
/// \label{eq:sep} \\
/// &+ \frac{\zeta_{max}}{\sqrt{-\alpha}(\zeta_{max}-\zeta_{min})}
/// \frac{\sqrt{\zeta'-\zeta_{min}}}{\sqrt{\zeta_{max}-\zeta'}}
/// \nonumber
/// \f}
/// Then we can seperate integral into two parts.
/// Substitute \f$ u = \sqrt{\zeta'-\zeta_{min}} \f$ and
/// \f$ v = \sqrt{\zeta_{max}-\zeta'} \f$ to
/// the first and second part of integral respectively, then
/// \f{align}{
/// \mathrm{integral}_i &= 
/// \frac{2\zeta_{min}}{\sqrt{-\alpha}(\zeta_{max}-\zeta_{min})}
/// \int_{\sqrt{\zeta_{i-1}-\zeta_{min}}}^{\sqrt{\zeta_{i}-\zeta_{min}}}
/// \sqrt{\zeta_{max}-\zeta_{min}-u^2} \mathrm{d} u \label{eq:substitue} \\
/// &+ \frac{2\zeta_{max}}{\sqrt{-\alpha}(\zeta_{max}-\zeta_{min})}
/// \int_{\sqrt{\zeta_{max}-\zeta_{i}}}^{\sqrt{\zeta_{max}-\zeta_{i-1}}}
/// \sqrt{\zeta_{max}-\zeta_{min}-v^2} \mathrm{d} v
/// \f}
/// by above equation \latexonly \eqref{eq:substitue} \endlatexonly,
/// We can deduce
/// \f{align}{
/// t_f - t_0 &= \pi \frac{\zeta_{min}+\zeta_{max}}{\sqrt{-\alpha}} 
/// \label{eq:tmax} \\
/// &= \pi a^{3/2} \nonumber \f}
/// ,where \f$ a = (\zeta_{min}+\zeta_{max})/2 \f$.
/// @section approx Approximation
/// Let \f$ f(u) = \sqrt{\zeta_{max}-\zeta_{min}-u^2} \f$ and define
/// \f$ u_i \f$ and \f$ v_i \f$ as following
/// \f{align}{
/// u_i &= \sqrt{\zeta_i - \zeta_{min}} \label{eq:u} \\
/// v_i &= \sqrt{\zeta_{max} - \zeta_i} \label{eq:v}
///\f}
/// then we have following relation
/// \f{align}{
/// u_i &= f(v_i) \label{eq:fv} \\
/// v_i &= f(u_i) \label{eq:fu}
///\f}
/// To exploit above relation \latexonly \eqref{eq:fv}--\eqref{eq:fu} \endlatexonly
/// and gain more accurate results, I use Simpson's Rule for unequally spaced 
/// ordinates (https://www.jstor.org/stable/2309244).
/// Let 
/// \f{align*}{
/// u_{i-1/2} &= \sqrt{(\zeta_{i-1}+\zeta_i)/2-\zeta_{min}} \\
/// v_{i-1/2} &= \sqrt{\zeta_{max}-(\zeta_{i-1}+\zeta_i)/2} \\
/// h_0^{u} &= u_{i-1/2} - u_{i-1} \\
/// h_0^{v} &= v_{i-1/2} - v_{i-1} \\
/// h_1^{u} &= u_{i} - u_{i-1/2} \\
/// h_1^{v} &= v_{i} - v_{i-1/2}
/// \f}
/// then,
/// \f{align}{
/// \mathrm{integral}_i &= c_1\frac{u_i-u_{i-1}}{6}
/// \left[\left(2-\frac{h_1^{u}}{h_0^{u}}\right)v_{i-1} +
/// \frac{(u_i-u_{i-1})^2}{h_0^uh_1^u}v_{i-1/2} +
/// \left(2-\frac{h_0^u}{h_1^u}\right)v_i
///\right] \label{eq:noneq_simple} \\
/// &- c_2 \cdot \text{(swap $u$ and $v$)} \\ \nonumber
///\f}
/// where
/// \f{align*}{
/// c_1 &= \frac{2\zeta_{min}}{\sqrt{-\alpha}(\zeta_{max}-\zeta_{min})} \\
/// c_2 &= \frac{2\zeta_{max}}{\sqrt{-\alpha}(\zeta_{max}-\zeta_{min})}
/// \f}
/// @section complexity Complexity
/// Complex is clearly 
/// \f{equation*}{O(n)\f}
/// @section accuracy Accuracy
/// Error bound is given by
/// \f{align*}{
/// \text{Error bound} &\leq M \sum_{i=1}^n (u_i-u_{i-1})^5  \\
/// &= M \left(\frac{\zeta_{max}-\zeta_{min}}{n}\right)^{5/2} 
/// \sum_{i=1}^n \left(\frac{1}{\sqrt{i}+\sqrt{i-1}}\right)^{5/2} \\
/// &< MC \left(\frac{\zeta_{max}-\zeta_{min}}{n}\right)^{5/2}
///\f}
/// So, the error bound is
/// \f{equation*} O(n^{-5/2}) \f}
/// @section conv Convergence
/// initial condition
/// \f{align*}{
/// t_0 &= 0 \\
/// \zeta_{min} &= 0.9
/// \f}
/// @image html plot_converge.png "Convergence plot"
/// @image latex plot_converge.eps "Convergence plot" width=10cm
//////////////////////////////////////////////////////////////////////////