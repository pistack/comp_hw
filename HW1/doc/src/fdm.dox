///////////////////////////////////////////////////////////////////////
///@page fdm Finite difference Method
///@tableofcontents
/// To solve Kepler problem, we need to solve
/// \f{equation}{\label{eq:ode_2nd}
/// \frac{\mathrm{d}^2 \zeta}{\mathrm{d}t^2} =
/// \frac{1}{\zeta^3} - \frac{1}{\zeta^2}
/// \f}
/// with initial condition
/// \f{align*}{
/// \zeta(t_0) &= \zeta_0 \\
/// \zeta'(t_0) &= \zeta'_0
/// \f}
/// To solve above 2nd order ordinary differential equation
/// \latexonly \eqref{eq:ode_2nd} \endlatexonly numerically,
/// we need to approximate 2nd derivative as finite difference.
/// Suppose that the solution \f$\zeta(t)\f$ has  
/// continuous 4th order derivative
/// in the Domain \f$[t_0, t_f]\f$. then 
/// \f{equation}{\label{eq:4thtaylor}
/// \zeta(x+h) = 
/// \zeta(x) + \zeta'(x)h + \frac{1}{2!}\zeta''(x)h^2 + 
/// \frac{1}{3!}\zeta'''(x)h^3 + \frac{1}{4!}\zeta^{(4)}(\eta)h^4
/// \f}
/// for some \f$\eta(x, h) \in (t_0, t_f)\f$.
/// Using 4th order taylor approximation 
/// \latexonly \eqref{eq:4thtaylor} \endlatexonly,
/// we can get following equation
/// \f{equation}{\label{eq:finite_diff}
/// \zeta(x-h)-2\zeta(x)+\zeta(x+h) = h^2\zeta''(x) + O(h^4)
/// \f}
/// Next uniformly divide the domain \f$[t_0, t_f]\f$ 
/// into \f$ n \f$ sub intervals.
/// let \f$ x_i \f$ be the end points of the sub intervals then for 
/// \f$ 0 \leq i \leq n \f$,
/// \f{equation}{\label{eq:x}
/// x_i = t_0 + ih
///\f}
///, where \f$ h = (t_f-t_0)/n \f$.
/// Now for \f$ 0 \leq i \leq n \f$, define \f$\zeta_i \f$ as following
/// \f{equation}{\label{eq:zeta_i}
/// \zeta_i = \zeta(x_i)\f}
/// Then we can rewrite finite difference equation 
/// \latexonly \eqref{eq:finite_diff} \endlatexonly as following
/// \f{equation}{\label{eq:revised_fdiff}
/// \zeta_{i-1}-2\zeta_i+\zeta_{i+1} = h^2\zeta''_i + O(h^4) \f}
/// for \f$ 1 \leq i \leq n-1 \f$.
/// Plug this equation \latexonly \eqref{eq:revised_fdiff} \endlatexonly
/// into 2nd order ode \latexonly \eqref{eq:ode_2nd} \endlatexonly,
/// the we have following recurrence relation
/// \f{equation}{\label{eq:recurr}
/// \zeta_{i-1} -2\zeta_i + \zeta_{i+1} = 
/// h^2\left(\frac{1}{\zeta_i^3} - \frac{1}{\zeta_i^2}\right) \f}
/// In above equation \latexonly \eqref{eq:recurr} \endlatexonly,
/// we turncate, so local turncation error is 
/// \f$ O(h^4)=O(n^{-4})\f$.
/// Therefore global turncation error can be roughly estimated to
/// \f$ O(n^{-3})\f$.
/// To solve recurrance relation, we need to know both 
/// \f$\zeta_0\f$ and \f$\zeta_1\f$. 
/// However only \f$\zeta_0\f$ is explictly given by the
/// initial condition.
/// To approximate \f$\zeta_1\f$ with \f$ O(n^{-3}) \f$
/// error bound, I use 2nd order talyor expension.
/// \f{equation}{\label{eq:2ndtalyor}
/// \zeta_1 \approx \zeta_0 + \zeta'_0 h + 
/// \frac{1}{2!} \zeta''_0 h^2
/// \f}
/// \f$\zeta''_0 \f$ can be derived by 2nd order ode 
/// \latexonly \eqref{eq:ode_2nd} \endlatexonly
/// \f{equation}{\label{eq:2nd_deriv_zeta0}
/// \zeta''_0 = \frac{1}{\zeta_0^3} - \frac{1}{\zeta_0^2}\f}
/// @section complex Complexity
/// Clearly 
/// \f{equation*}{ O(n). \f}
/// @section accuracy Accuracy
/// Global turncation error is roughly estimated by 
/// \f{equation*}{O(n^{-3}). \f}
/// @section conv Convergence
/// - Initial Condition
/// \f{align*}{
/// \zeta(0) &= 0.9 \\
/// \zeta'(0) &= 0
/// \f}
/// - Initial time: 0
/// - Final time: 10
/// @image html plot_converge_zeta.png "Convergence plot"
/// @image latex plot_converge_zeta.eps "Convergence plot" width=10cm
///////////////////////////////////////////////////////////////////////