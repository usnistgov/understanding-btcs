/** @file my_funcs.hpp
 *  @brief For calculating analytical solutions
 *
 *  Used to calculate analytical solutions or break-through curves
 *
 *  @author Robert F. DeJaco
 */

#ifndef MY_FUNCS
#define MY_FUNCS

/** @brief Calculate exponential term involving bessel function
 *  
 *  Evaluates \f$e^{-z-T}I_0(2\sqrt{zT})\f$, where \f$I_0\f$ is the modified Bessel function of zeroth order. 
 *
 *  @param T Scaled value for time
 *  @param z Scaled value for distance
 *  @return value of expression 
 */
double exp_I0(double T, double z);

/** @brief Integrate expI0 from T_jm1 to T_j
 *  
 *  Evaluates \f[\int\limits_{T_{j-1}}^{T_j} e^{-z-t}I_0\left(2\sqrt{zt}\right)\, dt\f]
 *
 *  @see exp_I0
 *  @param T_jm1 scaled time at previous index j-1
 *  @param T_j scaled time at current/next index j
 *  @param z scaled distance
 *  @return value of expression
 */
double integrate(double T_jm1, double T_j, double z);

/** @brief Updates solution from previous time step
 *
 *  @param c_j new solution at time j
 *  @param a_jm1 old solution at time j-1
 *  @param T_j current/new value of time
 *  @param T_jm1 previous value of time 
 *  @param X distance along column
 *  @param n last index of distance array
 *  @return void
 *
 *  \warning{not tested, probably uses wrong value of n}
 */
double update_solution(double *c_j, double *a_jm1, double T_jm1, double T_j, double *X, int n);

/** @brief Provides analytical break-through curve
 *
 *
 *  @param c fluid concentration
 *  @param T times
 *  @param X value for distance (1/epsilon)
 *  @param m number of times including 0
 *  @return void
 */
void analytical_btc(double *, double *, double X, int);

#endif