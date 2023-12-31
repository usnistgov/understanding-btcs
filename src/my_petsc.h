/** @file my_petsc.h
 *
 *  @brief Functions and structs for solving PDE
 *
 *  @author Robert F. DeJaco
 *
 */

#ifndef HEADER_FILE
#define HEADER_FILE

#include <stdbool.h>
#include <petsc.h>
#include "isotherm.h"

/** @struct VecAndArray
 *
 *  Used for converting between vector and array types
 *
 */
typedef struct {
    Vec vec; /*!< vector part*/
    PetscReal *arr; /*!< array part*/
} VecAndArray;


/** @struct AppCtx 
 *  @brief Useful parameters for solving, Application context
 */
typedef struct {
    PetscInt j; /*!< time step*/
    PetscInt m; /*!< total number of steps*/
    PetscReal errnorm_all; /*!< infinity norm*/
    PetscReal eps; /*!< value of \f$\varepsilon\f$ */
    ISO isotherm; /*!< structure for adsorption isotherm */
    PetscReal theta; /*!< fraction of BFD/FFD */
    PetscReal dx; /*!< grid spacing */
    PetscReal dtau;  /*!< time step*/
    PetscReal *btc; /*!< array for break-through curve*/
    PetscReal tau_star; /*!< total simulation time \f$\tau_{\star}\f$ */
    VecAndArray f_jm1; /*!< Isotherm function at previous time step*/
    VecAndArray S_jm1; /*!< Solid concentration at previous time step*/
    bool compute_analytical; /*!< Whether or not to compute analytical solution*/
    bool make_movie; /*!< whether or not to make movie*/
    bool save_btc; /*!< flag to save btc*/
    bool print_header; /*!< whether to print header when displaying results of simulation*/
    PetscViewer     viewer; /*!< petsc viewer*/
    PetscInt movie_step_interval; /*!< step interval for saving movie frames*/
    PetscReal omega; /*!< closeness to steady state at end of simulation*/
    PetscReal b; /*!< scale factor for beta */
    PetscReal K; /*!< scale factor for kappa */
    PetscReal beta; /*!< constant calculated */
    PetscReal kappa; /*!< constant calculated*/
    PetscReal dtau_est; /*!< estimated dtau for setting it when tau_star unknown */
} AppCtx;

/** @brief print user params associated with application
 */
void print_user_params(AppCtx u);


/** @struct Variables_j
 *
 *  @brief Analytical and numerical variables at time j (and jm1 if analytical)
 *
 *  \todo Change c_j and c_jm1 to a_j and a_jm1 to reflect notation change in paper
 */
typedef struct {
    VecAndArray W_j;  /*!< Numerical values of fluid concentration at time \f$j\f$*/
    VecAndArray c_j;  /*!< Analytical values of fluid concentration at time \f$j\f$*/
    VecAndArray c_jm1;  /*!< Analytical values of fluid concentration at time \f$j-1\f$*/
    Vec error; /*!< vector comprising errors (presumably b/t analytical and numerical)*/
} Variables_j;

/** @brief Initialize Struct containing both vector and array
 *  @see VecAndArray
 *
 */
extern PetscErrorCode InitializeVecAndArray(
    DM da /**< [in] manages an abstract grid object and its interactions with the algebraic solvers*/, 
    VecAndArray *x /**< Vector and array to initialize*/);

/** @brief Destroy Struct containing both vector and array
 *  @see VecAndArray
 *
 */
extern PetscErrorCode DestroyVecAndArray(
    DM da /**< [in] manages an abstract grid object and its interactions with the algebraic solvers*/, 
    VecAndArray *x /**< Vector and array to destroy*/);

/** @brief Perform one time step
 *
 */
extern PetscErrorCode one_time_step(
    SNES snes /**< [in] Solver */, 
    DMDALocalInfo *p_info /**< [in] Information about some parameters needed*/, 
    AppCtx *p_user /**< [in,out] Application context*/, 
    Variables_j *p_vars /**< [out] Variables at time \f$j\f$, updated for next time step*/
    );

/** @brief Perform calculations at all time steps
 */
extern PetscErrorCode simulate(
    DM da /**! Abstract grid object*/,
    SNES snes /**! solver*/, 
    DMDALocalInfo *p_info /**! local info on grid*/, 
    AppCtx *p_user /**< User Application context*/);

/** @brief Main function called
 */
extern PetscErrorCode run(
    int argc /**< [in] number of command arguments*/, 
    char **args /**< [in] command arguments*/,
    int m /**< [in] number of time steps*/, 
    double eps /**< [in] value of epsilon \f$\varepsilon\f$*/, 
    double tau_star /**< [in] total simulation time \f$\tau_\star\f$*/, 
    double theta /**< [in] fraction of bfd/ffd*/, 
    ISO *isotherm /**< [in] adsorption isotherm struct. Needs to be allocated/initialized beforehand*/,
    bool compute_analytical /**< [in] whether or not to compute analytical solution*/, 
    double *btc /**< [in,out] values of break-through curve*/, 
    bool make_movie /**< [in] whether or not to make movie*/,
    bool print_header /**< [in] whether or not to print header*/
    );

/** @brief Updates analytical solution
 *
 *  Updates analytical solution for time \f$j\f$ based off of info for \f$j-1\f$.
 *
 */
extern PetscErrorCode UpdateAnalytical(
    DMDALocalInfo *info /**< [in] information about local grid*/, 
    PetscReal *c_j /**< [in,out] new solution computed */,
    PetscReal *c_jm1 /**< [in] old solution*/, 
    AppCtx *user /**< [in] Application context */);

/** @brief Update solid concentration 
 */
extern PetscErrorCode UpdateSolid(
    DMDALocalInfo *info /**< [in] information about grid*/, 
    PetscReal *W_j /**< [in] Numerical solution of fluid concentration*/, 
    AppCtx *user /**< [in] application context*/);

extern PetscErrorCode UpdateIsotherm(
    DMDALocalInfo *info /**< [in] local information about grid*/, 
    PetscReal *W /**< [in] Solution of PDE at some value of time*/, 
    AppCtx *user /**< [in,out] application context. new values of isotherm stored here*/);


/** @brief evaluate value of solid concentration
 *  @return value of solid concentration \f$S_{i}^j\f$
 */
extern PetscReal eval_S_i_j(
    PetscReal f_j /**< value of adsorption isotherm at i and j*/, 
    PetscReal f_jm1 /**< value of adsorption isotherm at i but j - 1*/, 
    PetscReal exp_mdT /**< value for \f$e^{-\Delta \tau/\varepsilon}\f$ */, 
    PetscReal S_jm1 /**< Value for solid concentration at i but j - 1*/, 
    PetscReal dT /**< value for \f$\Delta \tau / \varepsilon\f$ */);

/** @brief evaluates gradient in solid concentration wrt fluid
 *  @brief returns gradient, \f$\Delta \tau  F^\prime / 2 / \varepsilon\f$
 */
extern PetscReal eval_dS_i_j(
    PetscReal df_i_j /**< [in] gradient in isotherm*/, 
    PetscReal dT /**< [in] scaled time step, likely \f$\Delta \tau / \varepsilon\f$ */);


/** @brief does theta rule for discretization in space
 *  @return value of expression averaged wrt theta, \f$v_i(1-\theta) + v_{i-1}\theta\f$ for some \f$v_i, v_{i-1}\f$
 */
extern PetscReal theta_rule(
    PetscReal val_i /**< value at spatial index i*/, 
    PetscReal val_im1 /**< value at spatial index i - 1*/, 
    PetscReal theta /**< value for \f$0 \le \theta \le 1\f$*/);

/** @brief residuals of nonlinear equations
 *  @return void
 */
extern PetscErrorCode FormFunctionLocal(
    DMDALocalInfo *info /**< Abstract object with grid info*/, 
    PetscReal *W_j /**< Solution, used to evaluate residuals*/, 
    PetscReal *FF /**< residuals or form function*/, 
    AppCtx *user /**< application context*/);

/** @brief jacobian of nonlinear equations
 *  @return void
 */
extern PetscErrorCode FormJacobianLocal(
    DMDALocalInfo *info /**< Abstract object with grid info*/, 
    PetscReal *W_j /**< Solution, used to evaluate residuals*/, 
    Mat J /**< not used?*/, 
    Mat P /**< jacobian */, 
    AppCtx *user /**< application context*/);


/** @brief calculate \f$\omega\f$
 *
 *  Parameter used to determine if simulation is at steady state
 *  @return void
 */
PetscErrorCode CalculateOmega(
    DMDALocalInfo *info /**< grid object*/, 
    PetscReal *W /**< [in] solution*/, 
    AppCtx *user /**< [in,out] stores attribute omega*/);

#endif // !HEADER_FILE