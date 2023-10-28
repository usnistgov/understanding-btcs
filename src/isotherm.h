/** @file isotherm.h
 *
 *  @brief Functions and structs for adsorption isotherm
 *
 *  @author Robert F. DeJaco
 *
 */

#ifndef ISOTHERM 
#define ISOTHERM

/* Constant R in m**3*Pa/mol/K */
#define R 8.31446261815324
#include <petsc.h>

/** @brief F
 *
 *  Function pointer used for isotherm function
 *
 */
typedef double (*Fun)(double c, int num_params, double *params, double *data); 

/** @brief dF
 *
 *  Function pointer used for isotherm derivative
 *
 */
typedef double (*dFun)(double c,  int num_params, double *params, double *data);

/** @struct Isotherm
 *
 *  Used for adsorption isotherm 
 */
 typedef struct {
    Fun F; /*!< adsorption isotherm*/
    dFun dF; /*!< derivative of adsorption isotherm*/
    int num_params; /*!< number of parameters used */
    double *params; /*!< parameters fit */
    double *data; /*!< data used for splines */
 } ISO;

/** @brief allocate memory for isotherm parameters array
 */
void initISOParams(ISO *isotherm);

/** @brief free isotherm parameters array
 */
void freeISOParams(ISO *isotherm);

/** @brief evaluate adsorption isotherm 
 *  
 *  Evaluates \f[ \frac{\left(1 + \kappa\right)c}{1 + \kappa c}\f]
 *  where \f$\kappa = p[0]\f$ and \f$p\f$ is the parameters array.
 *
 *  @return value of adsorption isotherm 
 */
double F_langmuir_dimensionless(double c, int num_params, double *params, double *data);

/** @brief evaluate derivative adsorption isotherm 
 *  
 *  Evaluates \f[ \frac{1 + \kappa}{\left(1 + \kappa c\right)^2}\f]
 *  where \f$\kappa = p[0]\f$ and \f$p\f$ is the parameters array.
 *
 *  @return value of adsorption isotherm 
 */
double dF_langmuir_dimensionless(double c, int num_params, double *params, double *data);


/** @brief initialize langmuir struct
 *  
 *  Initializes langmuir struct. Sets number of params to be 1, 
 *  initializes parameters array, sets F and dF to point to F_langmuir and dF_Langmuir.
 *  
 */
void initialize_langmuir_dimensionless(ISO *langmuir, double kappa);

/** @struct isotherm adata
 */
typedef struct {
    double *c_star; /*!< scaled fluid concentration */
    double *f_star; /*!< scaled solid concentration */
    int n; /*!< number of points*/
    double c_scale; /*!< scale factor for concentration, multiplying c_scale by c_star gives units of mol/m3*/
    double f_scale; /*!< scale factor for loading, multiplying f_scale by f_star gives units of mol/m3*/
} EXP_ISO;


/** @brief isotherm function for linear interpolation 
 *
 */
double F_li(double c, int n, double *cs, double *fs);

/** @brief Derivative of isotherm function for linear interpolation 
 *
 */
double dF_li(double c, int n, double *cs, double *fs);

/** @brief free li after stored 
 */
void free_li(ISO *iso);

void initialize_li(
    ISO *iso,  /**< [out] dimensionless iso for pde*/
    EXP_ISO exp_iso, /**< [in] experimental isotherm data*/
    PetscReal c_f, /**< [in] feed concentration mol/m3 exp */
    PetscReal *f_f /**< [out] feed equil solid concentration mol/m3 exp */
);

#endif // !HEADER_FILE