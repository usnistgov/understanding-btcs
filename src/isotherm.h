/** @file isotherm.h
 *
 *  @brief Functions and structs for adsorption isotherm
 *
 *  @author Robert F. DeJaco
 *
 */

#ifndef ISOTHERM 
#define ISOTHERM

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

#endif // !HEADER_FILE