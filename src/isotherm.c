#include <petsc.h>
#include <malloc.h>
#include "isotherm.h"


void initISOParams(ISO *isotherm)
{
    isotherm -> params = (double *) malloc(sizeof(double)*(isotherm -> num_params));
}

void freeISOParams(ISO *isotherm)
{
    free(isotherm -> params);
}

double F_langmuir_dimensionless(double c, int num_params, double *params, double *data)
{
    double kappa = params[0];
    return (1 + kappa)*c/(1 + kappa*c);
}

double dF_langmuir_dimensionless(double c, int num_params, double *params, double *data)
{
    double kappa = params[0];
    return (1 + kappa)/(1 + kappa*c)/(1 + kappa*c);
}

void initialize_langmuir_dimensionless(ISO *langmuir, double kappa)
{
    langmuir -> num_params = 1;
    initISOParams(langmuir);
    (langmuir -> params)[0] = kappa;
    langmuir -> F = F_langmuir_dimensionless;
    langmuir -> dF = dF_langmuir_dimensionless;
}