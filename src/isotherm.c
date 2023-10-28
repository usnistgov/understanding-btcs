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


double F_li(double c, int n, double *cs, double *fs)
{
    int i;
    if (c <= cs[0])
    {
        return fs[0]/cs[0]*c;
    }
    for (i = 1; i < n; i++)
    {
        if (c <= cs[i])
        {
            return fs[i-1] + (fs[i] - fs[i - 1])/(cs[i] - cs[i - 1])*(c - cs[i-1]);
        }
    }

    return fs[n - 1];
}

double dF_li(double c, int n, double *cs, double *fs)
{
    int i;
    if (c <= cs[0])
    {
        return fs[0]/cs[0];
    }
    for (i = 1; i < n; i++)
    {
        if (c <= cs[i])
        {
            return (fs[i] - fs[i - 1])/(cs[i] - cs[i - 1]);
        }
    }

    return 0.;
}

void free_li(ISO *iso)
{
    free(iso -> params);
    free(iso -> data);
}

/* Langmuir Isotherm END */
void initialize_li(
    ISO *iso, EXP_ISO exp_iso, PetscReal c_f, PetscReal *f_f
    )
{
    int i;
    PetscReal f_i;
    iso -> num_params = exp_iso.n; /* NOTE: this is not number of params!*/

    iso -> params = (double *) malloc(sizeof(double)*(iso -> num_params));
    iso -> data = (double *) malloc(sizeof(double)*(iso -> num_params));

    /* store c and f*/
    for (i = 0; i < exp_iso.n; i++)
    {
        (iso -> params)[i] = exp_iso.c_star[i]*exp_iso.c_scale;
        (iso -> data)[i] = exp_iso.f_star[i]*exp_iso.f_scale;
    }

    /* calculate f(cf) */
    *f_f = F_li(c_f, iso -> num_params, iso -> params, iso -> data);
    f_i = F_li(0., iso -> num_params, iso -> params, iso -> data);

    if (f_i < 0.)
    {
        printf("Initial solid concentration is negative, not physical!!");
    }


    for (i = 0; i < exp_iso.n; i++)
    {
        (iso -> params)[i] /= c_f;
        (iso -> data)[i] /= (*f_f);
    }

    iso -> F = F_li;
    iso -> dF = dF_li;
}