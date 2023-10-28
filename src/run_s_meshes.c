static char help[] = "-Run.  \n\n";


#include <petsc.h>
#include "my_petsc.h"
#include <stdbool.h>
#include "my_funcs.hpp"
#include "isotherm.h"


int main(int argc, char **args) 
{

    PetscErrorCode ierr;
    PetscInt N=8 /*n - 1 min*/, M=2 /* m - 1 min*/, s=2;
    PetscReal theta=0.5, eps=1, tau_star=0.1, kappa=0;
    PetscInt i; // for counting

    PetscCall(PetscInitialize(&argc,&args,NULL,help));
    PetscOptionsBegin(PETSC_COMM_WORLD,"rct_","options for reaction",""); 
    ierr = PetscOptionsGetInt(NULL, NULL, "-M",&M,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-N",&N,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-s",&s,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-eps",&eps,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-tau_star",&tau_star,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-theta",&theta,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-kappa",&kappa,NULL); CHKERRQ(ierr);
    // ierr = PetscOptionsGetBool(NULL, NULL, "-pheader",&print_header,NULL); CHKERRQ(ierr);
    PetscOptionsEnd();

    // printf("Performing %d meshes refining by 2 starting at N=%d\n", s, N);

    Variables_j     *varss;
    AppCtx        *users;
    DM            *das;
    SNES          *sness;
    DMDALocalInfo *infos;
    PetscReal *errors_infty;
    users = (AppCtx *) malloc(sizeof(AppCtx)*s);
    das = (DM *) malloc(sizeof(DM)*s);
    sness = (SNES *) malloc(sizeof(SNES)*s);
    infos = (DMDALocalInfo *) malloc(sizeof(DMDALocalInfo)*s);
    varss = (Variables_j *) malloc(sizeof(Variables_j)*s);
    errors_infty = (PetscReal *) malloc(sizeof(PetscReal)*(s-1));

    ISO langmuir;
    initialize_langmuir_dimensionless(&langmuir, kappa);

    for (i = 0; i < s; i++)
    {
        users[i].make_movie = PETSC_FALSE;
        users[i].compute_analytical = PETSC_FALSE;
        users[i].save_btc = PETSC_FALSE;
        users[i].tau_star = tau_star;
        users[i].m = M + 1;
        users[i].eps = eps;
        users[i].theta = theta;
        users[i].isotherm = langmuir;
        users[i].dtau = tau_star / M;
        users[i].errnorm_all = 0.;

        PetscCall(DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,N*PetscPowInt(2, i) + 1,1,1,NULL,&(das[i])));
        PetscCall(DMSetUp(das[i]));
        PetscCall(DMSetApplicationContext(das[i],&(users[i])));

        PetscCall(SNESCreate(PETSC_COMM_WORLD,&(sness[i])));
        SNESSetTolerances(sness[i], 1e-8, 1e-8, 1e-8, 100000, 100000);
        SNESSetType(sness[i], SNESNEWTONLS);
        PetscCall(SNESSetDM(sness[i],das[i]));
        PetscCall(DMDASNESSetFunctionLocal(das[i],INSERT_VALUES,(DMDASNESFunction)FormFunctionLocal,&(users[i])));
        PetscCall(DMDASNESSetJacobianLocal(das[i],(DMDASNESJacobian)FormJacobianLocal,&(users[i])));
        PetscCall(SNESSetFromOptions(sness[i]));

        PetscCall(InitializeVecAndArray(das[i], &(users[i].f_jm1)));
        PetscCall(InitializeVecAndArray(das[i], &(users[i].S_jm1)));
        PetscCall(DMDAGetLocalInfo(das[i],&(infos[i])));
        users[i].dx = 1.0 / (infos[i].mx - 1);

        PetscCall(InitializeVecAndArray(das[i], &(varss[i].W_j)));
        PetscCall(VecDuplicate(varss[i].W_j.vec,&(varss[i].error)));

        // set initial guess
        VecSet(varss[i].W_j.vec, 1.0);
        VecSet(users[i].S_jm1.vec, 0.0);
        VecSet(users[i].f_jm1.vec, 0.0);
        if (i < (s - 1))   errors_infty[i] = 0;
    }

    PetscReal err_ij;
    PetscInt i_1, i_2, i_last;
    PetscInt j;
    // time stepping
    for (j = 0; j < M + 1; j++)
    {
        /* simulate each*/
        for (i = 0; i < s; i++)
        {
            users[i].j = j;
            one_time_step(sness[i], &(infos[i]), &(users[i]), &(varss[i]));
        }

        /* update error*/
        for (i = 0; i < s-1; i++)
        {
            for (i_1=infos[i].xs; i_1 < infos[i].xs + infos[i].xm; i_1++)
            {
                i_2 = i_1 * 2;
                err_ij = PetscAbs(varss[i].W_j.arr[i_1] - varss[i+1].W_j.arr[i_2]);
                if (err_ij > errors_infty[i])   errors_infty[i] = err_ij;

                i_last = i_1 * PetscPowInt(2, s-1 - i);
                err_ij = PetscAbs(varss[i].W_j.arr[i_1] - varss[s - 1].W_j.arr[i_last]);
                if (err_ij > users[i].errnorm_all)   users[i].errnorm_all = err_ij;
            }
        }
    }

    PetscPrintf(PETSC_COMM_WORLD, "n,r,m,eps,tau_star,theta,kappa,error_inf,error_star,W_N^M\n");
    for (i = 0; i < s-1; i++)
    {
        PetscPrintf(PETSC_COMM_WORLD, "%d,%d,%d,%g,%g,%g,%g,%g,%g,%g\n", 
                    infos[i].mx, infos[i+1].mx, users[i].m, 
                    users[i].eps, users[i].tau_star, users[i].theta, users[i].isotherm.params[0], 
                    errors_infty[i], users[i].errnorm_all,
                    varss[i].W_j.arr[N*PetscPowInt(2, i)]);
    }

    for (i = 0; i < s; i++)
    {
        PetscCall(DestroyVecAndArray(das[i], &(varss[i].W_j)));
        PetscCall(VecDestroy(&varss[i].error));
        PetscCall(DestroyVecAndArray(das[i], &(users[i].f_jm1)));
        PetscCall(DestroyVecAndArray(das[i], &(users[i].S_jm1)));
        PetscCall(SNESDestroy(&(sness[i])));
        PetscCall(DMDestroy(&(das[i])));
    }
    
    freeISOParams(&langmuir);
    free(users);
    free(das);
    free(sness);
    free(infos);
    free(varss);
    free(errors_infty);
    PetscCall(PetscFinalize());
    return 0;
}
