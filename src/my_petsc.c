
#include <petsc.h>
#include <stdbool.h>
#include "my_funcs.hpp"
#include "my_petsc.h"
#include "isotherm.h"
#include <stdio.h>


PetscErrorCode InitializeVecAndArray(DM da, VecAndArray *x)
{
    PetscCall(DMCreateGlobalVector(da,&(x->vec)));
    PetscCall(DMDAVecGetArray(da,x->vec,&(x->arr)));
    return 0;
}
PetscErrorCode DestroyVecAndArray(DM da, VecAndArray *x)
{
    PetscCall(DMDAVecRestoreArray(da,x->vec,&(x->arr)));
    PetscCall(VecDestroy(&(x->vec)));
    return 0;
}

PetscErrorCode one_time_step(SNES snes, DMDALocalInfo *p_info, AppCtx *p_user, Variables_j *p_vars)
{
    PetscReal errnorm_j;
    char f_name[PETSC_MAX_PATH_LEN];

    PetscCall(SNESSolve(snes,NULL,(p_vars -> W_j).vec));
    if (p_user-> j > 0) PetscCall(UpdateSolid(p_info, (p_vars -> W_j).arr, p_user)); // note do this first, S_ij depends on f_ijm1
    PetscCall(UpdateIsotherm(p_info, (p_vars -> W_j).arr, p_user));

    if (p_user -> compute_analytical)
    {
        PetscCall(UpdateAnalytical(p_info, (p_vars -> c_j).arr, (p_vars -> c_jm1).arr, p_user));
        PetscCall(VecCopy((p_vars -> W_j).vec, p_vars -> error));
        PetscCall(VecAXPY(p_vars -> error,-1.0,(p_vars -> c_j).vec));    // W_j <- W_j - c_j
        PetscCall(VecNorm(p_vars -> error,NORM_INFINITY,&errnorm_j));
        if (errnorm_j > (p_user -> errnorm_all))   (p_user -> errnorm_all) = errnorm_j;
        PetscCall(VecCopy((p_vars -> c_j).vec,(p_vars->c_jm1).vec));
    }
    if (p_user -> make_movie && (p_user -> j % p_user -> movie_step_interval == 0))
    {
        sprintf(f_name, "movie-%d.txt", p_user -> j);
        PetscViewerFileSetName(p_user -> viewer, f_name);
        VecView((p_vars->W_j).vec, p_user -> viewer);
    }
    if (p_user -> save_btc)
    {
        (p_user -> btc)[p_user -> j] = (p_vars->W_j).arr[p_info -> mx - 1];
    }
    return 0;
}

PetscErrorCode simulate(DM da, SNES   snes, DMDALocalInfo *p_info, AppCtx *p_user)
{
    Variables_j     vars;

    if (p_user -> make_movie)
    {
        p_user -> movie_step_interval = (p_user -> m)/10;
        PetscPrintf(PETSC_COMM_WORLD, "Making movie after each %dth step\n", p_user -> movie_step_interval);
        PetscViewerCreate(PETSC_COMM_WORLD, &(p_user -> viewer));
        PetscViewerSetType(p_user -> viewer, PETSCVIEWERASCII);
        PetscViewerFileSetMode(p_user -> viewer, FILE_MODE_WRITE);
    }
    PetscCall(DMDAGetLocalInfo(da,p_info));
    p_user -> dx = 1.0 / (p_info -> mx - 1); /* TODO: should set this before simulate*/
    PetscCall(InitializeVecAndArray(da, &(vars.W_j)));
    PetscCall(InitializeVecAndArray(da, &(vars.c_j)));
    PetscCall(InitializeVecAndArray(da, &(vars.c_jm1)));
    PetscCall(VecDuplicate(vars.W_j.vec,&vars.error));

    // set initial guess
    VecSet(vars.W_j.vec, 1.0);
    VecSet((p_user -> S_jm1).vec, 0.0);
    VecSet((p_user -> f_jm1).vec, 0.0);

    // time stepping
    for (p_user -> j = 0; p_user -> j < p_user -> m; (p_user -> j)++)
    {
        one_time_step(snes, p_info, p_user, &vars);
    }
    // calculate omega
    PetscCall(CalculateOmega(p_info, vars.W_j.arr, p_user));

    if (p_user->print_header)
    {
        PetscPrintf(PETSC_COMM_WORLD, "n,m,eps,tau_star,theta,error_inf,omega,W_N^M");
        for (int ip = 0; ip < (p_user->isotherm).num_params; ip++)
        {
            PetscPrintf(PETSC_COMM_WORLD, ",isoparam%d", ip);
        }
        PetscPrintf(PETSC_COMM_WORLD, "\n");
    }
    PetscPrintf(PETSC_COMM_WORLD, "%d,%d,%g,%g,%g,%g,%g,%g", p_info -> mx, p_user -> m, p_user->eps, p_user->tau_star, p_user->theta, p_user -> errnorm_all, p_user -> omega, vars.W_j.arr[p_info->mx - 1]);
    for (int ip = 0; ip < (p_user->isotherm).num_params; ip++)
    {
        PetscPrintf(PETSC_COMM_WORLD, ",%g", (p_user->isotherm).params[ip]);
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");

    PetscCall(DestroyVecAndArray(da, &(vars.W_j)));
    PetscCall(DestroyVecAndArray(da, &(vars.c_j)));
    PetscCall(DestroyVecAndArray(da, &(vars.c_jm1)));
    PetscCall(VecDestroy(&vars.error));
    if (p_user->make_movie)
    {
        PetscViewerDestroy(&(p_user -> viewer));
    }
    return 0;
}

PetscErrorCode run(
    int argc, char **args, int m, double eps, double tau_star, 
    double theta, ISO *isotherm, bool compute_analytical, double *btc, bool make_movie,
    bool print_header)
{
    DM            da;
    SNES          snes;
    DMDALocalInfo info;
    AppCtx        user;

    // set from options
    PetscCall(DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,10,1,1,NULL,&da));
    PetscCall(DMSetFromOptions(da));
    PetscCall(DMSetUp(da));
    PetscCall(DMSetApplicationContext(da,&user));

    PetscCall(SNESCreate(PETSC_COMM_WORLD,&snes));
    PetscCall(SNESSetDM(snes,da));
    PetscCall(DMDASNESSetFunctionLocal(da,INSERT_VALUES,(DMDASNESFunction)FormFunctionLocal,&user));
    PetscCall(DMDASNESSetJacobianLocal(da,(DMDASNESJacobian)FormJacobianLocal,&user));
    SNESSetTolerances(snes, 1e-8, 1e-8, 1e-8, 100000, 100000);
    SNESSetType(snes, SNESNEWTONLS);
    PetscCall(SNESSetFromOptions(snes));

    PetscCall(InitializeVecAndArray(da, &(user.f_jm1)));
    PetscCall(InitializeVecAndArray(da, &(user.S_jm1)));

    user.dtau = tau_star / (m - 1);
    user.tau_star = tau_star;
    user.errnorm_all = 0.;
    user.m = m;
    user.eps = eps;
    user.isotherm = (*isotherm);
    user.theta = theta;
    user.compute_analytical = compute_analytical;
    user.make_movie = make_movie;
    user.save_btc = PETSC_TRUE;
    user.print_header = print_header;
    user.btc = btc;
    user.omega = 0.;
    simulate(da, snes, &info, &user);


    // destroy
    PetscCall(DestroyVecAndArray(da, &(user.f_jm1)));
    PetscCall(DestroyVecAndArray(da, &(user.S_jm1)));
    PetscCall(SNESDestroy(&snes));
    PetscCall(DMDestroy(&da));

    return 0;
} // end of run()

PetscErrorCode UpdateAnalytical(
  DMDALocalInfo *info, PetscReal *c_j, 
  PetscReal *c_jm1, AppCtx *user
) {
    PetscInt   i;
    double X, T_j, T_jm1, dX;
    dX = (user->dx)/(user->eps);
    // PetscPrintf(PETSC_COMM_WORLD, "analytical, j = %d\n", user->j);
    if (user -> j == 0)
    {
        for (i=info->xs; i<info->xs+info->xm; i++)
        {
            X = (double) i*dX;
            c_j[i] = exp_I0(0., X);
        }
    }
    else
    {
        T_j = (double) (user -> j) * (user -> dtau) / (user -> eps);
        T_jm1 = (double) (user -> j - 1) * (user -> dtau) / (user -> eps) ;
        PetscReal tmp;

        for (i=info->xs; i<info->xs+info->xm; i++)
        {
            X = (double) i*dX;
            tmp = (PetscReal) exp_I0(T_j, X) - exp_I0(T_jm1, X) + integrate(T_jm1, T_j, X);
            c_j[i] = c_jm1[i] + tmp;
        }
    }

    return 0;
}

PetscErrorCode UpdateSolid(DMDALocalInfo *info, PetscReal *W_j, AppCtx *user){
    PetscInt   i;
    PetscReal  *f_jm1 = (user->f_jm1).arr, *S_jm1 = (user->S_jm1).arr;
    PetscReal dT = (user->dtau)/(user->eps);
    PetscReal exp_mdT = PetscExpReal(-(user->dtau)/(user->eps));
    PetscReal S_i_j;
    PetscReal f_i_j;

    for (i=info->xs; i<info->xs+info->xm; i++) {
        f_i_j = (user->isotherm).F(W_j[i], (user->isotherm).num_params, (user->isotherm).params, (user->isotherm).data);
        S_i_j = eval_S_i_j(
            f_i_j, f_jm1[i], exp_mdT, S_jm1[i], dT
        );
        // PetscPrintf(PETSC_COMM_WORLD, "S_i_j for i = %d is %g\n", i, S_i_j);
        (user->S_jm1).arr[i] = S_i_j;
    }
    return 0;

}

PetscErrorCode UpdateIsotherm(DMDALocalInfo *info, PetscReal *W, AppCtx *user) {
    PetscInt   i;
    PetscReal f_i_j;

    for (i=info->xs; i<info->xs+info->xm; i++) {
        f_i_j = (user->isotherm).F(W[i], (user->isotherm).num_params, (user->isotherm).params, (user->isotherm).data);
        (user->f_jm1).arr[i] = f_i_j;
    }
    return 0;
}



PetscReal eval_S_i_j(PetscReal f_j, PetscReal f_jm1, PetscReal exp_mdT, PetscReal S_jm1, PetscReal dT){
    // PetscPrintf(PETSC_COMM_WORLD, " eval_S_i_j: f_j = %g, f_jm1 = %g, exp_mdT = %g, S_jm1 = %g, dT = %g\n", f_j, f_jm1, exp_mdT, S_jm1, dT);
    return dT/2*(f_j + exp_mdT*f_jm1) + exp_mdT*S_jm1;
}

PetscReal eval_dS_i_j(PetscReal df_i_j, PetscReal dT){
    return dT/2*df_i_j;
}

PetscReal theta_rule(PetscReal val_i, PetscReal val_im1, PetscReal theta){
    return val_i*(1.0 - theta) + val_im1*theta;
}

// //STARTFUNCTIONS
PetscErrorCode FormFunctionLocal(
  DMDALocalInfo *info, PetscReal *W_j, PetscReal *FF, AppCtx *user
  ) {
    PetscInt   i;
    PetscReal  e = user->eps, theta=user->theta, dx = user->dx;
    PetscReal  *f_jm1 = (user->f_jm1).arr, *S_jm1 = (user->S_jm1).arr;
    PetscReal dT = (user -> dtau) / (user-> eps);
    PetscReal  exp_mdT = PetscExpReal(-(user->dtau)/(user->eps));
    PetscReal f_i_j, f_im1_j, S_i_j, S_im1_j;

    for (i=info->xs; i<info->xs+info->xm; i++) {
        if (i == 0) {
            FF[i] = W_j[i] - 1.0;
        } 
        else {  // i = 1 to mx - 1
            f_i_j = (user->isotherm).F(W_j[i], (user->isotherm).num_params, (user->isotherm).params, (user->isotherm).data);
            f_im1_j = (user->isotherm).F(W_j[i-1], (user->isotherm).num_params, (user->isotherm).params, (user->isotherm).data);
            FF[i] = W_j[i] - W_j[i-1] + dx/e*theta_rule(f_i_j, f_im1_j, theta);
            if (user->j > 0){
                S_i_j = eval_S_i_j(f_i_j, f_jm1[i], exp_mdT, S_jm1[i], dT);
                S_im1_j = eval_S_i_j(f_im1_j, f_jm1[i-1], exp_mdT, S_jm1[i-1], dT);
                FF[i] -= dx/e*theta_rule(S_i_j, S_im1_j, theta);
            }
        }
    }
    return 0;
}

PetscErrorCode FormJacobianLocal(DMDALocalInfo *info, PetscReal *W_j,
                                 Mat J, Mat P, AppCtx *user) {
    PetscInt   i, col[2];
    PetscReal  e = user->eps, theta=user->theta, v[2], dx=user->dx;
    PetscReal dT = (user -> dtau) / (user-> eps);
    PetscReal df_i_j, df_im1_j, dS_i_j, dS_im1_j;

    for (i=info->xs; i<info->xs+info->xm; i++) {
        if (i == 0) {
            v[0] = 1.0;
            PetscCall(MatSetValues(P,1,&i,1,&i,v,INSERT_VALUES));
        } else {
            df_i_j = (user->isotherm).dF(W_j[i], (user->isotherm).num_params, (user->isotherm).params, (user->isotherm).data);
            df_im1_j = (user->isotherm).dF(W_j[i-1], (user->isotherm).num_params, (user->isotherm).params, (user->isotherm).data);
            col[0] = i;
            v[0] = 1.0 + df_i_j*theta*dx/e;
            col[1] = i - 1;
            v[1] = -1.0 + df_im1_j*(1-theta)*dx/e;
            if (user->j > 0) {
                dS_i_j = eval_dS_i_j(df_i_j, dT);
                dS_im1_j = eval_dS_i_j(df_im1_j, dT);
                v[0] -= dS_i_j*theta*dx/e;
                v[1] -= dS_im1_j*(1 - theta)*dx/e;
            }
            PetscCall(MatSetValues(P,1,&i,2,col,v,INSERT_VALUES));
        }
    }

    PetscCall(MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY));
    if (J != P) {
        PetscCall(MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY));
    }
    return 0;
}


PetscErrorCode CalculateOmega(DMDALocalInfo *info, PetscReal *W, AppCtx *user) {
    PetscInt   i;
    PetscReal error;

    for (i=info->xs + 1; i<info->xs+info->xm; i++) {
        error = PetscAbsReal(W[i] - W[i-1]);
        if (error > (user->omega))   (user->omega) = error;
    }
    return 0;
}


void print_user_params(AppCtx u)
{
    printf("----User Parameters--------\n");
    printf("\tTotal number of steps is %d\n", u.m);
    printf("\tTotal simulation time is %f\n", u.tau_star);
    printf("\tEpsilon is %f\n", u.eps);
    printf("\tTheta is %f\n", u.theta);
    printf("\tGrid spacing is %f\n", u.dx);
    printf("\tTime step is %f\n", u.dtau);
    printf("\tMovie step interval is %d\n", u.movie_step_interval);
}
//ENDFUNCTIONS