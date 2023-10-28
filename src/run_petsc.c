static char help[] = "-Run.  \n\n";


#include <petsc.h>
#include "my_petsc.h"
#include <stdbool.h>
#include "my_funcs.hpp"
#include "isotherm.h"


int main(int argc, char **args) 
{
    PetscReal tau_star=0.1, eps=1, kappa=0., theta=0.5;
    PetscInt m=3;
    PetscBool compute_analytical=PETSC_TRUE, flg, make_movie=PETSC_FALSE, print_header=PETSC_FALSE;
    PetscErrorCode ierr;
    char btc_file_name[PETSC_MAX_PATH_LEN];
    PetscViewer file_btc;

    PetscCall(PetscInitialize(&argc,&args,NULL,help));

    PetscOptionsBegin(PETSC_COMM_WORLD,"rct_","options for reaction",""); 
    ierr = PetscOptionsGetInt(NULL, NULL, "-m",&m,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-eps",&eps,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-tau_star",&tau_star,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-theta",&theta,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-kappa",&kappa,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL, NULL, "-analy",&compute_analytical,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL, NULL, "-f",btc_file_name,sizeof(btc_file_name), &flg); CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL, NULL, "-movie",&make_movie,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL, NULL, "-pheader",&print_header,NULL); CHKERRQ(ierr);
    PetscOptionsEnd();

    // PetscPrintf(PETSC_COMM_WORLD, "%s\n", btc_file_name);

    ISO langmuir;
    initialize_langmuir_dimensionless(&langmuir, kappa);
    
    Vec btc;
    PetscReal *array_btc;
    VecCreate(PETSC_COMM_WORLD, &btc);
    VecSetSizes(btc, PETSC_DECIDE, m);
    VecSetUp(btc);
    VecSet(btc, 1.0);
    VecGetArray(btc, &array_btc);
    // VecView(btc, PETSC_VIEWER_STDOUT_WORLD);
    
    run(argc, args, m, eps, tau_star, theta, &langmuir, compute_analytical, array_btc, make_movie, print_header);

    if (flg)
    {
        // PetscPrintf(PETSC_COMM_WORLD, "printing file\n");
        PetscViewerCreate(PETSC_COMM_WORLD, &file_btc);
        PetscViewerSetType(file_btc, PETSCVIEWERASCII);
        PetscViewerFileSetMode(file_btc, FILE_MODE_WRITE);
        PetscViewerFileSetName(file_btc, btc_file_name);
        VecView(btc, file_btc);
        PetscViewerDestroy(&file_btc);
    }
    VecRestoreArray(btc, &array_btc);

    freeISOParams(&langmuir);
    PetscCall(PetscFinalize());

    return 0;
}
