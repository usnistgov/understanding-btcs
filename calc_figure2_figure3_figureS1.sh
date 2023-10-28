export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/src

cd src && make run_petsc.exe && cd ..

# rarefaction
for e in 0.02 0.04 0.08 0.16 0.32 0.64
do 
    echo "e = " ${e}
    ./src/run_petsc.exe -analy False -m 50000 -tau_star 4.0 -da_grid_x 5000 -eps $e -kappa -0.5 -f out/rarefaction-${e}.dat
done

# shock
for e in 0.02 0.04 0.08 0.16 0.32 0.64
do 
    echo "e = " ${e}
    ./src/run_petsc.exe -analy False -m 25000 -tau_star 2.0 -da_grid_x 5000 -eps $e -kappa 1. -f out/shock-${e}.dat
done

# rarefaction extra
for e in 0.02 0.04 0.08 0.16 0.32 0.64
do 
    for k in -0.1 -0.2 -0.4 -0.8
    do
        echo "e = " ${e}
        ./src/run_petsc.exe -analy False -m 50000 -tau_star 4.0 -da_grid_x 5000 -eps $e -kappa ${k} -f out/rarefaction-${k}-${e}.dat
    done
done