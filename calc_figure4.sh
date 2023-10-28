export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/src

cd src && make run_petsc.exe && cd ..

# shock
e="0.1666666666666"
for k in "36.0" "6.0" "1.0"
do 
    mkdir -p movie-${k}
    ./src/run_petsc.exe -analy False -m $((2**13 + 1)) -tau_star 0.6 -da_grid_x $((2**13 + 1)) -eps $e -kappa ${k} -movie True
    mv movie-*.txt movie-${k}/
done