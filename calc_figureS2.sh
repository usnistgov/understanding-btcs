export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/src

cd src && make run_s_meshes.exe && make run_s_timesteps.exe && cd ..

theta=0.5
e=1
constant_args="-eps ${e} -theta ${theta}"

# spatial convergence
for kappa in "0" "1" 
do
    file_space="out/spatial-refinement-theta=${theta}-kappa=${kappa}-e=${e}.csv"
    echo ${file_space}
    args="-kappa $kappa ${constant_args}"
    ./src/run_s_meshes.exe -M 131072 -N 4 -s 7 -tau_star 4 ${args} > ${file_space}
done

for kappa in "8" "64" "256"
do
    file_space="out/spatial-refinement-theta=${theta}-kappa=${kappa}-e=${e}.csv"
    echo ${file_space}
    args="-kappa $kappa ${constant_args}"
    ./src/run_s_meshes.exe -M 131072 -N 16 -s 7 -tau_star 4 ${args} > ${file_space}
done

for kappa in "-0.4" "-0.8"
do
    file_space="out/spatial-refinement-theta=${theta}-kappa=${kappa}-e=${e}.csv"
    echo ${file_space}
    args="-kappa $kappa ${constant_args}"
    ./src/run_s_meshes.exe -M 262144 -N 4 -s 7 -tau_star 8 ${args} > ${file_space}
done

# temporal 
for kappa in "0" "1" "8" "64" "256"
do
    file_time="out/temporal-refinement-theta=${theta}-kappa=${kappa}-e=${e}.csv"
    echo ${file_time}
    args="-kappa $kappa ${constant_args}"
    ./src/run_s_timesteps.exe -M 8 -N 8192 -s 7 -tau_star 4 ${args} > ${file_time}
done

# temporal 
for kappa in "-0.8" "-0.4"
do
    file_time="out/temporal-refinement-theta=${theta}-kappa=${kappa}-e=${e}.csv"
    echo ${file_time}
    args="-eps $e -kappa $kappa ${constant_args}"
    ./src/run_s_timesteps.exe -M 8 -N 8192 -s 7 -tau_star 8 ${args} > ${file_time}
done
