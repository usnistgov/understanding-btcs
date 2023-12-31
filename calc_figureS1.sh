export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/src

cd src && make run_petsc.exe && cd ..

file=out/spatial.csv
# spatial
e=0.5
n=5
./src/run_petsc.exe -m 2049 -tau_star 2 -da_grid_x $n -eps $e -pheader True > ${file}
for n in 9 17 33 65
do 
    ./src/run_petsc.exe -m 2049 -tau_star 2 -da_grid_x $n -eps $e >> ${file}
done

# temporal
file=out/temporal1.csv
e=0.5  # 1/2
n=$((2**(6 + 1) + 1))
echo $n
m=9
./src/run_petsc.exe -m $m -tau_star 2 -da_grid_x $n -eps $e -pheader True > ${file}
for m in 17 33 65 129
do 
    ./src/run_petsc.exe -m $m -tau_star 2 -da_grid_x $n -eps $e >> ${file}
done

file=out/temporal2.csv
e=0.125  # 1/8
n=$((2**(8 + 1) + 1))
echo $n
m=33
./src/run_petsc.exe -m $m -tau_star 2 -da_grid_x $n -eps $e -pheader True > ${file}
for m in 65 129 257 513
do 
    ./src/run_petsc.exe -m $m -tau_star 2 -da_grid_x $n -eps $e >> ${file}
done

file=out/temporal3.csv
e=0.03125  # 1/32
n=$((2**(10 + 1) + 1))
echo $n
m=129
./src/run_petsc.exe -m $m -tau_star 2 -da_grid_x $n -eps $e -pheader True > ${file}
for m in 257 513 1025 2049
do 
    ./src/run_petsc.exe -m $m -tau_star 2 -da_grid_x $n -eps $e >> ${file}
done