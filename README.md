# Compile
g++ -stdlib=libc++ -std=c++11 -O3 main.cpp newmanziff.cpp -o nz.x
./nz.x mu0.forces mu0.nz 10000 -1 12 10 3

iso=( 3 2 2 1 1 ); i=0; for f in mu*forces; do ./nz.x $f nz_$f 10000 -1 12 10000 ${iso[i]}; let i=$i+1; done
