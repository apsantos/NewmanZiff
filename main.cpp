#include "newmanziff.h"

int main(int argc, char *argv[]) {
  if (argc != 6) {
    std::cout << "Usage: ./exe.out forces_file output_file npoints timeframe max_NN" << "\n";
    std::cout << "npoints is the number of lattice sites or particles\n";
    std::cout << "set timeframe to -1 if you want the program to find the last frame\n";
    std::cout << "max_NN is the maximum number of nearest neighbors any particle has or could have\n";
  }
  else {
    // Parse input
        // file with the neighbors of each particle/lattice site
		std::string forfile = argv[1];
		std::string outfile = argv[2];
        // maximum nearest neighbors
		int N = atoi(argv[3]);
		int nstep = atoi(argv[4]);
		int max_NN = atoi(argv[5]);
		int Nruns = 1000;
		//int Nruns = 1000000;
        
        int err;
        newmanziff NZ;
std::cout << "get\n";
        err = NZ.get_NN( forfile, max_NN, N, nstep );
        if ( err >= 0) {
          // loop over 1million iteration
          for (int i=0; i<Nruns; i++) {
            if (i % 1000 == 0) std::cout << i << "\n";
            NZ.permutation( N );
            NZ.percolate( N, max_NN );
            NZ.statistics( N );
          }
        }
        NZ.write( outfile );
	}
}
