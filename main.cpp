#include "newmanziff.h"

int main(int argc, char *argv[]) {
  if (argc != 8) {
    std::cout << "Usage: ./exe.out forces_file output_file npoints timeframe max_NN Ninstances local_isostatic_condition" << "\n";
    std::cout << "npoints is the number of lattice sites or particles\n";
    std::cout << "set timeframe to -1 if you want the program to find the last frame\n";
    std::cout << "max_NN is the maximum number of nearest neighbors any particle has or could have\n";
    std::cout << "local_isostatic_condition: 0 for keep rattlers; 3 frictionless, 2 sliding, 1 sliding+roll+twist\n";
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
		int Nruns = atoi(argv[6]);
		int local_isostatic_condition = atoi(argv[7]);
		//int Nruns = 1000000;
        
        int err;
        newmanziff NZ;
        std::cout << "read\n";
        err = NZ.get_NN( forfile, max_NN, N, nstep, Nruns, local_isostatic_condition );
        std::cout << "loop\n";
        if ( err >= 0) {
          // loop over 1million iteration
          for (int i=0; i<Nruns; i++) {
            if (i % 1000 == 0) std::cout << "iteration " << i << "\n";
            NZ.permutation( N );
            NZ.percolate( N, max_NN, i );
          }
        }
        std::cout << "write\n";
        NZ.write( outfile, N, Nruns ); // average and std and write
	}
}
