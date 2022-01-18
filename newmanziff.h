#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include "string.h"
#include <string>
#include <math.h> 
#include <vector>

class newmanziff {

 public:
  int get_NN( const std::string&, const int&, const int&, int, const int&, const float&);
  void remove_rattlers(std::vector<int>&, std::vector<int>&, std::vector<int>&, const float&);
  void permutation( const int& );
  int findroot( int );
  void percolate( const int&, const int&, const int& );
  void write( const std::string&, const int&, const int& );

 protected:
  int **nn;
  int *ptr;
  int *order;

  int **max_cluster;
  double **avg_cluster;
  int **S;
};
