#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include "string.h"
#include <string>

class newmanziff {

 public:
  int get_NN( const std::string&, const int&, const int&, int);
  void permutation( const int& );
  int findroot( int );
  void percolate( const std::string&, const int&, const int& );

 protected:
  int **nn;
  int *ptr;
  int *order;
};
