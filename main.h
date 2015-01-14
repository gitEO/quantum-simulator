#include <math.h>
#include <iostream> 
#include <stdlib.h>
#include <complex>
#include <fstream>
#include <blitz/array.h>
#include <time.h>



using namespace blitz;



void enPartikkel( int i, int j, double E, double dt);
void toPartikkel( int i, int j, int k, int l, double V, double dt);
void Utest(int phaseMultiplier);
void UHubbard(int phaseMultiplier); 
