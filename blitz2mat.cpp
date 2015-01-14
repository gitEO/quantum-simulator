#include <math.h>
#include <iostream> 
#include <stdlib.h>
#include <blitz/array.h>
#include <complex>
#include <fstream>


using namespace blitz;
using namespace std;

Array<complex<double>,1> A;


void writeToFile(const char*);

void readFromFile(const char*);



int main( int argc, char* argv[] ) {
  


 cout << "Argv[1]  " << argv[1] << "  Argv[2]  " << argv[2] << endl;

 readFromFile(argv[1]);
 writeToFile(argv[2]);

}//main



void writeToFile(const char* filename) {

  ofstream ofs(filename);
  if (ofs.bad())
    {
      cerr << "Unable to write to file: " << filename << endl;
      exit(1);
    } 
  
  
  
  for(int j=0; j<A.rows(); j++){
    ofs << real(A(j)) << "+" << imag(A(j)) << "i  ";
  }//for
  ofs << endl;

}//writeToFile



void readFromFile(const char* filename) {

  ifstream ifs(filename);
  if (ifs.bad())
    {
      cerr << "Unable to open file: " << filename << endl;
      exit(1);
    }
  

  ifs >> A;
 

  //cout << "Arrays restored from file: " << A <<endl;

  
  
}//readFromFile
