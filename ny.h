#include <math.h>
#include <iostream> 
#include <stdlib.h>
#include <complex>
#include <fstream>
#include <blitz/array.h>
#include <time.h>



using namespace blitz;



extern char* vFile;
extern char* mFile;
extern char* simFile;

extern double pi;



//Uavhengig av Hamiltonoperator
extern int operationsPerformed;/*We add one for each gate simulated*/
extern int dimension;/*Dimension of state, 2^(noOfQubits)*/
extern int workDim;/*Dimension of workQubit space*/
extern int simDim;/*Dimension of stateQubit space*/
extern int noOfQubits;/*Total number of qubits*/
extern int workQubits;/*Number of workqubits*/
extern int simQubits;/*Number of qubits used in the direct simulation*/
extern int firstSimQubit;/*Here is the first non-work qubit, before this
		    they all are work, after and including we have
		    simulation qubits.*/


//Simuleringsinformasjon
extern double deltaT; 
extern double Emax;
extern int noOfTimesteps;
extern int noOfMeasurements;
extern int noOfRuns;


//Hamiltonfunksjonsparametrene, åpenbart avhengig av problem
//pairing
extern double g; //koblingskonstanten
extern double dd; //levelspacinga

//Heisenberg
extern double h;



// andre ting vi trenger

// sparse matrix U;
// vector state;
// vector tempState;





extern Array<complex<double>,2>* Umatrix;
extern Array<complex<double>,2>* helpMatrix;
extern Array<complex<double>,1>* state;
extern Array<complex<double>,1>* helpVector1;
extern Array<complex<double>,1>* helpVector2;


/*An array with probabilities of eigenstates, and one with the eigenstates*/
extern Array<double,1>* probability;
extern Array<double,1>* totProb;
extern Array<double,1>* probs; /*Measured probabilities*/
extern Array<complex<double>,1>** eigenstates;



/*An array containing the measured eigenvalues*/
extern Array<double,1>* measured;

/*An array with the number of times each eigenvalue is measured, to calculate
  the statistical probability of each eigenvalue used in the simulation of 
  measurements*/
extern Array<int,1>* timesMeasured;


void readData();
void writeToFile(const char* filename, Array<double,1>* in);
void writeToFile(const char* filename, Array<int,1>* in);
void writeToFile(const char* filename, Array<complex<double>,1>* in);
Array<complex<double>,2> mult(Array<complex<double>,2> A, 
			      Array<complex<double>,2> B);
Array<complex<double>,2> outerProduct(Array<complex<double>,2> A, 
				      Array<complex<double>,2> B);
double innerProduct(Array<complex<double>,1> A);
void singleQubitGateNoMatrix(Array<complex<double>,2> unitary, int k);
void conditionalTwoQubitGateNoMatrix(Array<complex<double>,2> unitary, int k, int l);
Array<complex<double>,2>  singleQubitGate(Array<complex<double>,2> UU, 
					  int qubitOperatedOn);
Array<complex<double>,2> gCU(Array<complex<double>,2> UU, int c, int u);
Array<complex<double>,2> gUC(Array<complex<double>,2> UU, int u, 
			     int c);
void ZZ ( int i, int j, double a);
void UPairing( int phaseMultiplier);
void UpairingNonDegenerate(int phaseMultiplier);
void UHeisenberg(int phaseMultiplier);
void CU(int c);
void inverseFourierTransform();
void findingProbability();
void measurement();
void writeOutEnergies();
void initializeState(); 
int main( int argc, char* argv[] ); 
void hoved();
