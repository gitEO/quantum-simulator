#include <math.h>
#include <iostream> 
#include <stdlib.h>
#include <complex>
#include <fstream>
#include <blitz/array.h>
#include <time.h>
#include "main.h"
#include "ny.h"

using namespace blitz;

char* vFile;
char* mFile;
char* simFile;
double pi;







//Uavhengig av Hamiltonoperator
int operationsPerformed;/*We add one for each gate simulated*/
int dimension;/*Dimension of state, 2^(noOfQubits)*/
int workDim;/*Dimension of workQubit space*/
int simDim;/*Dimension of stateQubit space*/
int noOfQubits;/*Total number of qubits*/
int workQubits;/*Number of workqubits*/
int simQubits;/*Number of qubits used in the direct simulation*/
int firstSimQubit;/*Here is the first non-work qubit, before this
		    they all are work, after and including we have
		    simulation qubits.*/


//Simuleringsinformasjon
double deltaT; 
double Emax;
int noOfTimesteps;
int noOfMeasurements;
int noOfRuns;


//Hamiltonfunksjonsparametrene, åpenbart avhengig av problem
//pairing
double g; //koblingskonstanten
double dd; //levelspacinga

//Heisenberg
double h;



// andre ting vi trenger

// sparse matrix U;
// vector state;
// vector tempState;


// matrices H, I, x,y,z;
Array<complex<double>,2> x(2,2), y(2,2), z(2,2), H(2,2), I(2,2);


Array<complex<double>,2>* Umatrix;
Array<complex<double>,2>* helpMatrix;
Array<complex<double>,1>* state;
Array<complex<double>,1>* helpVector1;
Array<complex<double>,1>* helpVector2;


/*An array with probabilities of eigenstates, and one with the eigenstates*/
Array<double,1>* probability;
Array<double,1>* totProb;
Array<double,1>* probs; /*Measured probabilities*/
Array<complex<double>,1>** eigenstates;

/*An array containing the measured eigenvalues*/
Array<double,1>* measured;

/*An array with the number of times each eigenvalue is measured, to calculate
  the statistical probability of each eigenvalue used in the simulation of 
  measurements*/
Array<int,1>* timesMeasured;







/*Her starter metodene*/







void readData() {

  ifstream ifs(simFile);
  if (ifs.bad())
    {
      cerr << "Unable to open file: " << simFile << endl;
      exit(1);
    }
  
  ifs >> noOfQubits >> firstSimQubit >> dd >> g >> deltaT >> Emax >> noOfTimesteps >> 
    noOfMeasurements >> noOfRuns;
  
  cout << "Number of qubits:" << noOfQubits 
       << "\nFirst simulation qubit " << firstSimQubit 
       << "\nd  " << dd 
       << "\ng  " << g 
       << "\nTimestep " 
       << deltaT << "\nEmax "<< Emax << "\nNo. of timesteps "<< noOfTimesteps  
       << "\nNo. of measurements " << noOfMeasurements 
       << "\nNo. of runs " << noOfRuns << endl;
  
}//end readData






void writeToFile(const char* filename, Array<double,1>* in){

  ofstream ofs(filename);
  if (ofs.bad())
    {
      cerr << "Unable to write to file: " << filename << endl;
      exit(1);
    } 
  
  ofs << *in << endl;
  
}//write



void writeToFile(const char* filename, Array<int,1>* in){

  ofstream ofs(filename);
  if (ofs.bad())
    {
      cerr << "Unable to write to file: " << filename << endl;
      exit(1);
    } 
  
  ofs << *in << endl;
  
}//write


void writeToFile(const char* filename, Array<complex<double>,1>* in){

  ofstream ofs(filename);
  if (ofs.bad())
    {
      cerr << "Unable to write to file: " << filename << endl;
      exit(1);
    } 
  
  ofs << *in << endl;
  
}//write





Array<complex<double>,2> mult(Array<complex<double>,2> A, 
				     Array<complex<double>,2> B){
  
  if(A.rows() != B.rows()){
    cout << "Matrices do not match, multiplication impossible." << endl;
    
  }//if

  else {
    Array<complex<double>,2> temp(A.rows(), A.rows());
    
    
    firstIndex i;
    secondIndex j;
    thirdIndex k;
    
    temp = sum(A(i,k)*B(k,j), k);
    
    
    
    return temp;
  }

}//mult


Array<complex<double>,2> outerProduct(Array<complex<double>,2> A, 
				     Array<complex<double>,2> B){

  firstIndex i;
  secondIndex j;
  
  int dimA = A.rows();
  int dimB = B.rows();
  int dim = dimA*dimB;
  Array<complex<double>,2> C(dim, dim), D(dimB, dimB);

  
  for (int a=0; a < dimA; a ++){
    for (int b=0; b < dimA; b ++){
     
      C(Range(a*dimB, (a+1)*dimB - 1), Range(b*dimB, (b+1)*dimB -1)) = A(a,b)*B(i,j);
   
    }
  }
 
 
  return C;


}//finito

double innerProduct(Array<complex<double>,1> A){


  double c=0.;
  
  for(int i=0; i < A.rows(); i++){
   
    c += real(A(i))*real(A(i));
    c += imag(A(i))*imag(A(i));

  }



  return c;

}//innerProduct
void singleQubitGateNoMatrix(Array<complex<double>,2> unitary, int k) {

  
  *helpVector1 = 0.0;
  
  int a=(int)pow(2.0, (double)(k-1));
  int b=(int)pow(2.0, (double)(noOfQubits-k));  

  for(int i=0; i < a; i++) {
    for(int j=0; j< b; j++) {
            
      (*helpVector1)( i*2*b +j) +=
	unitary(0,0)* (*state)(i*2*b +j) 
	+ unitary(0,1)* (*state)(i*2*b +b +j);
      
      (*helpVector1)(i* 2*b+b +j) +=
	unitary(1,0)* (*state)(i*2*b +j) 
	+ unitary(1,1)* (*state)(i*2*b +b +j);
      
    }//for
  }//for
  

  operationsPerformed ++;
  *state=*helpVector1;

}//U


void conditionalTwoQubitGateNoMatrix(Array<complex<double>,2> unitary, int k, int l) {

  *helpVector1 = *state;
  
  //First part of the CU
  int a=(int)pow(2.0, (double)(k-1));
  int b=(int)pow(2.0, (double)(noOfQubits-k));  

  for(int i=0; i<a; i++) {
    for(int j=0; j<b; j++) {
      
      (*helpVector1)(i*2*b+b+j)=0.0;

    }//for
  }//for

  //cout << " h1 "<< helpVector1<< endl;



  if( k < l ) {  
    //Second part of the CU
    int c=(int)pow(2.0, (double)(l-1-k));
    int d=(int)pow(2.0, (double)(noOfQubits-l));
    
    *helpVector2=0.0;
    
    for(int i=0; i < a; i++) {
      for(int j=0; j < c; j++) {
	for(int w=0; w < d; w++) {
	  
	  (*helpVector2)( i*2*b +b +j*2*d  +w) +=
	    unitary(0,0)* (*state)( i*2*b +b +j*2*d +w)
	    +unitary(0,1)* (*state)( i*2*b +b +j*2*d +d +w);
	  
	  (*helpVector2)( i*2*b +b +j*2*d +d +w) +=
	    unitary(1,0)* (*state)( i*2*b +b +j*2*d +w)
	    +unitary(1,1)* (*state)( i*2*b +b +j*2*d +d +w);
	  
	}//for
      }//for
    }//for
    
  }//if


  else {
    int c=(int)pow(2.0, (double)(k-1-l));
    int d=(int)pow(2.0, (double)(noOfQubits-l));
    int e=(int)pow(2.0, (double)(l-1));

    *helpVector2=0;
    
    for(int i=0; i < e; i++) {
      for(int j=0; j < c; j++) {
	for(int w=0; w < b; w++) {
	  
	  (*helpVector2)( i*2*d  +j*2*b +b  +w) +=
	    unitary(0,0)* (*state)( i*2*d +j*2*b +b +w)
	    +unitary(0,1)* (*state)( i*2*d +d +j*2*b +b +w);
	  
	  (*helpVector2)( i*2*d +d +j*2*b +b +w) +=
	    unitary(1,0)* (*state)( i*2*d  +j*2*b +b +w)
	    +unitary(1,1)* (*state)( i*2*d +d +j*2*b +b +w);
	  
	}//for
      }//for
    }//for
    
  }//else




  //cout << " h2 " << helpVector2 << endl;

  operationsPerformed += 6;
  *state= *helpVector1 + *helpVector2;

}// conditionalTwoQubitGate()




/*This is for use in generating U for the simulation qubits*/
Array<complex<double>,2>  singleQubitGate(Array<complex<double>,2> UU, 
		  int qubitOperatedOn) {
  
		  
  Array<complex<double>,2> temp(simDim, simDim);
  
   
  int dima = (int) pow(2.0,(double)(qubitOperatedOn - 1)); //eksponent
  int dimb = (int) pow(2.0,(double)(simQubits - qubitOperatedOn));

  Range a(firstDim, dima -1), b(simDim - dimb , simDim -1), 
    tmp(simDim - 2*dimb , simDim -1), dim(firstDim, simDim-1);

  
  for (int i=0; i < simDim; i++){  
    temp(i,i)=1;                      //unity matrix for whole system
  }//for

 

  if( qubitOperatedOn == secondDim) { 
    temp = outerProduct( UU, temp(b, b));
    
  }//if
  //Only need one outerProduct because first qubit is operated on

  else{
    if( qubitOperatedOn == simQubits){
    
      temp = outerProduct( temp(a, a), UU);

    }//Here it's the last qubit we operate on
    


    else{
      temp(tmp, tmp) = outerProduct( UU, temp(b, b));
      temp = outerProduct(temp(a,a), temp(tmp,tmp));
      //This is OP tensor identity matrix for the rest of the qubits
    }//else
  }//else


  return temp;

}//singleQubitGate for Umatrix generation




/*This is the method for calculating the matrix for a general CU operation 
  on the state space*/
Array<complex<double>,2> gCU(Array<complex<double>,2> UU, int c, int u){

  /*In the case of multiqubit U operators, u is defined as the first of
    the qubits the Gate operate on, only one qubit Gates or next neighbour
    Gates are available. That is why we have two sets of dimensions for the
    help matrices*/

  int a, b, d, e, uDim;
  b = 1;
  uDim = UU.rows();
 


  a = (int) pow(2.0,(double)(u - c - 1)); //dimensions of help matrices
  d = (int) pow(2.0,(double)(simQubits - c + 1));
  e = (int) pow(2.0,(double)(c - 1));
  
  Array<complex<double>,2> temp(simDim, simDim), 
    temp2(simDim, simDim), temp3(simDim, simDim),
    upper(thirdDim, thirdDim),
    lower(thirdDim, thirdDim);// help matrices
  
  
  temp3 = 0;
  temp2 = 0;
  temp = 0;
  for (int i=0; i < simDim; i++){  
    temp2(i,i)=1;                      //unity matrix for whole system
  }//for
  
  upper = 1,0,
    0,0;//This is part of the Control operation
  
  lower = 0,0,
    0,1;
  
  temp3 = singleQubitGate(upper, c);
    
  
  //Ranges for help matrices
  
  Range two(firstDim, 2*a - 1), I2(firstDim, a-1), three(firstDim, d - 1),
    I3(firstDim, e - 1);
 
  
  b = (int) pow(2.0,(double)(simQubits - u + 1));
  
  Range one(simDim - b, simDim -1), I1(simDim - b/uDim, simDim -1);
  
  temp(one,one) = outerProduct(UU, temp2(I1, I1));
  temp(two,two) = outerProduct(lower, temp2(I2,I2));
  temp(three, three) = outerProduct( temp(two, two), temp(one, one));
  temp = outerProduct(temp2(I3, I3), temp(three, three));
  temp3 += temp;
    
    
  return temp3;
}//CU 



Array<complex<double>,2> gUC(Array<complex<double>,2> UU, int u, 
			    int c){
 /*In the case of multiqubit U operators, u is defined as the first of
    the qubits the Gate operate on, only one qubit Gate or next neighbour
    Gates are available. That is why we have two sets of dimensions for the
    help matrices*/

  int a, b, d, e, uDim;
  b = 1;
  uDim = UU.rows();

 
 
  //dimensions of help matrices, independent of size of U
  b = (int) pow(2.0,(double)(simQubits - c));
  d = (int) pow(2.0,(double)(simQubits - u + 1));
  e = (int) pow(2.0,(double)(u - 1));
  
  Array<complex<double>,2> temp(simDim, simDim), 
    temp2(simDim, simDim), temp3(simDim, simDim),
    upper(thirdDim, thirdDim),
    lower(thirdDim, thirdDim);// help matrices
  
  
  temp3 = 0;
  temp2 = 0;
  temp = 0;
  for (int i=0; i < simDim; i++){  
    temp2(i,i)=1;                      //unity matrix for whole system
  }//for
  
  upper = 1,0,
    0,0;//This is part of the Control operation
  
  lower = 0,0,
    0,1;
  
  temp3 = singleQubitGate(upper, c);
    
  
  //Ranges for help matrices
  
  
  Range one(simDim - 2*b, simDim -1), I1(simDim - b, simDim -1), 
    three(firstDim, d - 1), I3(firstDim, e - 1);
     
  
  if(uDim==2){
    a = (int) pow(2.0,(double)(c - u - 1));
    Range two(firstDim, 2*a - 1), I2(firstDim, a-1);
    
    temp(one,one) = outerProduct(lower, temp2(I1, I1));
    temp(two,two) = outerProduct(UU, temp2(I2,I2));
    temp(three, three) = outerProduct( temp(two, two), temp(one, one));
    temp = outerProduct(temp2(I3, I3), temp(three, three));
    temp3 += temp;
    
    
  }//if
  else{
    a = (int) pow(2.0,(double)(c - u - 2));
    Range two(firstDim, 4*a - 1), I2(firstDim, a-1);
    
    
    temp(one,one) = outerProduct(UU, temp2(I1, I1));
    temp(two,two) = outerProduct(lower, temp2(I2,I2));
    temp(three, three) = outerProduct( temp(two, two), temp(one, one));
    temp = outerProduct(temp2(I3, I3), temp(three, three));
    
    temp3 += temp;
    
  }//else
  
  return temp3;

}//gUC


void ZZ ( int i, int j, double a) {

  Array<complex<double>,2> D(2,2);
  D = exp(complex<double> (0, -2*a)), 0,
    0, exp(complex<double> (0, -2*a));
  Array<complex<double>,2> E(2,2);
  E = exp(complex<double> (0, a)), 0,
    0, exp(complex<double> (0, a));

  if ( i < j) {
      
    *Umatrix = mult(singleQubitGate(E,i),*Umatrix);
    
    *Umatrix = mult(singleQubitGate(x,i),*Umatrix);
    *Umatrix = mult(gCU(x, i, j), *Umatrix);
    *Umatrix = mult(singleQubitGate(x,i),*Umatrix);
    
    *Umatrix = mult(singleQubitGate(x,j),*Umatrix);
    *Umatrix = mult(gUC(D, i, j), *Umatrix);
    *Umatrix = mult(singleQubitGate(x,j),*Umatrix);

    *Umatrix = mult(singleQubitGate(x,i),*Umatrix);
    *Umatrix = mult(gCU(x, i, j), *Umatrix);
    *Umatrix = mult(singleQubitGate(x,i),*Umatrix);

  }

   if ( j < i ) {

     *Umatrix = mult(singleQubitGate(E,i),*Umatrix);

     *Umatrix = mult(singleQubitGate(x,i),*Umatrix);
     *Umatrix = mult(gUC(x, j, i), *Umatrix);
     *Umatrix = mult(singleQubitGate(x,i),*Umatrix);

     *Umatrix = mult(singleQubitGate(x,j),*Umatrix);
     *Umatrix = mult(gCU(D, j, i), *Umatrix);
     *Umatrix = mult(singleQubitGate(x,j),*Umatrix);

     *Umatrix = mult(singleQubitGate(x,i),*Umatrix);
     *Umatrix = mult(gUC(x, j, i), *Umatrix);
     *Umatrix = mult(singleQubitGate(x,i),*Umatrix);
     
   }
   
}// end ZZ 

/*Her kommer en-partikkel metoden, ta en vilkårlig Eij element og gjøre alle 
  operasjonene */








void UPairing( int phaseMultiplier) {
  /*This method shall generate Umatrix for ONE run, then someone else 
    will take the matrix generated here and multiply it by itself
    a gazillion times*/

  /*The Hamiltonian parameters are d and g. The total time for the whole
    simulation is deltaT, divide by the noOfTimesteps and we have 
    what we are going to use here*/
  
  
  *Umatrix=0;
  for(int i=0; i< simDim; i++)
    (*Umatrix)(i,i)=1.0;
  
  double t= deltaT/(double)noOfTimesteps*(double)phaseMultiplier*g;
  
  //first we initialize all help matrices needed
  Array<complex<double>,2> AB(2,2);
  AB=I*exp(complex<double> (0,(Emax/g+simQubits/2.0)*t));
  
  
  //Til Z, må forandres når nivåene ikke er degenererte
  Array<complex<double>,2> Rz(2,2);//Because it is the rotation operator
  Rz = exp(complex<double> (0,-1)*t/2.0),0,
    0,exp(complex<double> (0,1)*t/2.0);
  
  //Til XX og YY
  Array<complex<double>,2> A(2,2);
  A = cos(-t/4.0), complex<double> (0, -sin(-t/4.0)),
    complex<double> (0, -sin(-t/4.0)), cos(-t/4.0);
  
  
  //Til YY aleine
  Array<complex<double>,2> C(2,2);
  C = cos(-t/4.0), complex<double> (0, sin(-t/4.0)),
    complex<double> (0,  sin(-t/4.0)), cos(-t/4.0);
  
  
  
  
  
  *Umatrix = mult(singleQubitGate(AB,1),*Umatrix);
  
  for(int i=1; i <= simQubits; i ++){
    //Z on each qubit
    *Umatrix = mult(singleQubitGate(Rz,i),*Umatrix);
    
    for(int j=1; j<=simQubits; j++){
      //cout << "hey hey!, i = " << i << ", and j = " << j << endl;
      if( i < j) {
	
	//YY
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gCU(x, i, j), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	*Umatrix = mult(gUC(A, i, j), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	
	*Umatrix = mult(gUC(C, i, j), *Umatrix);
	
	//XX
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	*Umatrix = mult(gUC(A, i, j), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	
	*Umatrix = mult(gUC(A, i, j), *Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gCU(x, i, j), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
      }// if 
      
      if ( j < i ) {
	
	//YY
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gUC(x, j, i), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	*Umatrix = mult(gCU(A, j, i), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	
	*Umatrix = mult(gCU(C, j, i), *Umatrix);
	
	//XX
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	*Umatrix = mult(gCU(A, j, i), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	
	*Umatrix = mult(gCU(A, j, i), *Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gUC(x, j, i), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	
      }//if
      
      
      
    }// for j
  }//for i
  
  *helpMatrix = *Umatrix;
  
  for(int a=0; a < noOfTimesteps -1; a++){
    *helpMatrix = mult(*Umatrix, *helpMatrix);
  }// for a
  
  *Umatrix = *helpMatrix;
  
  
}// UPairing




void UpairingNonDegenerate(int phaseMultiplier) {

  /*This method shall generate Umatrix for ONE run, then someone else 
    will take the matrix generated here and multiply it by itself
    a gazillion times*/

  /*The Hamiltonian parameters are d and g. The total time for the whole
    simulation is deltaT, divide by the noOfTimesteps and we have 
    what we are going to use here*/


  *Umatrix=0;
  for(int i=0; i< simDim; i++)
    (*Umatrix)(i,i)=1.0;
  
  double t= deltaT/(double)noOfTimesteps*(double)phaseMultiplier;
  
  //first we initialize all help matrices needed
  Array<complex<double>,2> AB(2,2);
  AB=I*exp(complex<double> (0,(Emax+simQubits*g/2.0 
			       -simQubits*(simQubits-1)/2.0*dd)*t));
  
  

  //Til Z, må forandres når nivåene ikke er degenererte
  Array<complex<double>,2> Rz(2,2);//Because it is the rotation operator


  //Til XX og YY
  Array<complex<double>,2> A(2,2);
  A = cos(-t*g/4.0), complex<double> (0, -sin(-t*g/4.0)),
    complex<double> (0, -sin(-t*g/4.0)), cos(-t*g/4.0);
  
  
  //Til YY aleine
  Array<complex<double>,2> C(2,2);
  C = cos(-t*g/4.0), complex<double> (0, sin(-t*g/4.0)),
    complex<double> (0,  sin(-t*g/4.0)), cos(-t*g/4.0);



  *Umatrix = mult(singleQubitGate(AB,1),*Umatrix);
  
  for(int i=1; i <= simQubits; i ++){
    //Z on each qubit
   //  //here we have different Rz for each i
    Rz = exp(complex<double> (0,-1)*t*(0.5*g-(i-1)*dd)),0,
      0,exp(complex<double> (0,1)*t*(0.5*g-(i-1)*dd));


    


    *Umatrix = mult(singleQubitGate(Rz,i),*Umatrix);
    
    for(int j=1; j<=simQubits; j++){
      //cout << "hey hey!, i = " << i << ", and j = " << j << endl;
      if( i < j) {
	
	//YY
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gCU(x, i, j), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	*Umatrix = mult(gUC(A, i, j), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	
	*Umatrix = mult(gUC(C, i, j), *Umatrix);
	
	//XX
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	*Umatrix = mult(gUC(A, i, j), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	
	*Umatrix = mult(gUC(A, i, j), *Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gCU(x, i, j), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
      }// if 
      
      if ( j < i ) {

	//YY
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gUC(x, j, i), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	*Umatrix = mult(gCU(A, j, i), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	
	*Umatrix = mult(gCU(C, j, i), *Umatrix);
	
	//XX
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	*Umatrix = mult(gCU(A, j, i), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	  
	*Umatrix = mult(gCU(A, j, i), *Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gUC(x, j, i), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	
      }//if
      
      
      
    }// for j
  }//for i
  
  *helpMatrix = *Umatrix;

  for(int a=0; a < noOfTimesteps -1; a++){
    *helpMatrix = mult(*Umatrix, *helpMatrix);
  }// for a

  *Umatrix = *helpMatrix;


}// UpairingNonDegenerate



void UHeisenberg(int phaseMultiplier) {
  /*This method shall generate Umatrix for ONE run, then someone else 
    will take the matrix generated here and multiply it by itself
    a gazillion times*/

  /*The Hamiltonian parameters are d and g. The total time for the whole
    simulation is deltaT, divide by the noOfTimesteps and we have 
    what we are going to use here*/

  *Umatrix=0;
  for(int i=0; i< simDim; i++)
    (*Umatrix)(i,i)=1.0;
  
  double t= deltaT/(double)noOfTimesteps*(double)phaseMultiplier;
  
  //first we initialize all help matrices needed
  Array<complex<double>,2> AB(2,2);
  AB=I*exp(complex<double> (0,Emax*t)); //remove
  //AB=I*exp(complex<double> (0,(Emax+simQubits/2.0)*t));

 
  //remove
  t= t*2.0;
  t= t*2.0;
  
  //Til Z, må forandres når nivåene ikke er degenererte
  Array<complex<double>,2> Rz(2,2);//Because it is the rotation operator
  Rz = exp(complex<double> (0,-1)*t/4.0),0,
    0,exp(complex<double> (0,1)*t/4.0);
  //remove
  t=t/2.0;


  //remove
  t=-t;

  //Til XX og YY
  Array<complex<double>,2> A(2,2);
  A = cos(-t/2.0), complex<double> (0, -sin(-t/2.0)),
    complex<double> (0, -sin(-t/2.0)), cos(-t/2.0);
  
  
  //Til YY aleine
  Array<complex<double>,2> C(2,2);
  C = cos(-t/2.0), complex<double> (0, sin(-t/2.0)),
    complex<double> (0,  sin(-t/2.0)), cos(-t/2.0);

  //Til ZZ 
  Array<complex<double>,2> D(2,2), E(2,2);
  D = exp(complex<double> (0,2*(-t/2.0))), 0,
    0, exp(complex<double> (0,2*(-t/2.0)));

  E =  exp(complex<double> (0,-(-t/2.0)/2.0)), 0,
    0, exp(complex<double> (0,-(-t/2.0)/2.0));
  

  //remove
  t=-t;

  //remove
  t=t/2.0;

  //cout << "Rz = " << Rz << endl;
  //cout << "AB = " << AB << endl;

  //cout << "Umatrix = " << *Umatrix << endl;

  for(int a=0; a < noOfTimesteps; a++){

    *Umatrix = mult(singleQubitGate(AB,1),*Umatrix);
    //cout << "AB som fullstendig matrise = " << singleQubitGate(AB,1) << endl;
    

    for(int i=1; i <= simQubits; i ++){
      //Z on each qubit
      //cout << " i " << i << endl;
      *Umatrix = mult(singleQubitGate(Rz,i),*Umatrix);
      //cout << "Rz som fullstendig matrise = " << singleQubitGate(Rz,i) << endl;
      


      // xx og yy på i<j
      
 
      int j = i+1;
 
    //   //ZZ, remove
      if (i < simQubits ) {

	//ZZ
	*Umatrix = mult(singleQubitGate(E,i),*Umatrix);
	*Umatrix = mult(singleQubitGate(E,j),*Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gCU(x, i, j), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	*Umatrix = mult(gUC(D, i, j), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gCU(x, i, j), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);


	//YY
      
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gCU(x, i, j), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	*Umatrix = mult(gUC(A, i, j), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	
	*Umatrix = mult(gUC(C, i, j), *Umatrix);
	
	//XX
	
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	*Umatrix = mult(gUC(A, i, j), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,j),*Umatrix);
	
	*Umatrix = mult(gUC(A, i, j), *Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gCU(x, i, j), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	
      }




      else {
	//ZZ
	*Umatrix = mult(singleQubitGate(E,1),*Umatrix);
	*Umatrix = mult(singleQubitGate(E,i),*Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gUC(x, 1, i), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,1),*Umatrix);
	*Umatrix = mult(gCU(D, 1, i), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,1),*Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gUC(x, 1, i), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);


	//YY
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gUC(x, 1, i), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,1),*Umatrix);
	*Umatrix = mult(gCU(A, 1, i), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,1),*Umatrix);
	
	*Umatrix = mult(gCU(C, 1, i), *Umatrix);

	//XX
	*Umatrix = mult(singleQubitGate(x,1),*Umatrix);
	*Umatrix = mult(gCU(A, 1, i), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,1),*Umatrix);
	
	*Umatrix = mult(gCU(A, 1, i), *Umatrix);
	
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	*Umatrix = mult(gUC(x, 1, i), *Umatrix);
	*Umatrix = mult(singleQubitGate(x,i),*Umatrix);
	
      }

      
    }// simQubits for
    
  }// for noOfTimesteps
  
  //cout << "Umatrix  after Rz= " << *Umatrix << endl;
  

}//end U






void CU(int c) {

  /*This method takes the whole U matrix for the simulation qubits and
    performs the contolled U operatin on the total state without 
    a matrix*/

  /* c is the controlling work qubit*/
  *helpVector1 = *state;
  *helpVector2 = 0;


  int rngI = (int) pow(2.0,(double)(c-1));
  int rngJ = (int) pow(2.0,(double)(workQubits - c));
  int rngJDobl = (int) pow(2.0,(double)(workQubits - c +1));
  int halfRng = (int) pow(2.0,(double)(noOfQubits - c));


  for( int i = 0; i<rngI; i++)
    for( int j = 0; j<halfRng; j++)
      (*helpVector1)((i*2+1)*halfRng + j) = 0.0;



  for( int i = 0; i< rngI; i++) {
    for( int j = 0; j < rngJ; j++) {
      
      //now the U on x matrix vector multiplications
      for( int k = 0; k < simDim; k++) {
	

	for( int l = 0; l < simDim; l++) {
	  (*helpVector2)((2*i+1)*halfRng + j*simDim + k) +=
	    (*Umatrix)(k,l)*(*state)((2*i+1)*halfRng + j*simDim + l);
		      
	  
	}//for
      }//for
    }//for
  }//for
  
  *state = *helpVector1 + *helpVector2;


}//end CU Denna metoden tar matriseløs multiplikasjon på tilstandsvektoren







void inverseFourierTransform(){

  Array<complex<double>,2> Rk(2,2);

  for(int i=1; i < firstSimQubit; i++){

    /*Rks for the below work qubits*/
    for(int j=1; j < i; j++){
      Rk = 1, 0,
	0, exp(complex<double> (0, -2*pi/pow(2.0,(double)(i-j + 1))));
      //cout << "her tenker jeg feilen ligger" << endl;
      //qState->operate(CU(qGate, Rk, j, i));  
      conditionalTwoQubitGateNoMatrix(Rk, j, i);

    }//for
    
    //qState->operate(singleGate(qGate, H, i));
    singleQubitGateNoMatrix(H,i);

  }//for 
  
  
}//IFT







void findingProbability(){

  int slot = 0;
  
  /*The probabilities*/
  *probability = 0;
 
  
  for(int k = 0; k < workDim; k++){
    Range nonZero(k*simDim, (k+1)*simDim - 1);
    /*Finding the probabilities for each eigenvalue*/
    (*probability)(k) = innerProduct((*state)(nonZero));

  }//for



  /*Here is a test to see if the probabilitis sum up to one*/
  double sum = 0;
  for (int i = 0; i < workDim; i++)
    sum += (*probability)(i);

  cout << "Er sannsynet 1? " << sum << endl;

}//probability




void measurement(){

  double m=0;
  
  /*This is where each measurement goes*/
  *measured = 0;
 
  /*Times each eigenvalue is measured*/
 
  *timesMeasured = 0;

 
  /*Generate a random number to pick an eigenvalue*/
  for(int i=0; i < noOfMeasurements; i++){
    m = (double) rand()/RAND_MAX;
    
    int j =0;
    double totalProb = (*probability)(firstDim);
    while( (m > totalProb) && (j < workDim) ){
      
      j++;
      totalProb += (*probability)(j);

    }//while
    (*measured)(i) = j;
    (*timesMeasured)(j) += 1;

  }//for

  writeToFile("everyMeasurement.data", measured);
  writeToFile("frequency.data", timesMeasured);

  
  *probs = 0;
  
  for(int i=0; i < workDim; i++)
    (*probs)(i) = (double)(*timesMeasured)(i)/(double) noOfMeasurements;
  
  writeToFile("measuredProb.data", probs);

}//measurement


void writeOutEnergies(){
  /*This method writes out the energies corresponding to 
    the binary values of the work qubits*/

  double t =deltaT;
  double fie;
  int tamp;
  Array<int,1> wb(workQubits);
  Array<double,1> energies(workDim), phi(workDim), energyScale(workDim), 
    probOfEnergy(workDim);
 
 
  for(int i=0; i<workDim; i++){
    wb=0;
    fie=0;
    tamp=i;
    for(int j=0; j<workQubits; j++){
      wb(j)=tamp%2;
      tamp = tamp/2;
      fie+=wb(j)*pow(2.0,(double)(workQubits-j -1));
    }
    //cout<<"wb    " << wb<<endl;
    phi(i)=(double)fie/(double)workDim;
    energies(i) = phi(i)*(-2.0*pi/t); 
    energyScale((int)fie)=energies(i);
   
    probOfEnergy((int)fie)=(*probability)(i);
    
  }  
  
  for(int k=0; k<workDim; k++)
    energyScale(k)=energyScale(k)+Emax;
  

  writeToFile("scale.data", &energyScale);
  writeToFile("EnergyProb.data", &probOfEnergy);
  writeToFile("energies.data", &energies);
  writeToFile("phis.data", &phi);

}//writeOutEnergies



void initializeState(){
  /*This method uses the pairing N/2 state*/

  //initializing the state vector
  
  Array<complex<double>,1> temp(simDim);

  *state = 0;
  *helpVector1 = 0;
  *helpVector2 = 0;

  double m, k;
  for(int i=0; i < simDim; i++){
    m = (double) rand()/RAND_MAX;
    k = (double) rand()/RAND_MAX;
    temp(i)=complex<double>(m,k);
  }

  int y=0;
  int countZeros;

  

  for(int i=0; i< simDim; i++) { //for every entry of vec
    
    countZeros =0;
    //y is the binary zeros and ones of the index of vec
    for(int j=0; j<simQubits; j++) {
      y=(i % (int) pow(2.0,(double)j+1)) / (int) pow(2.0,(double)j);
      //cout << "\ni = " << i << " j = " << j << " y = " << y <<endl;
      if (y == 0)
	countZeros +=1;
      
    }//for
    
    //cout <<"\nAss, i= " << i << endl;
 

    // remove uncomment
    // if(countZeros != (simQubits/2))
//       temp(i)=0.0;
    
  }//for


  //normalization
  m=innerProduct(temp);
  for(int i=0; i < temp.rows(); i++){
    temp(i)/=sqrt(m);
  }
  
  
  

  //alternative initialization remove
  // *state = 0;
//   ifstream ifsRe("vectorReal.data");
//   ifstream ifsIm("vectorImag.data");
//   if (ifsRe.bad() || ifsIm.bad())
//     {
//       cerr << "Unable to open file: " << " vector*.data"  << endl;
//       exit(1);
//     }

//   double a, b;
//   temp=0;
//   for ( int i=0; i<simDim; i++) {  
//     ifsRe >> a;
//     ifsIm >> b;
//     temp(i)= complex<double>(a,b);
//   }
  
  //end remove
  

  Range I(0, simDim - 1);
  (*state)(I) = temp;
  cout << "The randomly generated simulation vector = " << temp << endl;


  //done initializing the state vector

} // done initializeState()


//her kommer alle kommandoene
void hoved() {

 


  



  //initializing the U matrix
 
  *Umatrix = 0;
  *helpMatrix = 0;
  for(int i=0; i<simDim; i++)
    (*Umatrix)(i,i)=1.0;


  initializeState();
    
  // //enPartikkel(1,1,0.63486,1);
//   enPartikkel(1,2,1,0.1);
//   cout << "E11=1, dt=1, U= " << *Umatrix<< endl;



  for(int i = 1; i <= workQubits; i++) 
    singleQubitGateNoMatrix(H,i);


  // cout << "\nState etter Hene = " << *state << endl;

  //phaseEst:
  //løkke over arbeidsqubitene
  //  lag U, send til parallellisering, få igjen U^n
  //  CU på vektoren


  for (int i=1; i<=workQubits; i++){
    // 1= heisenberg, 2=pairing, 3=pairinNonDeg, 4=Utest, 5=Hubbard
    //UHeisenberg((int) pow(2.0,(double)(i-1)));
    //UPairing((int) pow(2.0,(double)(i-1)));
    //UpairingNonDegenerate((int) pow(2.0,(double)(i-1)));    
    Utest((int) pow(2.0,(double)(i-1)));
    //UHubbard((int) pow(2.0,(double)(i-1)));
    
    //cout << "This is U, for exponent = " << i << "   " << *Umatrix << endl;
    CU(workQubits -i +1);
    //cout << "\nState = " << *state << endl;
  }//for


  //(*Umatrix)=singleQubitGate(x,1);
  //cout << "This is Umatrix " << *Umatrix << endl;
  //CU(2);


  inverseFourierTransform();
  //cout << "\nState = " << *state << endl;

  
  findingProbability();
 

  


}


int main( int argc, char* argv[] ) {
  

  vFile = "vector.data";
  mFile = "matrix.data";
  simFile = "simulation.data";
  pi = 3.141592653589793238;

  x = 0,1,
    1,0;
  y = 0, complex<double> (0, -1),
    complex<double> (0, 1), 0;
  z = 1, 0,
    0, -1;
  H = 1, 1,
    1, -1;
  H /= sqrt(2.0);
  I = 1,0,
    0,1;
  

  //seeding the random number generator
  srand(time(0)); 
  //srand(1);
  int t = time(0);

  
  //fetching the data for this simulation
  readData();
  dimension = (int) pow(2.0,(double)(noOfQubits));
  simQubits = noOfQubits - firstSimQubit + 1;
  workQubits = firstSimQubit - 1;
  simDim = (int) pow(2.0,(double)(simQubits)); 
  workDim = (int) pow(2.0,(double)(workQubits)); 


  //creating all arrays
  state = new Array<complex<double>,1>(dimension);
  helpVector1 = new Array<complex<double>,1>(dimension);
  helpVector2 = new Array<complex<double>,1>(dimension);
  Umatrix = new Array<complex<double>,2>(simDim, simDim);
  helpMatrix = new Array<complex<double>,2>(simDim, simDim);
  measured = new Array<double,1>(noOfMeasurements);
  timesMeasured = new Array<int,1>(workDim);
  probs = new Array<double,1>(workDim);


  probability = new Array<double,1>(workDim);
  totProb = new Array<double,1>(workDim);
  *probability = 0;
  *totProb = 0;
  for (int counter =1; counter <= noOfRuns; counter++) {

    hoved();

    for (int j = 0; j<workDim; j++)
      (*totProb)(j) += (*probability)(j);
  }

  for (int j = 0; j<workDim; j++)
    (*totProb)(j) /= (double) noOfRuns;

  *probability = *totProb;


  measurement();
  writeOutEnergies();



//tester det nye programmet
  

  t = time(0) - t;
  cout << "\nProgram running time:    " << t << " seconds.\n"
       << " Is " << t/3600 << " hours, "
       << (t%3600)/60 << " minutes and " << t%60 << " seconds." << endl;

  return 0;

}// end main






  


