#include <math.h>
#include <iostream> 
#include <stdlib.h>
#include <complex>
#include <fstream>
#include <blitz/array.h>
#include <time.h>
#include "ny.h"



using namespace blitz;

void enPartikkel( int i, int j, double E, double dt)
{
  /*kommentar: i > j funker ikke, gir uendelig sannsyn*/
 
  double tall = 1.0/sqrt(2.0);
  Array<complex<double>,2> A(2,2);

  if (i==j) 
    {
     
      A = exp(complex<double> (0, -E*dt)), 0,
	0, 1;
      
      *Umatrix = mult(singleQubitGate(A,i),*Umatrix);
    }//end i == j
  
  else 
    { 
      dt=-0.5*dt;

      //Først xzzz---zzzx
      //if j-i = odde
      if ( (j-i) %2 ==1 )
	{
	  A = complex<double> (tall, -tall) ,0, 
	    0, complex<double> ( tall, tall);
	  *Umatrix = mult(singleQubitGate(A,i),*Umatrix);

	  A= tall, tall,
	    -tall, tall;
	  *Umatrix = mult(singleQubitGate(A,j),*Umatrix);

	  if ( ((j-i)/2) % 2 == 1)
	    ZZ(i,j, -pi/4.0);
	  else
	    ZZ(i,j, pi/4.0);
	    
	}//end j-i = odde
      else //(j-i) = partall
	{

	  A = tall, tall,
	    -tall, tall;
	  *Umatrix = mult(singleQubitGate(A,j),*Umatrix);

	  if ( ((j-i-1)/2) % 2 == 1)
	    ZZ(i,j, -pi/4.0);
	  else
	    ZZ(i,j, pi/4.0);

	}//end (j-i) = partall


      for( int k = j-i; k >=2; k--)
	ZZ(i, k+i-1, pi/4);

      //U1
      A = tall, tall,
	-tall, tall;
      *Umatrix = mult(singleQubitGate(A,i),*Umatrix);

      //U
      A = complex<double> (cos(E*dt), -i*sin(E*dt)), 0, 
	0, complex<double> (cos(E*dt), i*sin(E*dt));
      *Umatrix = mult(singleQubitGate(A,i),*Umatrix);

      //så begynner vi baklengs....

      //U1*
      A = tall, -tall,
	tall, tall;
      *Umatrix = mult(singleQubitGate(A,i),*Umatrix);

       for( int k = 2; k <=j-i; k++)
	 ZZ(i, k+i-1, -pi/4);


       //if j-i = odde
       if ( (j-i) %2 ==1 )
	 {	   
	   if ( ((j-i)/2) % 2 == 1)
	     ZZ(i,j, pi/4.0);
	   else
	     ZZ(i,j, -pi/4.0);
	    
	   A= tall, -tall,
	     tall, tall;
	   *Umatrix = mult(singleQubitGate(A,j),*Umatrix);

	   A = complex<double> (tall, tall) ,0, 
	     0, complex<double> ( tall, -tall);
	   *Umatrix = mult(singleQubitGate(A,i),*Umatrix);

	 }//end j-i = odde
       else //(j-i) = partall
	 {

	 

	   if ( ((j-i-1)/2) % 2 == 1)
	     ZZ(i,j, pi/4.0);
	   else
	     ZZ(i,j, -pi/4.0);

	   A = tall, -tall,
	     tall, tall;
	   *Umatrix = mult(singleQubitGate(A,j),*Umatrix);

	 }//end (j-i) = partall






       //så yzz----zzy

       //if j-i = odde
       if ( (j-i) %2 ==1 )
	{
	  A = complex<double> (tall, -tall) ,0, 
	    0, complex<double> ( tall, tall);
	  *Umatrix = mult(singleQubitGate(A,i),*Umatrix);

	  A= tall, complex<double>(0,tall),
	    complex<double> (0, tall), tall;
	  *Umatrix = mult(singleQubitGate(A,j),*Umatrix);

	  if ( ((j-i-1)/2) % 2 == 1)
	    ZZ(i,j, -pi/4.0);
	  else
	    ZZ(i,j, pi/4.0);
	    
	}//end j-i = odde
      else //(j-i) = partall
	{
	  A= tall, complex<double>(0,tall),
	    complex<double> (0, tall), tall;
	  *Umatrix = mult(singleQubitGate(A,j),*Umatrix);

	  if ( ((j-i-1)/2) % 2 == 1)
	    ZZ(i,j, -pi/4.0);
	  else
	    ZZ(i,j, pi/4.0);

	}//end (j-i) = partall


      for( int k = j-i; k >=2; k--)
	ZZ(i, k+i-1, pi/4);

      //U1
      A= tall, complex<double>(0,tall),
	complex<double> (0, tall), tall;
      *Umatrix = mult(singleQubitGate(A,i),*Umatrix);

      //U
      A = complex<double> (cos(E*dt), -i*sin(E*dt)), 0, 
	0, complex<double> (cos(E*dt), i*sin(E*dt));
      *Umatrix = mult(singleQubitGate(A,i),*Umatrix);

      //så begynner vi baklengs....

      //U1*
      A= tall, complex<double>(0,-tall),
	    complex<double> (0, -tall), tall;
      *Umatrix = mult(singleQubitGate(A,i),*Umatrix);

       for( int k = 2; k <=j-i; k++)
	 ZZ(i, k+i-1, -pi/4);

       //if j-i = odde
       if ( (j-i) %2 ==1 )
	 {	   
	   if ( ((j-i-1)/2) % 2 == 1)
	     ZZ(i,j, pi/4.0);
	   else
	     ZZ(i,j, -pi/4.0);
	    
	   A= tall, complex<double>(0,-tall),
	    complex<double> (0, -tall), tall;
	   *Umatrix = mult(singleQubitGate(A,j),*Umatrix);

	   A = complex<double> (tall, tall) ,0, 
	     0, complex<double> ( tall, -tall);
	   *Umatrix = mult(singleQubitGate(A,i),*Umatrix);

	 }//end j-i = odde
       else //(j-i) = partall
	 {

	   if ( ((j-i-1)/2) % 2 == 1)
	     ZZ(i,j, pi/4.0);
	   else
	     ZZ(i,j, -pi/4.0);

	   A= tall, complex<double>(0,-tall),
	    complex<double> (0, -tall), tall;
	   *Umatrix = mult(singleQubitGate(A,j),*Umatrix);

	 }//end (j-i) = partall





    }//end i ~= j


}//end enPartikkel




/*Her kommer toPartikkel metoden, den skal ta et vilkårlig topartikkel ledd i H 
  og utføre alle operasjonene som trengs for å gi U av det leddet på 
  toqubitoperasjons form.*/

void toPartikkel( int i, int j, int k, int l, double V, double dt) 
{

  //først sjekke om vi er på diagonalen
  if ( i == l && j == k) 
    {




    }//end Vijij



  if (i == k && j == l)
    {




    }//end Vijji



}//end toPartikkel
 


void UHubbard(int phaseMultiplier) {

 
  if (simQubits%2 ==1)
    cout << "Moron! We need two qubits per level in the spin 1/2 Hubbard model" 
	 << endl;

  double t= deltaT/(double)noOfTimesteps*(double)phaseMultiplier;

  *Umatrix=0;
  for(int i=0; i< simDim; i++)
    (*Umatrix)(i,i)=1.0;

  //first we initialize all help matrices needed
  Array<complex<double>,2> AB(2,2);
  AB= exp(complex<double> (0, Emax*t)), 0,
    0, exp(complex<double> (0, Emax*t));

  *Umatrix = mult(singleQubitGate(AB,1),*Umatrix);

  

  //selve Hamiltonfunksjonen
  double epsilon = 1;
  double T = -1;
  double U=1;



  for( int i=1; i <=simQubits; i++)
    enPartikkel(i,i, epsilon, t);

  for( int i=1; i <= simQubits -2; i++)
    enPartikkel(i, i+2, T, t);
  
  
  

  //U^n
  *helpMatrix = *Umatrix;
  for(int a=0; a < noOfTimesteps -1; a++)
    *helpMatrix = mult(*Umatrix, *helpMatrix);
  
  *Umatrix = *helpMatrix;


}




void Utest(int phaseMultiplier) {
  
  double t= deltaT/(double)noOfTimesteps*(double)phaseMultiplier;

  *Umatrix=0;
  for(int i=0; i< simDim; i++)
    (*Umatrix)(i,i)=1.0;

  //first we initialize all help matrices needed
  Array<complex<double>,2> AB(2,2);
  AB= exp(complex<double> (0, Emax*t)), 0,
    0, exp(complex<double> (0, Emax*t));

  *Umatrix = mult(singleQubitGate(AB,1),*Umatrix);

  //enPartikkel(1,1, 1, t);
  //enPartikkel(1,2, 0.343, t);
  //enPartikkel(2,2, 1, t);
  enPartikkel(1,3, 1, t);
  //enPartikkel(2,4, 1, t);

  *helpMatrix = *Umatrix;
  //cout << "Umatrix before multiplications " << *Umatrix << endl;
  for(int a=0; a < noOfTimesteps -1; a++){
    *helpMatrix = mult(*Umatrix, *helpMatrix);
  }// for a
  
  *Umatrix = *helpMatrix;
  //cout << "Umatrix after multiplications " << *Umatrix << endl;


}


