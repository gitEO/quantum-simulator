# Quantum computer simulator
## Intro

This is a simulation of a quantum computer simulating many-body quantum physics.

The model has been moved to github to share with students and researchers around the world. Any contributions will of course also be highly appreciated.

For more information see the
[paper](http://arxiv.org/abs/0705.1928v1)  discussing the model or my PhD 
 [thesis](http://folk.uio.no/ovrum/PhDthesis.pdf) with an introduction to quantum simulators as well.


### The physics of the model
From my thesis: "We used the Jordan-Wigner transformation to develop a program
that can take any two-body fermionic Hamiltonian and output the qubit gates a
quantum computer needs to simulate the time evolution of that Hamiltonian. We
then incorporated this into a quantum computer simulator, to show how this
general fermionic simulator finds the eigenvalues of the Hubbard and pairing models."

What the code does is two-fold:

1. The first is to run a "quantum simulator compiler", that takes any two-body fermionic Hamiltonian and outputs a quantum simulator "code", in terms of two-qubit gates needed to perform the simulation of that Hamiltonian.

2. The next part is to simulate a quantum simulator (special purpose quantum computer) that uses the qubit gates found in the previous step to find the eigenvalues of the two-body Hamiltonian using the phase-estimation algortithm.


The code could apropriately be used by students wishing to understand quantum computers and quantum simulators especially.


### Programming languages
The main code is written in C++, using the [blitz++](https://en.wikipedia.org/wiki/Blitz%2B%2B) library for vectors and matrices. There are also additional analysis scripts written using matlab in .m files, these can all be run by using [octave](http://www.gnu.org/software/octave/).

### File descriptions

+ main.cpp: Main code bulk.
+ main.h: Corresponding library file.
+ ny.cpp: New function definitions.
+ ny.h: Corresponding library file.
+ blitz2mat.cpp: Converts the output files from blitz format to octave readable.
+ analysis.m: Performs analysis of simulation, calls fig.m.
+ fig.m: Plots the resulting energy spectrum.
+ EogT.m: Calculates delta t and Emax from desired start and end of spectrum.

### Compile commands
See Makefile.



