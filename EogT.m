function [Emax, deltaT]=EogT(top, bottom)
%Her regner vi ut delta t og Emax utifra �nska start 
%og slutt p� spekteret. 
%Syntax [Emax, deltaT]=EogT( toppenAvSpekteret,  b�nnAvSpekteret);
    Emax = top;
    Erange= top - bottom
    Emin = Emax - Erange;
    deltaT=2*pi/Erange;
    
    
   