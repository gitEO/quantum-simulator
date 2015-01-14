function [Emax, deltaT]=EogT(top, bottom)
%Her regner vi ut delta t og Emax utifra ønska start 
%og slutt på spekteret. 
%Syntax [Emax, deltaT]=EogT( toppenAvSpekteret,  bånnAvSpekteret);
    Emax = top;
    Erange= top - bottom
    Emin = Emax - Erange;
    deltaT=2*pi/Erange;
    
    
   