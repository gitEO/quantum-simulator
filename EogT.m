function [Emax, deltaT]=EogT(top, bottom)
% Here we calculate delta t and Emax from desired start and end of spectrum
% Syntax [Emax, deltaT]=EogT( topOfSpectrum,  bottomOfSpectrum);
    Emax = top;
    Erange= top - bottom    Emin = Emax - Erange;
    deltaT=2*pi/Erange;
    
    
   
