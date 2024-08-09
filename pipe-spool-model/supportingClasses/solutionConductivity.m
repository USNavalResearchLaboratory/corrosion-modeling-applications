function k = solutionConductivity(T,c)
% solutionConductivity - Calculates the conductivity of a solution
%
% Function to calculate the conductivity of an NaCl solution.
%
% Syntax:  k = solutionConductivity(T,c)
%
% Inputs: 
% c = solution concentration (M) 
% T = solution temperature (C)
%  
%
% Outputs: Conductivity (S/m)
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 
%
%==========================================================================
% Author:   Steve Policastro, Ph.D., Materials Science
% Center for Corrosion Science and Engineering, U.S. Naval Research
% Laboratory
% email address: steven.policastro@nrl.navy.mil  
% Website: 
% June 2022; Last revision: 24 June 2022
%==========================================================================

    a11 = 0.0480;
    a12 = 0.0034;
    a21 = 6.7545;
    a22 = 0.2392;
    a31 = 0.2065;
    a32 = 0.0013;

    num = (a11 + a12.*T) + ((a21 + a22.*T).*c);
    denom = 1.0 + ((a31 + a32.*T).*c);

    k = num./denom; %S/m
    
end