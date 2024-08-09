function [zmod,zr,zi,zp,goodCalc,outputFitParams] = PolynomialTestPlain(b,f)
%matrixImpedance_ModifiedUndamagedCoating - calculates the impedance of an
%equivalent circuit
%
% This function calculates the equivalent impedance of the modified
% undamaged coating equivalent circuit.
%
% Syntax:  [zmod,zr,zi,zp,goodCalc,outputFitParams] = 
% matrixImpedance_ModifiedUndamagedCoating(b,f)
%
% Inputs:
%   b = vector adjustable parameters
%   f = vector of frequencies
%
% Outputs: 
%   zmod = vector containing the impedance modulus per frequency
%   zr = vector containing the real component of the impedance
%   zi = vector containing the imaginary component of the impedance
%   zp = vector containing the impedance phase
%   goodCalc = boolean variable to indicate convergence
%   outputFitParams = vector of the calculated circuit element parameters
%
%
% Other m-files required: solutionConductivity, Zcpe
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
% June 2022; Last revision: 13 August 2022
%==========================================================================
    %==============================
    c1 = b(1);
    c2 = b(2);
    %==============================
    %==============================
    iT = zeros(size(f));
    % zr = zeros(size(f));
    % zi = zeros(size(f));
    % zp = zeros(size(f));
    zt = zeros(size(f));
    goodCalc = true;
    %==============================
    
    n = numel(f);
    for i = 1:n
        zt(i) = c1 + c2*(f(i)^2);
    end
    zr = real(zt);
    zi = zeros(size(zr));
    zp = atan(zi./zr).*(180.0/pi);     
    % zmod = sqrt(zr.^2 + zi.^2); 
    % zpc = conj(zt);
    zmod = zr;    
    outputFitParams = [c1,c2];
end