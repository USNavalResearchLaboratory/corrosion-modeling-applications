function [zmod,zr,zi,zp,goodCalc,outputFitParams] = matrixImpedance_REAPPlain(b,f)
%matrixImpedance_REAP - calculates the impedance of an
%equivalent circuit
%
% This function calculates the equivalent impedance of the REAP equivalent
% circuit.
%
% Syntax:  [zmod,zr,zi,zp,goodCalc,outputFitParams] = 
% matrixImpedance_REAP(b,f)
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
% Other m-files required: solutionConductivity, Zc, Zcpe
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
    Rs0 = b(1);
    R1 = b(2);
    Cc = b(3);
    Rp = b(4);
    Y0oxide = b(5);
    aOxide = b(6);
    %==============================
    %==============================
    iT = zeros(size(f));
    zr = zeros(size(f));
    zi = zeros(size(f));
    zp = zeros(size(f));
    goodCalc = true;
    %==============================
    %==============================
    A = zeros(3,3);
    vBase = 1.0;
    v = [vBase; vBase; vBase];

    [~,zCR,zCI,~] = Zc(f,Cc); 
    [~,zcpeDLR,zcpeDLI,~] = Zcpe(f,aOxide,Y0oxide);  
    
    n = numel(f);
    
    A(1,3) = Rs0;    
    A(3,1) = Rs0;   
    A(2,3) = Rs0;
    A(3,2) = Rs0;

    A(1,1) = Rs0 + R1 + Rp;
    A(1,2) = Rs0 + R1;
    A(2,1) = Rs0 + R1;

    for i = 1:n
        A(2,2) = complex(Rs0 + R1 + zcpeDLR(i), zcpeDLI(i));
        A(3,3) = complex(Rs0 + zCR(i), zCI(i));
       
        iM = linsolve(A,v);
        iSum = (iM(1) + iM(2) + iM(3));
        iSumc = conj(iSum);
        iMod = sqrt(iSum*iSumc);
        iT(i) = iMod;

        zt = vBase/iSum;
        zr(i) = real(zt);
        zi(i) = imag(zt);
        zp(i) = atan(zi(i)/zr(i))*(180.0/pi);       
    end
    zmod = sqrt(zr.^2 + zi.^2); 
    outputFitParams = [Rs0, R1, Cc, Rp, Y0oxide, aOxide];
end