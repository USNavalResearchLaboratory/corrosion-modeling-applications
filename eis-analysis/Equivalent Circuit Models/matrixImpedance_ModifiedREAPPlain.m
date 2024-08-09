function [zmod,zr,zi,zp,goodCalc,outputFitParams] = matrixImpedance_ModifiedREAPPlain(b,f)
%matrixImpedance_ModifiedREAP - calculates the impedance of an
%equivalent circuit
%
% This function calculates the equivalent impedance of the modified REAP
% equivalent circuit.
%
% Syntax:  [zmod,zr,zi,zp,goodCalc,outputFitParams] = 
% matrixImpedance_ModifiedREAP(b,f)
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
% June 2022; Last revision: 24 June 2022
%==========================================================================
    %==============================
    Rs0 = b(1);
    R1 = b(2);
    Y0c = b(3);
    aC = b(4);
    Rpre = b(5);
    Y0oxide = b(6);
    aOxide = b(7);
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

    [~,zcpeDLR,zcpeDLI,~] = Zcpe(f,aOxide,Y0oxide);   
    [~,zcpeCR,zcpeCI,~] = Zcpe(f,aC,Y0c);  
    
    n = numel(f);
    
    A(1,3) = Rs0;    
    A(3,1) = Rs0;   
    A(2,3) = Rs0;
    A(3,2) = Rs0;

    A(1,1) = Rs0 + R1 + Rpre;
    A(1,2) = Rs0 + R1;
    A(2,1) = Rs0 + R1;

    for i = 1:n
        A(2,2) = complex(Rs0 + R1 + zcpeDLR(i), zcpeDLI(i));
        A(3,3) = complex(Rs0 + zcpeCR(i), zcpeCI(i));
       
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
    outputFitParams = [Rs0, R1, Y0c, aC, Rpre, Y0oxide, aOxide];
end