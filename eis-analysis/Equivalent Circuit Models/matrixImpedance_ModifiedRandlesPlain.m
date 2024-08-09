function [zmod,zr,zi,zp,goodCalc,outputFitParams] = matrixImpedance_ModifiedRandlesPlain(b,f)
%matrixImpedance_ModifiedRandles - calculates the impedance of an
%equivalent circuit
%
% This function calculates the equivalent impedance of the modified Randles
% equivalent circuit.
%
% Syntax:  [zmod,zr,zi,zp,goodCalc,outputFitParams] = 
% matrixImpedance_ModifiedRandles(b,f,physicalParams)
%
% Inputs:
%   b = vector adjustable parameters
%   f = vector of frequencies
%   physicalParams = vector of unchanging parameters
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
% Other m-files required: solutionConductivity, Zcpe, ZwB
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
    Y0 = b(3);
    alpha = b(4);
    sigma = b(5);
    B = b(6);
    % Deff = b(6);
    % diffusionPathLength = b(7);
    %==============================
    %==============================
    % n = numel(f);
    % iT = zeros(size(f));
    zt = zeros(size(f));  
    iSum = zeros(size(f));    
    goodCalc = true;
    %==============================
    %==============================
    vBase = 1.0;
    mSize = 2;
    n = numel(f);    
    A = zeros(mSize,mSize); %,n
    v = [vBase;vBase];
    %==============================
    %==============================
 
    [~,zCPETR,zCPETI,~] = Zcpe(f,alpha,Y0);  
    
    flag = 0;    
    if flag == 0
        [~,zWTR,zWTI,~] = ZwB(f,sigma,B); %Deff,diffusionPathLength
    elseif flag == 1
        [~,zWTR,zWTI,~] = ZwR(f,sigma,B); %Deff,diffusionPathLength
    elseif flag == 2
        [~,zWTR,zWTI,~] = Zw(f,sigma);
    end

    A(1,2) = Rs0;
    A(2,1) = Rs0;
    for i = 1:n
        A(1,1) = complex(Rs0 + R1 + zWTR(i), zWTI(i));
        A(2,2) = complex(Rs0 + zCPETR(i),zCPETI(i));

        iM = linsolve(A,v);
        iSum = (iM(1) + iM(2));       
        zt(i) = vBase/iSum;      
    end
    % zt = vBase./iSum;
    zr = real(zt);
    zi = imag(zt);
    zp = atan(zi./zr).*(180.0/pi);      
    zpc = conj(zt);
    zmod = sqrt(zt.*zpc);       
    outputFitParams = [Rs0,R1,Y0,alpha,sigma,B]; %Deff,diffusionPathLength
end