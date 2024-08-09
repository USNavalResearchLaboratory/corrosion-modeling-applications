function [zmod,zr,zi,zp,goodCalc,outputFitParams] = matrixImpedance_ExtendedModifiedREAPPlain(b,f)
%matrixImpedance_ExtendedModifiedREAP - calculates the impedance of an
%equivalent circuit
%
% This function calculates the equivalent impedance of the extended
% modified REAP equivalent circuit.
%
% Syntax:  [zmod,zr,zi,zp,goodCalc,outputFitParams] = 
% matrixImpedance_ExtendedModifiedREAP(b,f)
%
% Inputs:
%   b = vector of adjustable parameters
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
% Other m-files required: solutionConductivity, ZwB, Zcpe
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
    RTopcoat = b(2);
    Y0topcoat = b(3);
    aTopcoat = b(4);
    sC_Topcoat = b(5);
    Deff_Topcoat = b(6);
    wetTopcoatDepth = b(7);
    RPrimer = b(8);
    Y0primer = b(9);
    aPrimer = b(10);
    sC_Primer = b(11);
    Deff_Primer = b(12);
    wetPrimerDepth = b(13);
    Rpre = b(14);
    Y0oxide = b(15);
    aOxide = b(16);
    %==============================
    %==============================
    nf = numel(f);
    nzFs = f > 0.0;
    fCorr = f(nzFs);
    n = numel(fCorr); 
    iSum = zeros(size(fCorr));   
    goodCalc = true;
    %==============================
    %==============================
    vBase = 1.0;
    mSize = 4;           
    A = zeros(mSize,mSize,n);    
    v(1:mSize) = vBase;
    %==============================
    %==============================
    [~,zcpeTR,zcpeTI,~] = Zcpe(fCorr,aTopcoat,Y0topcoat); 
    [~,zcpePR,zcpePI,~] = Zcpe(fCorr,aPrimer,Y0primer);   
    [~,zcpePreR,zcpePreI,~] = Zcpe(fCorr,aOxide,Y0oxide);

    flag = 0;
    
    if flag == 0
        [~,zWTopcoatR,zWTopcoatI,~] = ZwB(fCorr,sC_Topcoat,Deff_Topcoat,wetTopcoatDepth);
    elseif flag == 1
        [~,zWTopcoatR,zWTopcoatI,~] = ZwR(fCorr,sC_Topcoat,Deff_Topcoat,wetTopcoatDepth);
    elseif flag == 2
        [~,zWTopcoatR,zWTopcoatI,~] = Zw(fCorr,sC_Topcoat);
    end

    if flag == 0
        [~,zWPrimerR,zWPrimerI,~] = ZwB(fCorr,sC_Primer,Deff_Primer,wetPrimerDepth);
    elseif flag == 1
        [~,zWPrimerR,zWPrimerI,~] = ZwR(fCorr,sC_Primer,Deff_Primer,wetPrimerDepth);
    elseif flag == 2
        [~,zWPrimerR,zWPrimerI,~] = Zw(fCorr,sC_Primer);        
    end
    
    A(1,4,1:n) = Rs0;    
    A(2,4,1:n) = Rs0;    
    A(3,4,1:n) = Rs0;
    A(4,1,1:n) = Rs0;
    A(4,2,1:n) = Rs0;
    A(4,3,1:n) = Rs0;

    A(1,1,1:n) = complex(Rs0 + RTopcoat + zWTopcoatR(1:n) + RPrimer + zWPrimerR(1:n) + Rpre, zWTopcoatI(1:n) + zWPrimerI(1:n));
    A(1,2,1:n) = complex(Rs0 + RTopcoat + zWTopcoatR(1:n) + RPrimer + zWPrimerR(1:n), zWTopcoatI(1:n) + zWPrimerI(1:n));        
    A(1,3,1:n) = complex(Rs0 + RTopcoat + zWTopcoatR(1:n), zWTopcoatI(1:n));
    
    A(2,1,1:n) = A(1,2,1:n);
    A(2,2,1:n) = complex(Rs0 + RTopcoat + zWTopcoatR(1:n) + RPrimer + zWPrimerR(1:n) + zcpePreR(1:n), zWTopcoatI(1:n) + zWPrimerI(1:n) + zcpePreI(1:n));
    A(2,3,1:n) = A(1,3,1:n);
    
    A(3,1,1:n) = A(1,3,1:n);
    A(3,2,1:n) = A(2,3,1:n);
    A(3,3,1:n) = complex(Rs0 + RTopcoat + zWTopcoatR(1:n) + zcpePR(1:n), zWTopcoatI(1:n) + zcpePI(1:n));

    A(4,4,1:n) = complex(Rs0 + zcpeTR(1:n), zcpeTI(1:n)); 

    parfor i = 1:n
        Ainv2 = pinv(A(:,:,i));
        iM =  Ainv2 * v'; %v/A(:,:,i)
        iSum(i) = (iM(1) + iM(2) + iM(3) + iM(4));
    end

    zt = vBase./iSum;   
    zr = real(zt);
    zi = imag(zt);
    zp = atan(zi./zr).*(180.0/pi);      
    zpc = conj(zt);
    zmod = sqrt(zt.*zpc);      
    outputFitParams = [Rs0, RTopcoat, Y0topcoat, aTopcoat, sC_Topcoat, Deff_Topcoat, wetTopcoatDepth, RPrimer, Y0primer, aPrimer, sC_Primer,Deff_Primer,wetPrimerDepth, Rpre, Y0oxide, aOxide];
    
    if n < nf
        zr(n+1:nf) = 0.0;
        zi(n+1:nf) = 0.0;
        zp(n+1:nf) = 0.0;
        zmod(n+1:nf) = 0.0;
    end
end
    % for i = 1:n
    %     val1 = zcpePreR(i);
    %     val2 = zcpePreI(i);
    %     val3 = zcpePR(i);
    %     val4 = zcpePI(i);
    %     val5 = zcpeTR(i);
    %     val6 = zcpeTI(i);
    % 
    %     if  isreal(val1) && isreal(val2) && isreal(val3) && isreal(val4) && isreal(val5) && isreal(val6)
    % 
    %                 A(1,1) = complex(Rs0 + RTopcoat + zWTopcoatR(i) + RPrimer + zWPrimerR(i) + Rpre, zWTopcoatI(i) + zWPrimerI(i));
    %                 A(1,2) = complex(Rs0 + RTopcoat + zWTopcoatR(i) + RPrimer + zWPrimerR(i), zWTopcoatI(i) + zWPrimerI(i));        
    %                 A(1,3) = complex(Rs0 + RTopcoat + zWTopcoatR(i), zWTopcoatI(i));
    % 
    %                 A(2,1) = A(1,2);
    %                 A(2,2) = complex(Rs0 + RTopcoat + zWTopcoatR(i) + RPrimer + zWPrimerR(i) + zcpePreR(i), zWTopcoatI(i) + zWPrimerI(i) + zcpePreI(i));
    %                 A(2,3) = A(1,3);
    % 
    %                 A(3,1) = A(1,3);
    %                 A(3,2) = A(2,3);
    %                 A(3,3) = complex(Rs0 + RTopcoat + zWTopcoatR(i) + zcpePR(i), zWTopcoatI(i) + zcpePI(i));
    % 
    %                 A(4,4) = complex(Rs0 + zcpeTR(i), zcpeTI(i));   
    % 
    %                 iM = pinv(A) * v; %linsolve(A,v);
    %                 iSum = (iM(1) + iM(2) + iM(3) + iM(4));      
    %                 % iSumc = conj(iSum);
    %                 % iMod = sqrt(iSum*iSumc);
    %                 % iT(i) = iMod;
    % 
    %                 zt(i) = vBase/iSum;
    %                 % zr(i) = real(zt);
    %                 % zi(i) = imag(zt);
    %                 % zp(i) = atan(zi(i)/zr(i))*(180/3.14159265);                      
    %     else
    %         zt(i) = 0.0;
    %          % zr(i) = 0.0;
    %          % zi(i) = 0.0;
    %          % zp(i) = 0.0;
    %          goodCalc = false;
    %     end      
    % end