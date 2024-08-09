function MSE = objectiveFunction(yRD, yID, yRM, yIM, k)
% objectiveFunction - Calculates the MSE between a fit and data.
%
% Function to calculate the impedance of a capacitor.
%
% Syntax:  j = objectiveFunction(yRD, yID, yRM, yIM, k)
%
% Inputs: 
% yRD = vector of real data 
% yID = vector of imaginary data 
% yRM = vector of real model results
% yIM = vector of imaginary model results.
% k = number of fit parameters.
%  
%
% Outputs: j = scalar value of the MSE result
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
    nzyRD = yRD > 0.0;
    N = length(yRD(nzyRD));
    
    rError = ((yRD(1:N) - yRM(1:N))./yRD(1:N)).^2; 

    iError = zeros(size(rError));
    testIerror = yID(1:N) - yIM(1:N); 
    naniE = abs(testIerror) < 1.0e-100;
    iError(naniE) = 0.0;
    iError(~naniE) = (testIerror(~naniE)./yID(~naniE)).^2;

    errSum = sum(rError + iError);
    MSE = sqrt(errSum)/(N - k);  
end
