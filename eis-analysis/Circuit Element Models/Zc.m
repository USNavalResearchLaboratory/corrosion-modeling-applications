function [zm,zr,zi,zp] = Zc(f,C)
% Zc - Calculates the impedance of a capacitor.
%
% Function to calculate the impedance of a capacitor.
%
% Syntax:  [zm,zr,zi,zp] = Zc(f,C)
%
% Inputs: 
% f = vector of frequencies
% C = capacitance.
%  
%
% Outputs: 
% zm = vector with the modulus of impedance at each frequency
% zr = vector with the real component of impedance at each frequency
% zi = vector with the imaginary component of impedance at each frequency
% zp = vector with the phase of impedance at each frequency
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
    w = (2*pi).*f;
    zr = zeros(size(w));
    zi =zeros(size(w));
    zp = zeros(size(w));
    ztemp = zeros(size(w));

    n = numel(w);

    for i = 1:n
        zr(i) = 0.0;
        zi(i) = -1.0/(w(i)*C);    
        ztemp(i) = complex(zr(i),zi(i));
        zp(i) = atan(-zi(i)/zr(i));
    end
    zpc = conj(ztemp);
    zm = sqrt(ztemp.*zpc);     
end