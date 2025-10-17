function [zm,zr,zi,zp] = ZwR(f,sigma,B)
% ZwB - Calculates the impedance of a finite length Warburg element with a
% reflective boundary
%
% Function to calculate the impedance of a reflective Warburg element.
%
% Syntax:  [zm,zr,zi,zp] = ZwR(f,sigma,B)
%
% Inputs: 
% f = vector of frequencies.
% sigma = Warburg "capacitance"
% B = del/sqrt(D0)
% D0 = ionic diffusion coefficeient (cm^2/s)
% del = diffusion length (cm).
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
% August 2022; 
% Revision: 23 August 2022
% Last revision: 23 July 2025
%==========================================================================
    w = (2*pi).*f;
    zr = zeros(size(w));
    zi =zeros(size(w));
    zp = zeros(size(w));
    ztemp = zeros(size(w));

    n = numel(w);
    csub = sqrt(1i);

    for i = 1:n
        % Using the Wikipedia equations! https://en.wikipedia.org/wiki/Warburg_element
        a = (csub*sqrt(w(i)));
        b = sigma/a;
        c = B; %del/D0;
        d = coth(c*a);
        z = b*d;        
        
        zr(i) = real(z);
        zi(i) = imag(z);

        ztemp(i) = complex(zr(i),zi(i));
        zp(i) = atan(zi(i)/zr(i));
    end
    zpc = conj(ztemp);
    zm = sqrt(zp.*zpc); 
end