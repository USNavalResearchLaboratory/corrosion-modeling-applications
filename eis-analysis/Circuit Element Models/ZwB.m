function [zm,zr,zi,zp] = ZwB(f,s,B8) %D0,del
% ZwB - Calculates the impedance of a transmissive Warburg element.
%
% Function to calculate the impedance of a transmissive Warburg element.
%
% Syntax:  [zm,zr,zi,zp] = ZwB(x,s,D0,del)
%
% Inputs: 
% f = vector of frequencies.
% s = Warburg "capacitance"
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
% June 2022; Last revision: 24 June 2022
%==========================================================================
    w = (2*pi).*f;
    zr = zeros(size(w));
    zi =zeros(size(w));
    zp = zeros(size(w));
    ztemp = zeros(size(w));

    n = numel(w);
    csub = sqrt(1i);

    for i = 1:n    
        % a = s/sqrt(w(i));
        % b = -s/sqrt(w(i));
        % 
        % ff = ((del*sqrt(w(i)))/(sqrt(D0))) * csub;
        % 
        % c = tanh(ff);
        % 
        % d = complex(a,b);
        % e = c * d;
        % Using the Wikipedia equations! https://en.wikipedia.org/wiki/Warburg_element
        a = (csub*sqrt(w(i)));
        b = s/a;
        c = B8; %del/D0;
        d = tanh(c*a);
        e = b*d;
        zr(i) = real(e);
        zi(i) = imag(e);

        ztemp(i) = complex(zr(i),zi(i));
        zp(i) = atan(zi(i)/zr(i));
    end
    zpc = conj(ztemp);
    zm = sqrt(zp.*zpc); 
%     zm = vBase./iT;
end