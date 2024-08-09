function [zm,zr,zi,zp] = Zcpe(f,a,Y0)
% Zcpe - Calculates the impedance of a constant phase element.
%
% Function to calculate the impedance of a constant phase element.
%
% Syntax:  [zm,zr,zi,zp] = Zcpe(f,a,Y0)
%
% Inputs: 
% f = vector of frequencies.
% a = exponent
% Y0 = CPE constant
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
        zr(i) = 1.0/((w(i)^a)*Y0) * cos((a*pi)/2.0);
        zi(i) = -1.0/((w(i)^a)*Y0) * sin((a*pi)/2.0);   
        if zr(i) == isreal(zr(i)) && zi(i) == isreal(zi(i))
            ztemp(i) = complex(zr(i),zi(i));
        else
            ztemp(i) = complex(-1.0,-1.0);
        end
        zp(i) = atan(-zi(i)/zr(i));
    end
    zpc = conj(ztemp);
    zm = sqrt(ztemp.*zpc);    
end