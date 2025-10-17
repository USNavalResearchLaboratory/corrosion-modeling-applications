function [zm,zr,zi,zp] = Zw(f,sigma)
% Zw - Calculates the impedance of a Warburg element with semi-infinite
% diffusion to a large planar electrode
%
% Function to calculate the impedance of a bounded Warburg element.
%
% Syntax:  [zm,zr,zi,zp] = Zw(f,sigma)
%
% Inputs: 
% f = vector of frequencies.
% sigma = Warburg coefficient, aka "capacitance"
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
% email address: steven.a.policastro.civ@us.navy.mil  
% Website: 
% Created: August 2022
% Last revision: 23 July 2025
%==========================================================================
    w = (2*pi).*f;
    zr = zeros(size(w));
    zi =zeros(size(w));
    zp = zeros(size(w));
    ztemp = zeros(size(w));

    n = numel(w);
    for i = 1:n
        % Using the equations from Gamry for the "semi-infinite" Warburg
        % impedance
        % https://www.gamry.com/Framework%20Help/HTML5%20-%20Tripane%20-%20Audience%20A/Content/EIS/Theory/Physical%20Electrochemistry%20and%20Circuit%20Elements/Diffusion.htm
        term1 = sigma/sqrt(w(i));
        term2 = sigma/(1i*sqrt(w(i)));
        z = term1 + term2;
        
        zr(i) = real(z);
        zi(i) = imag(z);
        ztemp(i) = complex(zr(i),zi(i));
        zp(i) = atan(zi(i)/zr(i));
    end
    zpc = conj(ztemp);
    zm = sqrt(zp.*zpc); 
end