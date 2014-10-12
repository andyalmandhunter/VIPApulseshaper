% N_FUSED_SILICA   Calculate the index of refraction for fused
%    silica, given the wavelength, using Sellmeier coefficients.
%
% A. Almand-Hunter, Thu Apr 24 22:21:37 MDT 2014
function n = n_fused_silica(lambda)

B1 = 0.696166300;
B2 = 0.407942600;
B3 = 0.897479400;

% C coefficients in um^2
C1 = 4.67914826e-3;
C2 = 1.35120631e-2;
C3 = 97.9340025;

lambda = lambda*1e-3;
n = sqrt(1 + B1*lambda^2/(lambda^2-C1) + ...
         B2*lambda^2/(lambda^2-C2) + ...
         B3*lambda^2/(lambda^2-C3));

end
