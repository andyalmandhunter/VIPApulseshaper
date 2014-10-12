% RS Calculate the discrete Rayleigh-Sommerfeld integral.
%
%    The Rayleigh-Sommerfeld integral is a convolution integral
%    between the initial field U and a kernel H calculated from the
%    impulse response of the Helmholtz equation.  Use FFT to do
%    this integral efficiently (See Shen and Wang [2006]).
%
%    Arguments: Four arrays for x and y values of initial and result
%               matrices.  Matrix of the initial field.  Propagation
%               distance z.  Angular wavenumber k.  Physical grid
%               spacing ds and dn.  Angle theta of final plane.
%               Vertical position y0 of rotation axis.
%
%    Return:    Matrix: final scalar field.
%
% A. Hunter, Tue Sep 10 22:25:49 MDT 2013
function E1 = rs(E0, ds, dn, Si, Ni, Xi, Yi, z, k, theta, y0)

    % Weight matrix for Simpson's Rule
    w = 3 .\ [1 repmat([4 2],1,499) 4 1];
    [w1, w2] = meshgrid(w, w);
    W = w1 .* w2;
% $$$     E0 = E0 .* W;

    % Calculate some necessary constants, and initial and final
    % matrices
    N = length(Si);
% $$$     [n s] = meshgrid(Si, Ni);
% $$$     [y x] = meshgrid(Xi, Yi);
    Xj = [Xi(1) - Si(N + 1 - (1:(N-1))) Xi((N:(2*N-1)) - N + 1) - Si(1)];
    Yj = [Yi(1) - Ni(N + 1 - (1:(N-1))) Yi((N:(2*N-1)) - N + 1) - Ni(1)];
    [Y, X] = meshgrid(Xj, Yj);

    % Construct initial matrix U and convolution seed H.  Zero
    % padded
    U = [[E0 zeros(length(Yi), length(Xi)-1)]; ...
         zeros(length(Yi)-1, 2*length(Xi)-1)];
    Z = z + (Y - y0) * tan(theta);
    r = sqrt(X.^2 + Y.^2 + Z.^2);
    H = 1/(2*pi) * exp(1i*k*r)./r .* Z./r .* (1./r - 1i*k);

    % Calculate integral
    S = ifft2(fft2(U) .* fft2(H)) * ds * dn;

    % Translate to result coordinates
    E1 = S(N:end,N:end);
    
end
