% VIPAPS   Propagate light waves through the VIPA pulse shaper, so
%    that I can adjust parameters (positions of optics) in the
%    computer before I attempt optimization on the optics table.
%
%    For the initial attempt, solve this exactly, without paraxial
%    approximation.  If this is too hard, okay to assume small angles.
%
% A. Hunter, Fri Aug 30 22:34:26 MDT 2013
%
%    I can solve this exactly ... see FFT-DI method described in Shen
%    and Wang (2006).
%
% A. Hunter, Fri Sep  6 06:56:20 MDT 2013
%
%    Extended means beginning and end frames are bigger (so I can
%    ultimately include more virtual images to see what that does to
%    the system).
%
%    NOTE -- now the pixel spacing may be below the sampling limit.
%
%    Yeah, there are ailasing errors ... need to bump up resolution
%
%    Mode from 0-31
%
% A. Almand-Hunter, Mon Apr  1 2014
%
%    Back to only using the middle 1cm by 1cm piece for the
%    VIPA-to-lens step (for speed and to match reality).
%
%    (Checked -- see dispersion20 run today -- and this changes the
%    result hardly at all ... much less than going to the paraxial
%    approximation)
%
% A. Almand-Hunter, Fri Apr 25 15:55:58 MDT 2014



% $$$ function E7 = vipaps(mode)
clear;
mode = 0;


% Simulation parameters  ---------------------------------------------
%
%    all distances in um
%
% --------------------------------------------------------------------

% VIPA parameters
thickness = 4130;
height = 10000;
thetai = 1.21 * pi/180;  % Angle between VIPA normal and optical axis
T = 1 - 0.985;  % Transmittance of front face
foci = 1:56;  % Just specify this, rather than calculating it from the
              % height of the VIPA and the angle.
n = n_fused_silica(815);  % Index of refraction of VIPA material
                          % (fused silica at 815 nm)

% dy and dz between VIPA foci, calculated from above parameters
dy = 2*thickness * tan(thetai);
dz = 2*n*thickness / cos(thetai);

nfoci = floor(height*cos(thetai)/dy);
foci = 0:nfoci;

% Pulse-shaper parameters (choice of w0y and w0x is *influenced* by
% the VIPA geometry, but can be adjusted independently).
f = 750e3;  % 75 cm Fourier-transform lens
% $$$ w0y = 27.1;  % Chosen to maximize light coupling to VIPA
w0y = 27.1;  % Chosen to maximize light coupling to VIPA
w0x = 2350;

% Laser parameters
lambda = 0.815;
k = 2*pi/lambda + (mode * 1/32) * 2*pi/dz;  % 13.7 mm in z between foci


% Dimensions of physical space to consider  --------------------------
%
%    (S,N) are coordinates in initial plane; (X,Y) are coordinates in
%    final plane.  Spacing in um.  Window for doing convolution steps
%    is 1001x1001 elements.
% --------------------------------------------------------------------

% $$$ ds = 5;
% $$$ dn = ds;
% $$$ Si = -5000:ds:5005;
% $$$ Ni = -40035:dn:40040;
% $$$ [s n] = meshgrid(Si, Ni);

% $$$ % Blocks for convolutions
% $$$ xn = 0:1;
% $$$ yn = 0:15;
% $$$ bs = 1001;  % Block size

ds = 10;
dn = ds;
Si = -5000:ds:5000;
Ni = -15010:dn:15010;
[s n] = meshgrid(Si, Ni);

% Blocks for convolutions
xn = 0;
yn = 0:2;
bs = 1001;  % Block size

% Preallocate arrays;
Efinal = zeros(1, length(foci));
E0f = zeros(length(Ni), length(Si));
E1f = zeros(length(Ni), length(Si));
E2f = zeros(length(Ni), length(Si));
E3f = zeros(length(Ni), length(Si));
E6f = zeros(length(Ni), length(Si));


% --------------------------------------------------------------------
% Find field at FT lens
%
%    Iterate through all of the virtual beam focuses created by the
%    VIPA.  Add all of the resulting fields in the plane of the
%    Fourier-transform lens.  Beam parameters w0x and w0y at the
%    focuses are given above.  For now, assume vacuum external and
%    internal to the VIPA.  If necessary, I can modify this later.
%
%    In order to get a large enough window at the lens and SLM,
%    this must be done in pieces, then stitched together. 
% --------------------------------------------------------------------

for j = 1  % Iterate through initial foci
    
    E0 = zeros(length(Ni), length(Si));
    E1 = zeros(length(Ni), length(Si));
    L = zeros(length(Ni), length(Si));

    % Calculate propagation distance z and vertical focus position y0
    middle = nfoci/2;
    z = f - dz*middle + dz*j;
    y0 = -1*dy*middle + dy*j;
    
    % Construct initial field
    E0 = exp(-1*(s.^2) / w0x^2 - ((n - y0).^2) / w0y^2);
    
    E0f = E0f + E0;

    fprintf(1, 'Propagating to lens from focus %d of %d.\n', j, ...
            length(foci));
    
    % Iterate through (x,y) blocks (final plane)
    for xc = xn
        for yc = yn
            
            % Iterate through (s,n) blocks (initial plane)
            for sc = xn
                for nc = yn
                    
                    % Calculate index ranges
                    xr = (1:bs) + bs*xc;
                    yr = (1:bs) + bs*yc;
                    sr = (1:bs) + bs*sc;
                    nr = (1:bs) + bs*nc;

                    % Perform this part of the convolution
                    E1(yr,xr) = E1(yr,xr) + rs(E0(nr,sr), ds, dn, ...
                                               Si(sr), Ni(nr), Si(xr), ...
                                               Ni(yr), z, k, 0, 0);
                    fprintf(1, '(%d,%d) --> (%d,%d); ', sc,nc,xc,yc);
                    
                end
            end
            
        end
    end

    E0 = E1;
    E1 = zeros(length(Ni), length(Si));

    fprintf(1, '\n');
    
    E1f = E1f + E0;

    
    % --------------------------------------------------------------------
    % Propagate through an ideal thin lens
    %
    %    Note: Allow lens to have a y offset from the original optical
    %    axis.
    %
    %    Lens phase factor: exp[-i pi (x^2 + y^2) / (lambda f)]
    % --------------------------------------------------------------------

    % Distance above optical axis in um
    y0L = 0;

    % Lens phase
    L = exp(-1i * pi * (s.^2 + (n - y0L).^2) / (lambda * f));

    % Apply
    E0 = E0 .* L;

    E2f = E2f + E0;
    

    % --------------------------------------------------------------------
    % Propagate to Fourier plane of pulse shaper (SLM)
    %
    %     Allow for the possibility that the SLM is tilted at a small
    %     angle.  To be more exact, this means that I need to replace the
    %     term z/r in H with (n {dot} r)/r (see Jackson Eq. 10.85).  n is
    %     a normal to the plane of the initial field; thus (n {dot} r)
    %     is often just z, but now z can depend on x and y.
    % --------------------------------------------------------------------

    fprintf(1, 'Propagating to SLM.\n');

    % Propagation distance in um
    z = 7.5e5;

    % Angle of SLM
    theta = 0;

    % Iterate through (x,y) blocks (final plane)
    for xc = xn
        for yc = yn
            
            % Iterate through (s,n) blocks (initial plane)
            for sc = xn
                for nc = yn
                    
                    % Calculate index ranges
                    xr = (1:bs) + bs*xc;
                    yr = (1:bs) + bs*yc;
                    sr = (1:bs) + bs*sc;
                    nr = (1:bs) + bs*nc;

                    % Perform this part of the convolution
                    E1(yr,xr) = E1(yr,xr) + rs(E0(nr,sr), ds, dn, Si(sr), ...
                                               Ni(nr), Si(xr), Ni(yr), ...
                                               z, k, theta, y0L);
                    fprintf(1, '(%d,%d) --> (%d,%d); ', sc,nc,xc,yc);
                    
                end
            end
            
        end
    end
    
    E0 = E1;
    E1 = zeros(length(Ni), length(Si));
    
    
    fprintf(1, '\n');
    
    E3f = E3f + E0;
    
    % Save this electric field:
    outfile = sprintf(['/home/andy/Desktop/VIPA_Simulation/dispersion110/' ...
                       'Eslm%dreal.dat'], mode);
    out = real(E3);
    save(outfile, 'out', '-ascii', '-double');
    outfile = sprintf(['/home/andy/Desktop/VIPA_Simulation/dispersion110/' ...
                       'Eslm%dimag.dat'], mode);
    out = imag(E3);
    save(outfile, 'out', '-ascii', '-double');


    % --------------------------------------------------------------------
    % For testing, chop off top and bottom 500 rows of E3
    % --------------------------------------------------------------------

% $$$ E3(1:500,:) = 0;
% $$$ E3(end-500:end,:) = 0;


    % --------------------------------------------------------------------
    % Propagate back to lens
    %
    %    Again need to account for the fact that the SLM may be
    %    tilted.  This is the same as the previous step.
    % --------------------------------------------------------------------

    fprintf(1, 'Propagating back to lens.\n');

    % Propagation distance in um
    z = 7.5e5;

    % Angle of SLM
    theta = 0;

    % Iterate through (x,y) blocks (final plane)
    for xc = xn
        for yc = yn
            
            % Iterate through (s,n) blocks (initial plane)
            for sc = xn
                for nc = yn
                    
                    % Calculate index ranges
                    xr = (1:bs) + bs*xc;
                    yr = (1:bs) + bs*yc;
                    sr = (1:bs) + bs*sc;
                    nr = (1:bs) + bs*nc;

                    % Perform this part of the convolution
                    E1(yr,xr) = E1(yr,xr) + rs(E0(nr,sr), ds, dn, Si(sr), ...
                                               Ni(nr), Si(xr), Ni(yr), ...
                                               z, k, theta, y0L);
                    fprintf(1, '(%d,%d) --> (%d,%d); ', sc,nc,xc,yc);
                    
                end
            end
            
        end
    end

    E0 = E1;
    E1 = zeros(length(Ni), length(Si));

    fprintf(1, '\n');
    

    % --------------------------------------------------------------------
    % Propagate back throgh ideal thin lens
    % --------------------------------------------------------------------

    % Apply
    E0 = E0 .* L;


    % --------------------------------------------------------------------
    % Propagate back to virtual focuses
    %
    %    Again, do this piecewise, at the locations of the original
    %    virtual focuses.
    %
    %    MODIFY: Only need to calculate middle piece.  Iterate through
    %    and weight according to initial fields for each virtual
    %    focus.  I think I can call this "the field radiated from the
    %    pulse shaper back into the input mode."
    %
    %    As mentioned below, the next improvement should be to filter
    %    out just the k=0 part (think about this ... possibly the
    %    convolution with the Gaussian already does this).
    %
    %    DONE -- Just integrate over E instead of |E|^2.
    % --------------------------------------------------------------------

    for m = 1
        
        % Calculate propagation distance z and vertical focus position y0
        z = f - dz*middle + dz*m;
        
        fprintf(1, 'Propagating from lens to focus %d of %d.\n', m, length(foci));
        
        % Iterate through (x,y) blocks (final plane)
        for xc = xn
            for yc = yn
                
                % Iterate through (s,n) blocks (initial plane)
                for sc = xn
                    for nc = yn
                        
                        % Calculate index ranges
                        xr = (1:bs) + bs*xc;
                        yr = (1:bs) + bs*yc;
                        sr = (1:bs) + bs*sc;
                        nr = (1:bs) + bs*nc;

                        % Perform this part of the convolution
                        E1(yr,xr) = E1(yr,xr) + rs(E0(nr,sr), ds, dn, ...
                                                   Si(sr), Ni(nr), Si(xr), ...
                                                   Ni(yr), z, k, 0, 0);
                        fprintf(1, '(%d,%d) --> (%d,%d);', sc,nc,xc,yc);
                        
                    end
                end
                
            end
        end

        E0 = E1;
        E1 = zeros(length(Ni), length(Si));
        
        fprintf(1, '\n');
        
        % Filter field to input mode
        y0 = -1*dy*28.5 + dy*m;
        E1 = exp(-1*(s.^2) / w0x^2 - ((n - y0).^2) / w0y^2);
        Efinal(j) = Efinal(j) + sum(sum(E1 .* E0));
        
        E6f = E6f + E0;

    % Save final field, spatially filtered so I can later take
    % Fourier transforms to get different k components:
    outfile = sprintf(['/home/andy/Desktop/VIPA_Simulation/dispersion110/' ...
                       'Efinal%dfocus%dreal.dat'], mode, j);
    out = real(E0 .* E6);
    save(outfile, 'out', '-ascii', '-double');
    outfile = sprintf(['/home/andy/Desktop/VIPA_Simulation/dispersion110/' ...
                       'Efinal%dfocus%dimag.dat'], mode, j);
    out = imag(E0 .* E6);
    save(outfile, 'out', '-ascii', '-double');
        
    end

end


% Save final field  --------------------------------------------------

outfile = sprintf(['/home/andy/Desktop/VIPA_Simulation/dispersion200/' ...
                   'Efinal%dreal.dat'], mode);
out = real(E3f);
save(outfile, 'out', '-ascii', '-double');
outfile = sprintf(['/home/andy/Desktop/VIPA_Simulation/dispersion200/' ...
                   'Efinal%dimag.dat'], mode);
out = imag(E3f);
save(outfile, 'out', '-ascii', '-double');

outfile = sprintf(['/home/andy/Desktop/VIPA_Simulation/dispersion110/' ...
                   'Efinal%dreal.dat'], mode);
out = real(Efinal);
save(outfile, 'out', '-ascii', '-double');
outfile = sprintf(['/home/andy/Desktop/VIPA_Simulation/dispersion110/' ...
                   'Efinal%dimag.dat'], mode);
out = imag(Efinal);
save(outfile, 'out', '-ascii', '-double');

E7 = sum(Efinal);
