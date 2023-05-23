%% Integrate mie calculation for a single radii over a size distribution

% INPUTS:

%   (1) r_modal - the modal radius of the droplet distribution. This is
%   NOT the ratio of the third moment of the droplet distribution to
%   the second moment. This is the most common radius observed in a large
%   random sampling of n(r), the size distribution. It is used to define
%   the droplet distributions most likely observed outcome.

%   (2) dist_var - the variance of the droplet distribution. A typical
%   value used for liquid water clouds is 7

%   (3) wavelength - (nanometers) this is the wavelength range used to calculate
%   the mie properties across the size distribution. If wavelength is a
%   single value, then a monochromatic calculation is performed. If
%   wavelength is a vector with 3 entries, then the first two are the boundaries
%   and the 3rd value is the discrete step between the two boundaries.

%   (4) index_of_refraction - this is the index of refraction used in the
%       scattering calculations, effectively telling the code which substance
%       to use. Options are:
%       (a) 'water' - this is an optional string input for libRadTran
%       (b) 'ice' - this is an optional string input for libRadTran
%       (c) a + bi - this numeric option should include a real and
%       imaginary component. One mie calculation will be run for each row
%       of numeric input - thus a new mie file for each substance

%   (5) size_distribution - this is a string that tells the code which type
%   of size distribution to integrate over. The options are:
%       (a) 'gamma' - gamma droplet distribution

% OUTPUTS:
% (1) ssa_avg - single scattering albedo averaged over all drop sizes
% within the distribution defined

% (2) Qe_avg - extinction coefficient averaged over all drop sizes
% within the distribution defined

% (3) g_avg - asymmetry averaged over all drop sizes
% within the distribution defined




% By Andrew John Buggee
%%

function [ssa_avg, Qe_avg, g_avg] = average_mie_over_size_distribution(r_eff, dist_var, wavelength,...
    index_of_refraction, size_distribution)

% ---------------------------
% ----- CHECK INPUTS --------
% ---------------------------

% The length of r_eff defines the number of unique droplet distributions.
% The length of the distribution variance must be the same length as the
% effective radius


% For the remaining entries, make sure they all have the same length as the
% length of the effective radius vector
if length(r_eff)~=length(dist_var)
    error([newline, 'The first two inputs must be the same length',newline])
end


% Check the wavelength input. This can only be a single value, or a vector
% with three values
if length(wavelength)~=1 && length(wavelength)~=3
    error([newline, 'The wavelength input can only be a single value, implying a monochromatic calculation',...
        newline, ', or have 3 values: the starting value, the end value, and the sampling interval',newline])
end



% ==============================


if strcmp(size_distribution, 'gamma')==true


    % For some effective radius, we have defined a droplet size distribution
    % ----- for a gamma distribution -----

    ssa_avg = zeros(1,length(r_eff));
    Qe_avg = zeros(1,length(r_eff));
    g_avg = zeros(1, length(r_eff));




    % Lets define the radius vector the spans the full size distribution
    % By doing this outside the loop, we only have to write and run a
    % single mie calculation
    r = linspace(min(r_eff)/100, max(r_eff)*8, 300);

    % for each effective radius we create a radius vector that spans the
    % full size distribution. We need to calculate the extinction
    % efficiency and the single scattering albedo for each value in our
    % radius vector


    % ---------------------------------------------------------------
    % ---------------     RUN MIE CALCULATIONS    -------------------
    % ---------------------------------------------------------------

    % define the wavelength
    % The wavelength input is defined as follows:
    % [wavelength_start, wavelength_end, wavelength_step].
    if length(wavelength)==1
        % monochromatic calculation
        wavelength = [wavelength, wavelength, 0];          % nanometers
        % define the wavelength vector
        wl_vector = wavelength(1);               % nanometers
    elseif length(wavelength)==3
        % broadband calculation
        wavelength = [wavelength(1), wavelength(2), wavelength(3)];          % nanometers
        % define the wavelength vector
        wl_vector = wavelength(1):wavelength(3):wavelength(2);               % nanometers
    end


    % The first entry belows describes the type of droplet distribution
    % that should be used. The second describes the distribution width. If
    % running a mono-dispersed calculation, no entry for distribution width is
    % required.
    size_distribution = {'mono', []};           % droplet distribution

    % What mie code should we use to compute the scattering properties?
    mie_program = 'MIEV0';               % type of mie algorithm to run

    % Do you want a long or short error file?
    err_msg_str = 'verbose';


    % Define the size of the scatterer and its scattering properties
    % Assuming a pure homogenous medium composed of a single substance.
    % The radius input is defined as [r_start, r_end, r_step].
    % where r_step is the interval between radii values (used only for
    % vectors of radii). A 0 tells the code there is no step. Finally, the
    % radius values have to be in increasing order.
    mie_radius = [r(1), r(end), r(2) - r(1)];    % microns




    % Create a mie file
    [input_filename, output_filename, mie_folder] = write_mie_file(mie_program, index_of_refraction,...
        mie_radius,wavelength,size_distribution, err_msg_str);

    % run the mie file
    [~] = runMIE(mie_folder,input_filename,output_filename);

    % Read the output of the mie file
    [ds,~,~] = readMIE(mie_folder,output_filename);

    % -----------------------------------------------------------------



    % Step through each modal radius and compute the average over the size
    % distribution
    for rr = 1:length(r_eff)

        mu = dist_var(rr)+3;                                            % to ensure I have the correct gamma distribution



        b = mu/r_eff(rr);                                   % exponent parameter
        N = mu^(mu+1)/(gamma(mu+1) * r_eff(rr)^(mu+1));     % normalization constant

        n_r = N*r.^mu .* exp(-b*r);                            % gamma droplet distribution



        % step through each wavelength. At each wavelength we integrate over
        % the vector r, a range of droplet sizes

        for ww = 1:length(wl_vector)

            % Average single scattering albedo over a droplet distribution

            ssa_avg(ww,rr) = trapz(r, ds.Qsca(ww,:) .* n_r)./...
                trapz(r, ds.Qext(ww,:) .* n_r);

            % Compute the average asymmetry parameter over a size distribution


            g_avg(ww,rr) = trapz(r, ds.asymParam(ww,:) .* n_r)./...
                trapz(r, n_r);

            % Compute the average extinction efficiency over a droplet size
            % distribution

            Qe_avg(ww,rr) = trapz(r, ds.Qext(ww,:) .* n_r)./...
                trapz(r, n_r);


        end

    end





else

    error([newline, 'I can only integrate over gamma distributions at the moment...',newline])

end






end