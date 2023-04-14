%% Integrate mie calculation for a single radii over a size distribution


% By Andrew John Buggee
%%

function [ssa_avg, Qe_avg, g_avg] = average_mie_over_size_distribution(ssa, g, Qe, r_eff, dist_var, wavelength, index_of_refraction, size_distribution)

% ---------------------------
% ----- CHECK INPUTS --------
% ---------------------------

% Each of the 4 inputs must be the same length
N(1) = length(ssa);
N(2) = length(g);
N(3) = length(Qe);
N(4) = length(r_eff);
N(5) = length(dist_var);

if all(N==N(1))~=true
    error([newline, 'The first four inputs must be the same length',newline])
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
    wavelength = [wavelength, wavelength, 0];          % nanometers

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

    for ii = 1:length(r_eff)

        mu = dist_var(ii)+3;                                            % to ensure I have the correct gamma distribution



        b = mu/r_eff(ii);                                   % exponent parameter
        N = mu^(mu+1)/(gamma(mu+1) * r_eff(ii)^(mu+1));     % normalization constant

        n_r = N*r.^mu .* exp(-b*r);                            % gamma droplet distribution



        % Average single scattering albedo over a droplet distribution

        ssa_avg(ii) = trapz(r, ds.Qsca .* n_r)./...
                       trapz(r, ds.Qext .* n_r);

        % Compute the average asymmetry parameter over a size distribution
    
        
        g_avg(ii) = trapz(r, ds.asymParam .* n_r)./...
                    trapz(r, n_r);
        
        % Compute the average extinction efficiency over a droplet size
        % distribution

        Qe_avg(ii) = trapz(r, ds.Qext .* n_r)./...
                    trapz(r, n_r);


    end





else

    error([newline, 'I can only integrate over gamma distributions at the moment...',newline])

end






end