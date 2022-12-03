%% Testing 2-stream monte carlo

% By Andrew John Buggee

%% Define inputs

clear variables


% Define the boundaries of the medium
inputs.tau_y_lower_limit = 0;
inputs.tau_y_upper_limit = 2;

% define the solar zenith angle
% This is the angle of the incident radiation with respect to the medium
% normal direction
inputs.solar_zenith_angle = 0;                  % deg from zenith

% Define the albedo of the bottom boundary (tau upper limit)
inputs.albedo_maxTau = 0;


% Layers can differ by radius, by material, or both. If the layers differ
% by radius, create a vector describing each layer radius from top to
% bottom (tau = 0 to tau = tau0). If the layers differ by index of
% refraction, create a vector describing each layer indec of refraction
% from top to bottom. If both are true, create a cell array for both the
% changing radii and the changing index of refraction

inputs.layerRadii = 10;      % radius of spheres in each layer

% There HAS to be even spacing for the mie calcualtion to take place!


% Define the number of layers within the medium that differ
inputs.N_layers = length(inputs.layerRadii);

% Define the layer boundaries given the number of layers and the boundaries
% of the entire medium
inputs.layerBoundaries = linspace(inputs.tau_y_lower_limit, inputs.tau_y_upper_limit, inputs.N_layers +1);


% Define the number of photons to inject into the medium
inputs.N_photons = 1e5;

%% MIE CALCULATIONS

% ------------------------------------------------------------------
% Run Mie calculation to obtain single scatering albedo and asymmetry
% parameter
% ------------------------------------------------------------------

% define the wavelength
% The wavelength input is defined as follows:
% [wavelength_start, wavelength_end, wavelength_step].
inputs.mie.wavelength = [550, 550, 0];          % nanometers

% The first entry belows describes the type of droplet distribution
% that should be used. The second describes the distribution width. If
% running a mono-dispersed calculation, no entry for distribution width is
% required.
inputs.mie.distribution = {'mono', []};           % droplet distribution

% What mie code should we use to compute the scattering properties?
inputs.mie.mie_program = 'MIEV0';               % type of mie algorithm to run

% Do you want a long or short error file?
inputs.mie.err_msg_str = 'verbose';

% What is the index of refraction of the scatterer? If there is more than
% one scatterer, enter multiple values for the index of refraction
inputs.mie.indexOfRefraction = 'water';
%inputs.mie.indexOfRefraction = 1.33 + 9.32e-5i;

% Define the size of the scatterer and its scattering properties
% Assuming a pure homogenous medium composed of a single substance.
% The radius input is defined as [r_start, r_end, r_step].
% where r_step is the interval between radii values (used only for
% vectors of radii). A 0 tells the code there is no step. Finally, the
% radius values have to be in increasing order.
if inputs.N_layers==1
    inputs.mie.radius = [inputs.layerRadii, inputs.layerRadii, 0];    % microns
else
    % define the min, max, and step. Record the vector because these are
    % the exact values used in the mie calculations

    % Define the radius vector for our mie calculation
    inputs.mie.radiusVector = sort(inputs.layerRadii);
    % min value
    inputs.mie.radius(1) = min(inputs.layerRadii);
    % max value
    inputs.mie.radius(2) = max(inputs.layerRadii);
    % step value
    inputs.mie.radius(3) = inputs.mie.radiusVector(2) - inputs.mie.radiusVector(1);    % microns

end



% Create a mie file
[input_filename, output_filename, mie_folder] = write_mie_file(inputs.mie.mie_program, inputs.mie.indexOfRefraction,...
    inputs.mie.radius,inputs.mie.wavelength,inputs.mie.distribution, inputs.mie.err_msg_str);

% run the mie file
[~] = runMIE(mie_folder,input_filename,output_filename);

% Read the output of the mie file
[ds,~,~] = readMIE(mie_folder,output_filename);

% --------------------------------------------------
% Outputs vary by wavelength along the row dimension
% Outputs vary by radii along the column dimension
% --------------------------------------------------

% ------------------------------------------------------------
% Organize the scattering properties from medium top to bottom
% Use the definition of the each layer's radii
% This vector is define from layer top to layer bottom
% ------------------------------------------------------------
if inputs.N_layers>1
inputs.g = zeros(1, length(inputs.layerRadii));
inputs.ssa = zeros(1, length(inputs.layerRadii));
inputs.Qe = zeros(1, length(inputs.layerRadii));
for LL = 1:inputs.N_layers

    % define the asymmetry parameter
    [~,index] = min(abs((inputs.mie.radius(1):inputs.mie.radius(3):inputs.mie.radius(2))-inputs.layerRadii(LL)));
    inputs.g(:,LL) = ds.asymParam(:,index);

    % define the single scattering albedo
    inputs.ssa(:,LL) = ds.ssa(:,index);

    % define the extinction efficiency
    inputs.Qe(:,LL) = ds.Qext(:,index);

end

else
    inputs.g = ds.asymParam;

    % define the single scattering albedo
    inputs.ssa = ds.ssa;

    % define the extinction efficiency
    inputs.Qe = ds.Qext;
end


% We don't need the rest of the mie computations, so lets delete them
clear ds




%% Do you want to integrate over a size distribution?

% --------------------------------------------------
inputs.mie.integrate_over_size_distribution = false;
% --------------------------------------------------



if inputs.mie.integrate_over_size_distribution==true

    % Define the type of size distribution
    size_dist = 'gamma';

    % Define the distribution variance, depending on the distribution type used
    % Has to be the same length as the numer of layers in our medium
    dist_var = linspace(7,7,length(inputs.ssa));           % Typically value for liquid water clouds

    % Compute the average value for the single scattering albedo over a size
    % distribution
    [inputs.ssa_avg, inputs.Qe_avg, inputs.g_avg] = average_mie_over_size_distribution(inputs.ssa, inputs.g, inputs.Qe,...
        inputs.layerRadii,dist_var, inputs.mie.wavelength(1), inputs.mie.indexOfRefraction, size_dist);


end


%% Run 2 stream 2D monte carlo code
tic
[F_norm, final_state, photon_tracking, inputs] = twoStream_2D_monteCarlo(inputs);
toc


%% Compare the monte carlo solution with the analytical solution

plot_2strm_2D_monteCarlo(inputs,F_norm);


%% Do you want to save you results?

% save in the following folder
inputs.folder_name_2save = 'Monte_Carlo_Simulation_Results';
cd(inputs.folder_name_2save)
save(['2D_MC_',char(datetime('today')),'_Wavelength_',num2str(inputs.mie.wavelength(1)),...
    '_N-Photons_',num2str(inputs.N_photons),'_N-Layers_',num2str(inputs.N_layers),...
    '_Tau0_',num2str(inputs.tau_y_upper_limit),'_SZA_',num2str(inputs.solar_zenith_angle),'.mat'],...
    "inputs","F_norm", "final_state", "photon_tracking");
cd ..


%% Lets plot the max depth reached by each photon normalized by the total number of photons

%plot_probability_maxDepth(inputs, photon_tracking, 'probability')


%% Let's plot the conditional probability of an absorbed photon reach a max depth of tau

%plot_probability_absorbedDepth(inputs, final_state, photon_tracking, 'probability')

%% Let's plot the conditional probability of a photon that scattered out the cloud top
% reaching a max depth of tau

%plot_probability_scatterOutTop_maxDepth(inputs, final_state, photon_tracking, 'probability')


%% Let's plot two conditional probabilities on the same plot
% Plot the probability of a photon reaching a max depth of tau if it was
% scattered out the cloud top, and plot the probability of a photon being
% absorbed at a depth tau given that it was absorbed.

plot_probability_absANDscatTop_maxDepth(inputs, final_state, photon_tracking, 'pdf')


%% Plot a bar chart showing the probability of each final state

plot_probability_finalStates(final_state,inputs)

