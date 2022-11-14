%% Testing 2-stream monte carlo

% By Andrew John Buggee

%% Define inputs

clear variables


% Define the boundaries of the medium
inputs.tau_lower_limit = 0;
inputs.tau_upper_limit = 10;

% Define the albedo of the bottom boundary (tau upper limit)
inputs.albedo_maxTau = 0.05;

% Define the number of photons to inject into the medium
inputs.N_photons = 10000;

% define the wavelength
inputs.wavelength = 500;           % nanometers
inputs.wavelength = 1950;          % nanometers

% Define the size of the scatterer and its scattering properties
% Assuming a pure homogenous medium composed of a single substance

inputs.radius = 10;                     % microns

% ------------------------------------------------------------------
% Run Mie calculation to obtain single scatering albedo and asymmetry
% parameter
% ------------------------------------------------------------------

inputs.mie.distribution = 'mono';           % droplet distribution
inputs.mie.distribution_width = 0;
inputs.mie.mie_program = 'MIEV0';               % type of mie algorithm to run
inputs.mie.err_msg_str = 'verbose';
inputs.mie.indexOfRefraction = 'water';
%inputs.mie.indexOfRefraction = 1.33 + 9.32e-5i;




% Create a mie file
[input_filename, output_filename, mie_folder] = write_mie_file(inputs.mie.mie_program, inputs.mie.indexOfRefraction,inputs.radius,...
    inputs.wavelength,inputs.mie.distribution, inputs.mie.distribution_width, inputs.mie.err_msg_str);

% run the mie file
runMIE(mie_folder,input_filename,output_filename);

% Read the output of the mie file
[ds,~,~] = readMIE(mie_folder,output_filename);

% Define the single scattering albedo 
inputs.ssa = ds.ssa;

% Define the asymmetry parameter
inputs.g = ds.asymParam;


%% Run 2 stream monte carlo code

[F_norm, final_state, photon_tracking] = twoStream_monteCarlo(inputs);


%% Compare the monte carlo solution with the analytical solution

plot_2strm_monteCarlo_with_analytical(inputs,F_norm);


%% Lets plot the max depth reached by each photon normalized by the total number of photons

plot_probability_maxDepth(inputs, photon_tracking, 'probability')


%% Let's plot the conditional probability of an absorbed photon reach a max depth of tau

%plot_probability_absorbedDepth(inputs, final_state, photon_tracking, 'probability')

%% Let's plot the conditional probability of a photon that scattered out the cloud top
% reaching a max depth of tau

%plot_probability_scatterOutTop_maxDepth(inputs, final_state, photon_tracking, 'probability')


%% Let's plot two conditional probabilities on the same plot
% Plot the probability of a photon reaching a max depth of tau if it was
% scattered out the cloud top, and plot the probability of a photon being
% absorbed at a depth tau given that it was absorbed.

plot_probability_absANDscatTop_maxDepth(inputs, final_state, photon_tracking, 'probability')


%% Plot a bar chart showing the probability of each final state

plot_probability_finalStates(final_state,inputs)

