% ----- Planks Function -----

% this function takes in a vector of wavelength values, and a vector of
% temperature values and computes the plank function, which defines the
% amount of radiance

% ----- inputs -----

% (1) - wavelenghts - a row vector defining the wavelengths to calculate
% the radiance for. units are defined by the units input

% (2) - temperature - a row vector. For each temperature, we create a
% row vector of radiance across the wavelengths defined by the user.
% units are in kelvin

% (3) - wavelength_units - a string defining the units to use for both the
% input and the output. The options are:
%       (a) 'microns'
%       (b) 'nanometers'

% ----- outputs -----

% (1) radiance - units are in watts/m^2/nm/sr


% By Andrew J. Buggee
%%

function [radiance] = planks_function(wavelengths,temperatures,wavelength_units)

% define constants

h = 6.626e-34; % Joules/sec - planks constant
c = 299792458; % meters/s - speed of light
k = 1.3806e-23; % J/Kelvin - Boltzmann constant

radiance = zeros(length(temperatures),length(wavelengths));

if strcmp(wavelength_units,'microns')
    
    h = h * (1e6)^2; % converts h to microns
    c = c * (1e6); % converts c to microns per secon
    k = k * (1e6)^2; % converts k to microns
    
    for ii = 1:length(temperatures)
        
        radiance(ii,:) = 2*h*c^2./(wavelengths.^5 .*(exp(h*c./(wavelengths*k*temperatures(ii))) - 1));
        
    end
    
    
elseif strcmp(wavelength_units,'nanometers')
    
    h = h * (1e9)^2; % converts h to microns
    c = c * (1e9); % converts c to microns per secon
    k = k * (1e9)^2; % converts k to microns
    
    for ii = 1:length(temperatures)
        
        radiance(ii,:) = 2*h*c^2./(wavelengths.^5 .*(exp(h*c./(wavelengths*k*temperatures(ii))) - 1));
        
    end
    
end




end

