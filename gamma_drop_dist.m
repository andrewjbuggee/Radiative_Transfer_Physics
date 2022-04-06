%% This function computes a Gamma droplet distribution

% INPUTS:
%   (1) re - effective droplet radius (microns) - this is a single value
%   that defines the effective droplet radius of our size distribution.
%   This is the third moment of the distribution to the second moment.
%   Think of it like an average

%   (2) N0 - total droplet concentration (cm^(-3)) - this is the total
%   droplet number concentration including all sizes. 

%   (3) z - altitude above sea level (kilometers) - this is a vector with
%   the same length as re. The first value defines the base of
%   the cloud. If there are multiple values, each entry defines the start
%   of the next layer where the re value changes.

%   (4) H - geometric thickness of the cloud (kilometers) - this is a
%   single value that defines the total geometric cloud thickness.

%   (5) lambda - wavelength that defines the cloud optical depth
%   (nanometers) - This is the wavelength that defines the cloud optical
%   depth. There should only be one value!

%   (6) distribution_str - a string telling the code with droplet size 
%   distribution to use  - One can chose from two options:
%       (a) 'mono' - monodispersed distribution
%       (b) 'gamma' - gamma droplet distribution. By default this will use
%       a gamma distribution with an alpha value of 7, which is typical for
%       liquid water clouds.

% OUTPUTS:
%   (1) .Dat file saved in the libRadTran folder:
%   /.../libRadtran-2.0.4/data/wc

% All look up tables were computed using the opensoure radiative transfer
% code, LibRadTran

% By Andrew John Buggee
%%
function [outputArg1,outputArg2] = untitled4(inputArg1,inputArg2)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

