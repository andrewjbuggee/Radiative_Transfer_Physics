%% This function computes a Gamma droplet distribution

% INPUTS:
%   (1) re - effective droplet radius (microns) - this is a single value
%   that defines the effective droplet radius of our size distribution.
%   This is the third moment of the distribution to the second moment.
%   Think of it like an average

%   (2) N0 - total droplet concentration (cm^(-3)) - this is the total
%   droplet number concentration including all sizes. 

%   (3) alpha - paramter that is at least 1, and can be as large as you
%   want. The defines the spread of our distribution

% OUTPUTS:
%   (1) n(r) - the number concentration of droplets for a given radius r -
%   this output is a vector of the same length as r, which is a hard-coded
%   vector


% By Andrew John Buggee
%%

function [nr, r] = gamma_drop_dist(re,N0, alpha)


% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 3 inputs, droplet radius, cloud optical
% depth, and the altitude vector associated with this cloud


if nargin~=3
    error([newline,'Not enough inputs. Need 3: droplet effective radius, total number concentration, and the alpha parameter value.', newline])
end

% Check to make sure re is the same length as the altitude vector

if alpha<1

    error([newline,'The alpha parameter cannot be less than 1', newline])
end

if N0<0
    
    error([newline,'The number concentration must be greater than 0!', newline])
end

if re<0
    
    error([newline,'The effective droplet radius must be greater than 0!',newline])
end


%%


% ------------------------------------------------------------
% --------------- Create Droplet Distribution ----------------
% ------------------------------------------------------------





% create the r vector based on the input value for the effective droplet
% radius

r = linspace(re/10, re*10, 100);      % microns - r vector that spans our function

% ------------------------------------------------------------
% --------------- THIS COMPUTATION IS INCORRECT!!!!!! --------
% ------------------------------------------------------------
nr = N0*r.^alpha .* exp(-(alpha+3)*r./re);       % cm^-3 - number concentration at a given r value




end

