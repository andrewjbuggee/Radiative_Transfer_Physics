% Create a gamma particle size distribution

% INPUTS:
%   (1) r_modal - modal radius (microns) - this is a single value
%   that defines the radius that would be most commonly observed if we
%   randomly sampled our gamma particle size distribution a large number of
%   times

%   (2) alpha - paramter that is at least 1, and can be as large as you
%   want. This defines the spread of our distribution. 7 is often used for
%   liquid water clouds

%   (2) N0 - total droplet concentration (cm^(-3)) - this is the total
%   droplet number concentration including all sizes. If n(r) is integrated
%   over r, we get N0.

% OUTPUTS:
%   (1) n(r) - the number concentration of droplets for a given radius r -
%   this output is a vector of the same length as r, which is a hard-coded
%   vector

%   (2) r - the independent variable to defines our distribution n(r). 


% By Andrew John Buggee

%%

function [n_r,r] = gamma_size_distribution(r_modal,alpha, N0)

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

if r_modal<0
    
    error([newline,'The effective droplet radius must be greater than 0!',newline])
end

%%

% For some modal radius, we have define a droplet size distribution over a
% range of radii values

r = linspace(0.02*r_modal, 8*r_modal, 100);                  % microns - vector based on C.Emde (2016)

mu = alpha+3;                                            % to ensure I have the correct gamma distribution


b = mu/r_modal;                                   % exponent parameter
N = mu^(mu+1)/(gamma(mu+1) * r_modal^(mu+1));  % normalization constant

n_r = N0 * N*r.^mu .* exp(-b*r);                            % gamma droplet distribution



end