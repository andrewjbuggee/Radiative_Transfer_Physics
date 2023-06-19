% Fit a gamma particle size distribution to a droplet distribution data set

% INPUTS:
%   (1) nr_data - this is a vector of the droplet distribution data. The
%   units must by #/cm^(-3)/micron. It must be the same length as the r
%   vector

%   (2) r - this is the independent vector associated with the distribution
%   data nr_data. It must be in units of microns.

%   (3) 

% OUTPUTS:
%   (1) n(r) - the number concentration of droplets for a given radius r -
%   this output is a vector of the same length as r, which is a hard-coded
%   vector

%   (2) r - the independent variable to defines our distribution n(r). 


% By Andrew John Buggee

%%

function [nr_fit, mu, var_coef, min_rms] = fit_gamma_size_distribution_kokhanovsky(nr_data, r)

% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 3 inputs, droplet radius, cloud optical
% depth, and the altitude vector associated with this cloud


if nargin~=2
    error([newline,'Not enough inputs. Need 3: droplet effective radius, total number concentration, and the alpha parameter value.', newline])
end

% Check to make sure r is the same length as nr_data

if length(r)~=length(nr_data)

    error([newline,'The inputs r and nr_data must be the same length', newline])           % says who?
end



%% using the independent vector provide, compute the gamma distribution for many alpha and r0 values
% This definition follows the definition of A. Kokhonovsky in 'Cloud
% Optics' on page 2

% make sure r and nr_data are column vectors
r = reshape(r, [],1);
nr_data = reshape(nr_data, [], 1);

% lets find the modal radius
[~,idx] = max(nr_data);
r_modal = r(idx);

% Compute the total number of droplets
N0 = trapz(r, nr_data);

mu_fit = 1:0.1:50;            % values of mu to test


N = mu_fit.^(mu_fit+1)./(gamma(mu_fit+1) .* r_modal.^(mu_fit+1));  % normalization constant


% Compute the fit for the first 100 values of mu
nr_fits = zeros(numel(r), numel(mu_fit));

rms_diff = zeros(1, numel(mu_fit));

for nn = 1:numel(mu_fit)
    nr_fits(:, nn) = N0 * N(nn) .* r.^mu_fit(nn) .* exp(-mu_fit(nn)*r./r_modal);                            % gamma droplet distribution

    % Find the least squares difference between the fit true data

    rms_diff(nn) = sqrt(sum((nr_data - nr_fits(:,nn)).^2, 1));

end



% find the mu value associated witht he minimum rms difference
[min_rms, min_idx] = min(rms_diff);

% Grab the mu value and the nr_fit associated with the minimum rms 
mu = mu_fit(min_idx);
var_coef = 1/sqrt(1+mu);
nr_fit = nr_fits(:, min_idx);




end