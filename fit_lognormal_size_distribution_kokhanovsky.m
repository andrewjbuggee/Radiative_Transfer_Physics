% Fit a lognormal particle size distribution to a droplet distribution data set

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

%   (2) std_dev - the standard deviation of logarithmic values

%   (3) r_modal - the modal radius of the input data


% By Andrew John Buggee

%%

function [nr_fit, std_dev, min_rms] = fit_lognormal_size_distribution_kokhanovsky(nr_data, r)

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

std_dev_fit = logspace(-2,0,100);            % values of mu to test

N = 1./(sqrt(2*pi) * std_dev_fit);  % normalization constant


% Compute the fit for the first 100 values of mu
nr_fits = zeros(numel(r), numel(std_dev_fit));

rms_diff = zeros(1, numel(std_dev_fit));

for nn = 1:numel(std_dev_fit)
    
    nr_fits(:, nn) = N0 * N(nn) * 1./r .* exp(-log(r./r_modal).^2 ./(2*std_dev_fit(nn)^2));       % lognormal droplet distribution

    % Find the least squares difference between the fit true data

    rms_diff(nn) = sqrt(sum((nr_data - nr_fits(:,nn)).^2, 1));

end



% find the mu value associated witht he minimum rms difference
[min_rms, min_idx] = min(rms_diff);

% Grab the mu value and the nr_fit associated with the minimum rms 
std_dev = std_dev_fit(min_idx);
nr_fit = nr_fits(:, min_idx);




end